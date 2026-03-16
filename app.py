import os
import base64
import shutil
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

from scoring_engine.v5_opc_scoring_engine import score_v5_opc
from utils.metadata_utils import merge_metadata

from analysis.validation_metrics import compute_validation_report
from analysis.purification import build_purification_outputs
from analysis.surface_panel import build_surface_panel_outputs

from analysis.visualization import (
    plot_umap_by_cluster,
    plot_umap_by_v5_score,
    plot_umap_highlight_top_clusters,
    ensure_umap,
)

from reporting.pdf_report import generate_pdf_report


OPCODE_VERSION = "1.1.0"
OUTPUT_DIR = "outputs"
MAX_CLUSTER_PLOT = 10

# Large-dataset safeguards
HUGE_DATASET_CELL_THRESHOLD = 500_000
HUGE_DATASET_SAMPLE_FRAC = 0.01
HUGE_DATASET_MIN_SAMPLE_N = 20_000
HUGE_DATASET_MAX_SAMPLE_N = 100_000
RANDOM_SEED = 42

# Figure export DPI for faster PDF assembly
FIG_DPI = 200

st.set_page_config(layout="wide")
st.title(f"OPCODE — OPC Detection Engine v{OPCODE_VERSION}")
st.markdown(
    """
Upload a **.h5ad dataset** or provide a **local file path** to compute OPCODE V5 scores,
identify OPC clusters, discover surface markers, and generate a scientific PDF report.

For very large datasets, using a **local path** is recommended instead of browser upload.
If no cluster metadata are present, OPCODE will generate **Leiden clusters automatically**.
For extremely large datasets, OPCODE uses a **sampled analysis mode** for interactive
visualization and PDF generation.
"""
)

os.makedirs(OUTPUT_DIR, exist_ok=True)

if "pdf_bytes" not in st.session_state:
    st.session_state["pdf_bytes"] = None
if "pdf_name" not in st.session_state:
    st.session_state["pdf_name"] = None
if "pdf_progress_active" not in st.session_state:
    st.session_state["pdf_progress_active"] = False
if "pdf_progress_bar" not in st.session_state:
    st.session_state["pdf_progress_bar"] = None
if "pdf_status_box" not in st.session_state:
    st.session_state["pdf_status_box"] = None

analysis_progress_bar = st.progress(0)
analysis_status_box = st.empty()


def update_analysis_progress(value: int, stage: str):
    analysis_progress_bar.progress(value)
    analysis_status_box.info(f"{stage}...")


def update_pdf_progress(value: int, stage: str):
    pdf_bar = st.session_state.get("pdf_progress_bar")
    pdf_box = st.session_state.get("pdf_status_box")

    if pdf_bar is not None:
        pdf_bar.progress(value)
    if pdf_box is not None:
        pdf_box.info(f"{stage}...")


# =====================================================
# LOADERS / HELPERS
# =====================================================

@st.cache_data(show_spinner=False)
def get_file_size_gb(path: str) -> float:
    return os.path.getsize(path) / (1024 ** 3)


def load_dataset(path: str):
    file_size_gb = get_file_size_gb(path)

    if file_size_gb > 2:
        return sc.read_h5ad(path, backed="r")

    return sc.read_h5ad(path)


def detect_cluster_column(adata):
    for col in [
        "cluster_alias",
        "clusters",
        "cluster",
        "leiden",
        "cell_cluster",
        "louvain",
        "seurat_clusters",
    ]:
        if col in adata.obs.columns:
            return col
    return None


def _has_usable_neighbors(adata) -> bool:
    return (
        "neighbors" in adata.uns
        and "connectivities" in getattr(adata, "obsp", {})
        and "distances" in getattr(adata, "obsp", {})
    )


def _looks_like_ensembl_ids(index_values) -> bool:
    vals = pd.Index(index_values).astype(str)
    if len(vals) == 0:
        return False
    sample = vals[: min(200, len(vals))]
    score = np.mean(
        [
            x.startswith("ENSMUSG")
            or x.startswith("ENSG")
            or x.startswith("ENSMUST")
            or x.startswith("ENSMUSP")
            for x in sample
        ]
    )
    return score >= 0.5


def _matrix_max_value(adata):
    try:
        x_max = adata.X.max()
        if hasattr(x_max, "item"):
            x_max = x_max.item()
        return float(x_max)
    except Exception:
        return None


def _looks_logged(adata) -> bool:
    x_max = _matrix_max_value(adata)
    if x_max is None:
        return True
    return x_max < 50


def normalize_for_opcode(adata):
    """
    Normalize expression before OPCODE scoring so score magnitudes are more
    comparable across datasets.

    If the matrix does not appear log-transformed:
      1) normalize total counts per cell
      2) log1p transform
    Otherwise leave it unchanged.
    """
    adata = adata.copy()

    if not _looks_logged(adata):
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    return adata


def normalize_gene_names(adata):
    candidate_cols = [
        "gene_name",
        "gene_names",
        "gene_symbol",
        "gene_symbols",
        "symbol",
        "symbols",
        "feature_name",
        "feature_names",
        "mgi_symbol",
        "gene_short_name",
    ]

    if not _looks_like_ensembl_ids(adata.var_names):
        return adata, None, False

    symbol_col = None
    for col in candidate_cols:
        if col in adata.var.columns:
            values = adata.var[col].astype(str)
            nonempty = values.str.strip().replace({"nan": "", "None": ""})
            if (nonempty != "").sum() > 0:
                symbol_col = col
                break

    if symbol_col is None:
        return adata, None, False

    adata = adata.copy()
    adata.var["original_gene_id"] = adata.var_names.astype(str)

    symbols = adata.var[symbol_col].astype(str).fillna("").str.strip()
    cleaned = []
    original_ids = adata.var["original_gene_id"].astype(str).tolist()

    for sym, gid in zip(symbols.tolist(), original_ids):
        if sym == "" or sym.lower() == "nan" or sym.lower() == "none":
            cleaned.append(gid)
        else:
            cleaned.append(sym)

    adata.var_names = pd.Index(cleaned)
    adata.var_names_make_unique()

    return adata, symbol_col, True


def ensure_cluster_column(adata):
    cluster_column = detect_cluster_column(adata)
    if cluster_column is not None:
        return adata, cluster_column, False

    adata = ensure_umap(adata)

    if "X_pca" not in adata.obsm:
        sc.pp.pca(adata, n_comps=min(50, max(2, adata.shape[1] - 1)))

    if not _has_usable_neighbors(adata):
        sc.pp.neighbors(
            adata,
            n_neighbors=15,
            n_pcs=min(30, adata.obsm["X_pca"].shape[1]),
        )

    if "leiden" not in adata.obs.columns:
        sc.tl.leiden(adata, resolution=1.0)

    return adata, "leiden", True


def safe_read_highconf_ids(output_dir: str):
    for path in [
        os.path.join(output_dir, "opc_highconf.csv"),
        os.path.join(output_dir, "opc_high_confidence.csv"),
    ]:
        if os.path.exists(path):
            df = pd.read_csv(path)

            for col in ["cell_id", "Cell_ID", "index", "barcode"]:
                if col in df.columns:
                    return df[col].astype(str).tolist()

            if len(df.columns) > 0:
                return df.iloc[:, 0].astype(str).tolist()

    return None


def save_highconf_umap(adata, highconf_ids, dataset_name):
    if not highconf_ids:
        return

    if "X_umap" not in adata.obsm:
        adata = ensure_umap(adata)

    obs_names = pd.Index(adata.obs_names.astype(str))
    highlight = obs_names.isin(set(map(str, highconf_ids)))

    coords = adata.obsm["X_umap"]
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(
        coords[~highlight, 0],
        coords[~highlight, 1],
        s=8,
        c="lightgrey",
        linewidths=0,
    )
    ax.scatter(
        coords[highlight, 0],
        coords[highlight, 1],
        s=18,
        c="red",
        linewidths=0,
    )

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("")
    ax.legend(["Other", "High-Confidence OPC Cells"], loc="best", frameon=False)

    fig.savefig(
        os.path.join(OUTPUT_DIR, f"{dataset_name}_umap_highconf.png"),
        dpi=FIG_DPI,
        bbox_inches="tight",
    )
    plt.close(fig)


def _safe_remove_path(path: str):
    try:
        if os.path.isfile(path):
            os.remove(path)
        elif os.path.isdir(path):
            shutil.rmtree(path, ignore_errors=True)
    except Exception:
        pass


def cleanup_outputs_workspace():
    """
    Remove everything inside OUTPUT_DIR after the PDF has been loaded into memory.
    This keeps the pipeline file-based internally, but prevents clutter on disk.
    """
    if not os.path.exists(OUTPUT_DIR):
        return

    for name in os.listdir(OUTPUT_DIR):
        path = os.path.join(OUTPUT_DIR, name)
        _safe_remove_path(path)


def render_pdf_download():
    pdf_bytes = st.session_state.get("pdf_bytes")
    pdf_name = st.session_state.get("pdf_name")

    if pdf_bytes is None or pdf_name is None:
        return

    st.subheader("Download PDF Report")

    st.download_button(
        label="Download PDF Report",
        data=pdf_bytes,
        file_name=pdf_name,
        mime="application/pdf",
        key="pdf_download_button",
    )

    b64 = base64.b64encode(pdf_bytes).decode("utf-8")
    href = (
        f'<a href="data:application/pdf;base64,{b64}" '
        f'download="{pdf_name}" '
        f'style="display:inline-block;padding:0.5rem 0.75rem;'
        f'background:#1f77b4;color:white;text-decoration:none;'
        f'border-radius:0.5rem;margin-top:0.5rem;">'
        f'Direct Download Link'
        f'</a>'
    )
    st.markdown(href, unsafe_allow_html=True)


def _compute_sample_size(n_obs: int) -> int:
    pct_n = int(round(n_obs * HUGE_DATASET_SAMPLE_FRAC))
    pct_n = max(HUGE_DATASET_MIN_SAMPLE_N, pct_n)
    pct_n = min(HUGE_DATASET_MAX_SAMPLE_N, pct_n)
    pct_n = min(pct_n, n_obs)
    return pct_n


def maybe_materialize_or_sample(adata):
    original_n_obs = int(adata.n_obs)

    if not getattr(adata, "isbacked", False):
        return adata, False, original_n_obs, original_n_obs

    if original_n_obs > HUGE_DATASET_CELL_THRESHOLD:
        sample_n = _compute_sample_size(original_n_obs)

        update_analysis_progress(
            15,
            (
                f"Preparing sampled analysis mode for {original_n_obs:,} cells "
                f"at {HUGE_DATASET_SAMPLE_FRAC:.0%} with bounds "
                f"[{HUGE_DATASET_MIN_SAMPLE_N:,}, {HUGE_DATASET_MAX_SAMPLE_N:,}]"
            ),
        )

        rng = np.random.default_rng(RANDOM_SEED)
        sample_idx = np.sort(rng.choice(original_n_obs, size=sample_n, replace=False))

        update_analysis_progress(22, f"Loading sampled subset of {sample_n:,} cells into memory")
        adata = adata[sample_idx, :].to_memory()

        st.info(
            f"Sampled {sample_n:,} cells from {original_n_obs:,} total cells "
            "for scoring, visualization, and report generation."
        )
        return adata, True, original_n_obs, sample_n

    update_analysis_progress(18, "Loading large dataset into memory")
    adata = adata.to_memory()
    return adata, False, original_n_obs, original_n_obs


# =====================================================
# DATASET INPUT MODE
# =====================================================

st.header("Dataset Input")

input_mode = st.radio(
    "Choose dataset source",
    ["Upload file", "Use local file path"],
    horizontal=True,
)

uploaded_file = None
local_h5ad_path = None

if input_mode == "Upload file":
    uploaded_file = st.file_uploader("Upload .h5ad dataset", type=["h5ad"])
else:
    local_h5ad_path = st.text_input("Enter full local path to .h5ad file")

metadata_file = st.file_uploader("Optional metadata CSV", type=["csv"])

if input_mode == "Upload file":
    if uploaded_file is None:
        st.info("Waiting for dataset upload.")
        render_pdf_download()
        st.stop()

    dataset_name = os.path.splitext(uploaded_file.name)[0]
    temp_path = os.path.join(OUTPUT_DIR, uploaded_file.name)

    with open(temp_path, "wb") as f:
        f.write(uploaded_file.getbuffer())

    dataset_path = temp_path

else:
    if not local_h5ad_path:
        st.info("Enter a local .h5ad path to continue.")
        render_pdf_download()
        st.stop()

    if not os.path.exists(local_h5ad_path):
        st.error("The specified local file path does not exist.")
        render_pdf_download()
        st.stop()

    dataset_name = os.path.splitext(os.path.basename(local_h5ad_path))[0]
    dataset_path = local_h5ad_path


# =====================================================
# LOAD DATASET
# =====================================================

st.header("Dataset Loading")
try:
    update_analysis_progress(5, "Loading dataset")
    adata = load_dataset(dataset_path)
    st.write(adata)

    adata, sampled_mode, original_n_obs, sampled_n = maybe_materialize_or_sample(adata)
    update_analysis_progress(28, "Dataset ready for preprocessing")

except Exception as e:
    st.error(f"Dataset loading failed: {repr(e)}")
    st.exception(e)
    st.stop()


# =====================================================
# GENE NAME NORMALIZATION
# =====================================================

try:
    update_analysis_progress(32, "Normalizing gene identifiers")
    adata, gene_symbol_col, gene_names_normalized = normalize_gene_names(adata)

    if gene_names_normalized:
        st.success(
            f"Gene identifiers normalized using adata.var['{gene_symbol_col}'] "
            "so outputs prefer recognizable gene symbols."
        )
    else:
        st.info("Using existing gene identifiers as provided in the dataset.")

except Exception as e:
    st.warning(f"Gene-name normalization skipped: {e}")


# =====================================================
# METADATA
# =====================================================

if metadata_file is not None:
    try:
        update_analysis_progress(36, "Merging metadata")
        adata = merge_metadata(adata, metadata_file)
        st.success("Metadata merged successfully.")
    except Exception as e:
        st.warning(f"Metadata merge failed: {e}")


# =====================================================
# RUN OPCODE
# =====================================================

st.header("Running OPCODE V5")
try:
    update_analysis_progress(42, "Normalizing expression for OPCODE scoring")
    adata = normalize_for_opcode(adata)

    update_analysis_progress(45, "Running OPCODE V5 scoring")
    adata = score_v5_opc(adata)
    st.success("Scoring completed.")
    update_analysis_progress(55, "Scoring completed")
except Exception as e:
    st.error(f"Scoring failed: {e}")
    st.stop()

if "V5_OPC_score" not in adata.obs.columns:
    st.error("Scoring completed, but 'V5_OPC_score' was not found in adata.obs.")
    st.stop()

scores = adata.obs["V5_OPC_score"].astype(float).values


# =====================================================
# SCORE SUMMARY
# =====================================================

st.header("Score Summary")
col1, col2, col3 = st.columns(3)
col1.metric("Cells analyzed", len(scores))
col2.metric("Mean V5", round(float(np.mean(scores)), 3))
col3.metric("Max V5", round(float(np.max(scores)), 3))

if sampled_mode:
    st.caption(
        f"Sampled analysis mode: {sampled_n:,} cells analyzed from an original dataset of "
        f"{original_n_obs:,} cells."
    )

st.subheader("V5 Score Distribution")
fig_dist, ax = plt.subplots(figsize=(8, 5))
ax.hist(scores, bins=50)
ax.set_xlabel("V5 Score")
ax.set_ylabel("Cell Count")
ax.set_title("")
st.pyplot(fig_dist)
fig_dist.savefig(
    os.path.join(OUTPUT_DIR, f"{dataset_name}_score_distribution.png"),
    dpi=FIG_DPI,
    bbox_inches="tight",
)
plt.close(fig_dist)


# =====================================================
# CLUSTERS
# =====================================================

try:
    update_analysis_progress(62, "Preparing cluster structure")
    adata, cluster_column, generated_clusters = ensure_cluster_column(adata)

    if generated_clusters:
        st.info("No cluster metadata detected. OPCODE generated Leiden clusters automatically.")
    else:
        st.success(f"Using cluster column: {cluster_column}")

    update_analysis_progress(68, "Cluster structure ready")

except Exception as e:
    cluster_column = None
    st.warning(f"Cluster detection/generation failed: {e}")


# =====================================================
# VALIDATION METRICS
# =====================================================

st.header("Validation Metrics")
report = None
vis_score = None
try:
    report = compute_validation_report(
        adata,
        cluster_column=cluster_column or "leiden",
        class_column="class",
    )
    st.json(report)
    if isinstance(report, dict):
        vis_score = report.get("VIS", None)
except Exception as e:
    st.warning(f"Validation report unavailable: {e}")


# =====================================================
# CLUSTER RANKING
# =====================================================

if cluster_column is not None:
    st.subheader("Cluster Ranking")
    try:
        cluster_means = (
            adata.obs.groupby(cluster_column)["V5_OPC_score"]
            .mean()
            .sort_values(ascending=False)
        )

        n_clusters = len(cluster_means)
        if n_clusters > MAX_CLUSTER_PLOT:
            st.info(
                f"{n_clusters} clusters detected. Showing top {MAX_CLUSTER_PLOT} clusters by mean V5 score "
                f"in the plot for readability."
            )
            st.dataframe(cluster_means.head(MAX_CLUSTER_PLOT))
            cluster_plot = cluster_means.head(MAX_CLUSTER_PLOT)
        else:
            st.dataframe(cluster_means)
            cluster_plot = cluster_means

        fig_cluster, ax = plt.subplots(figsize=(8.5, 5.2))
        cluster_plot.plot(kind="bar", ax=ax)

        ax.set_xlabel("Cluster")
        ax.set_ylabel("Mean V5 Score")
        ax.set_title("")
        ax.set_xticklabels(cluster_plot.index.astype(str), rotation=45, ha="right")

        st.pyplot(fig_cluster)

        fig_cluster.savefig(
            os.path.join(OUTPUT_DIR, f"{dataset_name}_cluster_means.png"),
            dpi=FIG_DPI,
            bbox_inches="tight",
        )
        plt.close(fig_cluster)
    except Exception as e:
        st.warning(f"Cluster ranking failed: {e}")


# =====================================================
# UMAP
# =====================================================

st.header("UMAP Visualization")

try:
    update_analysis_progress(75, "Generating UMAP visualizations")
    adata = ensure_umap(adata)
except Exception as e:
    st.warning(f"UMAP computation failed: {e}")

if cluster_column is not None:
    try:
        fig_umap_cluster = plot_umap_by_cluster(adata, cluster_column)
        for axis in fig_umap_cluster.axes:
            axis.set_title("")
        st.pyplot(fig_umap_cluster)
        fig_umap_cluster.savefig(
            os.path.join(OUTPUT_DIR, f"{dataset_name}_umap_clusters.png"),
            dpi=FIG_DPI,
            bbox_inches="tight",
        )
        plt.close(fig_umap_cluster)
    except Exception as e:
        st.warning(f"Cluster UMAP failed: {e}")

try:
    fig_umap_v5 = plot_umap_by_v5_score(adata)
    for axis in fig_umap_v5.axes:
        axis.set_title("")
    st.pyplot(fig_umap_v5)
    fig_umap_v5.savefig(
        os.path.join(OUTPUT_DIR, f"{dataset_name}_umap_v5.png"),
        dpi=FIG_DPI,
        bbox_inches="tight",
    )
    plt.close(fig_umap_v5)
except Exception as e:
    st.warning(f"V5 score UMAP failed: {e}")

if cluster_column is not None:
    try:
        fig_umap_highlight = plot_umap_highlight_top_clusters(adata, cluster_column)
        for axis in fig_umap_highlight.axes:
            axis.set_title("")
        st.pyplot(fig_umap_highlight)
        fig_umap_highlight.savefig(
            os.path.join(OUTPUT_DIR, f"{dataset_name}_umap_opc_cluster.png"),
            dpi=FIG_DPI,
            bbox_inches="tight",
        )
        plt.close(fig_umap_highlight)
    except Exception as e:
        st.warning(f"Predicted OPC cluster UMAP failed: {e}")

update_analysis_progress(84, "UMAP visualizations completed")


# =====================================================
# PURIFICATION
# =====================================================

st.header("OPC Candidate Extraction")
highconf_ids = None

try:
    update_analysis_progress(89, "Extracting high-confidence OPC candidates")
    build_purification_outputs(
        adata=adata,
        output_dir=OUTPUT_DIR,
        cluster_col=cluster_column,
    )
    highconf_ids = safe_read_highconf_ids(OUTPUT_DIR)
    save_highconf_umap(adata, highconf_ids, dataset_name)
    st.success("High-confidence OPC cells exported.")
    update_analysis_progress(93, "High-confidence OPC candidates extracted")
except Exception as e:
    st.warning(f"Purification failed: {e}")


# =====================================================
# SURFACE MARKER DISCOVERY
# =====================================================

st.header("Surface Marker Discovery")
try:
    update_analysis_progress(96, "Running surface marker discovery")
    build_surface_panel_outputs(
        adata=adata,
        output_dir=OUTPUT_DIR,
        dataset_name=dataset_name,
        highconf_cell_ids=highconf_ids,
    )
    st.success("Surface marker ranking generated.")

    marker_file = os.path.join(OUTPUT_DIR, "surface_marker_rankings.csv")
    if os.path.exists(marker_file):
        marker_df = pd.read_csv(marker_file).head(10).copy()
        marker_df.index = marker_df.index + 1
        st.dataframe(marker_df)

    update_analysis_progress(98, "Surface marker discovery completed")
    update_analysis_progress(100, "Analysis completed")
except Exception as e:
    st.warning(f"Surface marker discovery failed: {e}")


# =====================================================
# PDF
# =====================================================

st.header("Generate Scientific Report")

pdf_progress_container = st.empty()
pdf_status_container = st.empty()

if st.session_state.get("pdf_progress_active"):
    st.session_state["pdf_progress_bar"] = pdf_progress_container.progress(0)
    st.session_state["pdf_status_box"] = pdf_status_container.empty()

if st.button("Generate PDF Report"):
    try:
        st.session_state["pdf_progress_active"] = True
        st.session_state["pdf_progress_bar"] = pdf_progress_container.progress(0)
        st.session_state["pdf_status_box"] = pdf_status_container.empty()
        st.session_state["pdf_status_box"].info("Starting PDF generation...")

        pdf_path = os.path.join(OUTPUT_DIR, f"{dataset_name}_V5_OPC_Report.pdf")

        generate_pdf_report(
            adata,
            "OPCODE",
            "Ansel Belani",
            "Delaware County Community College",
            "2026",
            dataset_name,
            pdf_path,
            OPCODE_VERSION,
            validation_metrics=report,
            vis_score=vis_score,
            progress_callback=update_pdf_progress,
        )

        with open(pdf_path, "rb") as f:
            pdf_bytes = f.read()

        st.session_state["pdf_bytes"] = pdf_bytes
        st.session_state["pdf_name"] = os.path.basename(pdf_path)

        cleanup_outputs_workspace()

        st.success("PDF generated successfully.")
        update_pdf_progress(100, "PDF report generated successfully")

        if sampled_mode:
            st.caption(
                f"This PDF report was generated from a sampled subset of {sampled_n:,} cells "
                f"drawn from {original_n_obs:,} total cells."
            )

    except Exception as e:
        st.session_state["pdf_bytes"] = None
        st.session_state["pdf_name"] = None
        st.error(f"PDF generation failed: {e}")

render_pdf_download()