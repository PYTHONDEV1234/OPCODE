import os
import tempfile
from typing import Any, Dict, List, Optional

import pandas as pd
import matplotlib.pyplot as plt
from xml.sax.saxutils import escape

from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Image,
    Table,
    TableStyle,
    PageBreak,
    KeepTogether,
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.utils import ImageReader

try:
    import yaml
except Exception:
    yaml = None


# =====================================================
# STYLES
# =====================================================

def _make_styles():
    styles = getSampleStyleSheet()

    styles.add(
        ParagraphStyle(
            name="CenterTitle",
            parent=styles["Heading1"],
            alignment=1,
            fontSize=16,
            leading=20,
            spaceAfter=10,
            textColor=colors.black,
        )
    )

    styles.add(
        ParagraphStyle(
            name="SectionTitle",
            parent=styles["Heading2"],
            alignment=1,
            fontSize=14,
            leading=18,
            spaceAfter=8,
            textColor=colors.black,
        )
    )

    styles.add(
        ParagraphStyle(
            name="Caption",
            parent=styles["Normal"],
            alignment=1,
            fontSize=9.5,
            leading=11.5,
            textColor=colors.black,
            spaceAfter=8,
        )
    )

    styles.add(
        ParagraphStyle(
            name="InterpretationText",
            parent=styles["Normal"],
            fontSize=10,
            leading=13,
            textColor=colors.black,
        )
    )

    styles.add(
        ParagraphStyle(
            name="Small",
            parent=styles["Normal"],
            fontSize=9,
            leading=11,
            textColor=colors.black,
        )
    )

    styles.add(
        ParagraphStyle(
            name="CellWrap",
            parent=styles["Normal"],
            fontSize=7.8,
            leading=9.2,
            textColor=colors.black,
        )
    )

    return styles


# =====================================================
# HELPERS
# =====================================================

def _safe_list(value):
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _safe_dict(value):
    if isinstance(value, dict):
        return value
    return {}


def _detect_cluster_column(adata):
    for col in ["cluster_alias", "clusters", "cluster", "leiden", "cell_cluster"]:
        if col in adata.obs.columns:
            return col
    return None


def _existing_path(candidates: List[str]) -> Optional[str]:
    for path in candidates:
        if path and os.path.exists(path):
            return path
    return None


def _paragraphize(text: str) -> str:
    lines = text.splitlines()
    out: List[str] = []

    for line in lines:
        stripped = line.strip()

        if stripped == "":
            out.append("<br/>")
            continue

        if stripped.startswith("•"):
            out.append(f"&bull; {escape(stripped[1:].strip())}<br/>")
            continue

        out.append(f"{escape(stripped)}<br/>")

    return "".join(out)


def _interpretation_box(text: str, styles) -> Table:
    title = Paragraph("<b>Interpretation</b>", styles["InterpretationText"])
    body = Paragraph(_paragraphize(text), styles["InterpretationText"])

    t = Table([[title], [body]], colWidths=[6.7 * inch])
    t.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#F3F7FB")),
                ("BOX", (0, 0), (-1, -1), 0.75, colors.HexColor("#9EB6CE")),
                ("INNERGRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#D5E2EF")),
                ("LEFTPADDING", (0, 0), (-1, -1), 8),
                ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                ("TOPPADDING", (0, 0), (-1, -1), 6),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ]
        )
    )
    return t


def _scaled_image(path: str, max_width: float = 6.5 * inch, max_height: float = 4.6 * inch) -> Image:
    img_reader = ImageReader(path)
    width, height = img_reader.getSize()

    if width <= 0 or height <= 0:
        return Image(path, width=max_width, height=max_height)

    scale = min(max_width / width, max_height / height)
    draw_width = width * scale
    draw_height = height * scale

    return Image(path, width=draw_width, height=draw_height)


def _summary_stats_table(adata, validation_metrics: Optional[Dict[str, Any]], vis_score: Any) -> Table:
    scores = adata.obs["V5_OPC_score"].astype(float)
    rows = [
        ["Metric", "Value"],
        ["Minimum V5 score", f"{scores.min():.3f}"],
        ["Maximum V5 score", f"{scores.max():.3f}"],
        ["Mean V5 score", f"{scores.mean():.3f}"],
        ["Median V5 score", f"{scores.median():.3f}"],
        ["Standard deviation", f"{scores.std():.3f}"],
    ]

    cluster_col = _detect_cluster_column(adata)
    if cluster_col is not None:
        cluster_means = (
            adata.obs.groupby(cluster_col)["V5_OPC_score"]
            .mean()
            .sort_values(ascending=False)
        )
        top_cluster = str(cluster_means.index[0])
        top_cluster_mean = float(cluster_means.iloc[0])
        top_cluster_n = int((adata.obs[cluster_col].astype(str) == top_cluster).sum())

        rows.extend(
            [
                ["Predicted OPC cluster", top_cluster],
                ["Top cluster mean V5", f"{top_cluster_mean:.3f}"],
                ["Cells in predicted OPC cluster", str(top_cluster_n)],
            ]
        )

    if isinstance(validation_metrics, dict):
        signal_metrics = _safe_dict(validation_metrics.get("signal_metrics"))
        if "Dynamic_Range" in signal_metrics:
            rows.append(["Dynamic range", f"{float(signal_metrics['Dynamic_Range']):.3f}"])
        if "Top_Cluster_Gap" in signal_metrics:
            rows.append(["Top cluster gap", f"{float(signal_metrics['Top_Cluster_Gap']):.3f}"])
        if "Positive_Fraction" in signal_metrics:
            rows.append(["Positive fraction", f"{float(signal_metrics['Positive_Fraction']):.3f}"])

    if vis_score is not None:
        try:
            rows.append(["VIS", f"{float(vis_score):.3f}"])
        except Exception:
            rows.append(["VIS", str(vis_score)])

    table = Table(rows, colWidths=[3.4 * inch, 1.8 * inch])
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#5B7C99")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("GRID", (0, 0), (-1, -1), 0.5, colors.black),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ]
        )
    )
    return table


def _marker_discovery_table(marker_csv: str) -> Optional[Table]:
    if not os.path.exists(marker_csv):
        return None

    df = pd.read_csv(marker_csv).copy()
    if df.empty:
        return None

    keep = [
        "gene",
        "mean_highconf",
        "mean_rest",
        "detect_highconf",
        "detect_rest",
        "log2fc_highconf_vs_rest",
    ]
    keep = [c for c in keep if c in df.columns]
    if not keep:
        return None

    df = df[keep].head(10).copy()

    for col in df.columns:
        if col != "gene":
            df[col] = df[col].map(lambda x: f"{float(x):.3f}" if pd.notnull(x) else "NA")

    headers = [
        "Gene" if c == "gene" else
        "Mean HighConf" if c == "mean_highconf" else
        "Mean Rest" if c == "mean_rest" else
        "Detect HighConf" if c == "detect_highconf" else
        "Detect Rest" if c == "detect_rest" else
        "log2FC" if c == "log2fc_highconf_vs_rest" else c
        for c in df.columns
    ]

    data = [headers] + df.values.tolist()
    col_widths = [1.2 * inch if c == "gene" else 0.95 * inch for c in df.columns]

    table = Table(data, colWidths=col_widths, repeatRows=1)
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#D9E6F2")),
                ("GRID", (0, 0), (-1, -1), 0.45, colors.black),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), 8),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    return table


def _make_top_marker_chart(marker_csv: str) -> Optional[str]:
    if not os.path.exists(marker_csv):
        return None

    df = pd.read_csv(marker_csv).copy()
    if df.empty or "gene" not in df.columns:
        return None

    score_col = None
    for candidate in [
        "ranking_score",
        "log2fc_highconf_vs_rest",
        "specificity_detect_delta",
        "mean_highconf",
    ]:
        if candidate in df.columns:
            score_col = candidate
            break

    if score_col is None:
        return None

    df = df.sort_values(score_col, ascending=False).head(10).reset_index(drop=True)
    df.index = df.index + 1

    labels = [f"{idx}. {gene}" for idx, gene in zip(df.index, df["gene"])]
    values = df[score_col].astype(float).tolist()

    fig, ax = plt.subplots(figsize=(7.0, 4.5))
    ax.barh(labels[::-1], values[::-1])
    ax.set_xlabel(score_col.replace("_", " ").title())
    ax.set_title("")
    plt.tight_layout()

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
    fig.savefig(tmp.name, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return tmp.name


def _parse_sorting_panel(panel_path: str) -> Dict[str, Any]:
    if not os.path.exists(panel_path) or yaml is None:
        return {}

    with open(panel_path, "r", encoding="utf-8") as f:
        return _safe_dict(yaml.safe_load(f))


def _simple_list_table(title: str, items: List[str], styles) -> Optional[Table]:
    items = [str(x) for x in _safe_list(items) if str(x).strip()]
    if not items:
        return None

    data = [[Paragraph(f"<b>{escape(title)}</b>", styles["CellWrap"])]]
    for item in items:
        data.append([Paragraph(escape(item), styles["CellWrap"])])

    t = Table(data, colWidths=[6.3 * inch], repeatRows=1)
    t.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#D9E6F2")),
                ("GRID", (0, 0), (-1, -1), 0.4, colors.black),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 2),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
            ]
        )
    )
    return t


def _clean_gating_steps(raw_steps: List[str]) -> List[str]:
    cleaned = []

    fallback_text = (
        "Positive gate: No surface markers passed evidence thresholds in this dataset. "
        "Consider testing established OPC anchors such as PDGFRα or NG2 experimentally."
    )

    for step in raw_steps:
        step_clean = str(step).strip()
        if not step_clean:
            continue

        if (
            "Try anchors_to_test experimentally" in step_clean
            or "expand surface DB" in step_clean
            or "No surface markers passed evidence thresholds" in step_clean
        ):
            step_clean = fallback_text

        while "Positive gate: Positive gate:" in step_clean:
            step_clean = step_clean.replace("Positive gate: Positive gate:", "Positive gate:")

        cleaned.append(step_clean)

    return cleaned


def _add_page_with_image(story, styles, title, interpretation, image_candidates, caption=None):
    content = [Paragraph(f"<b>{escape(title)}</b>", styles["CenterTitle"])]

    image_path = _existing_path(image_candidates)
    if image_path:
        content.append(_scaled_image(image_path))
        if caption:
            content.append(Spacer(1, 5))
            content.append(Paragraph(f"<i>{escape(caption)}</i>", styles["Caption"]))
    else:
        content.append(Paragraph("<b>Figure missing:</b> expected image file was not found.", styles["Normal"]))

    content.append(Spacer(1, 10))
    content.append(_interpretation_box(interpretation, styles))
    story.append(KeepTogether(content))
    story.append(PageBreak())


def _report_progress(callback, value: int, message: str):
    if callback is not None:
        try:
            callback(value, message)
        except Exception:
            pass


# =====================================================
# MAIN REPORT
# =====================================================

def generate_pdf_report(
    adata,
    project_name,
    author,
    institution,
    year,
    dataset_name,
    output_path,
    opcode_version,
    validation_metrics=None,
    vis_score=None,
    progress_callback=None,
):
    styles = _make_styles()
    story = []

    _report_progress(progress_callback, 5, "Preparing PDF layout...")

    page1_text = """The V5 score measures lineage dominance of oligodendrocyte transcriptional programs relative to competing neuronal and immune gene expression.

• Higher positive V5 scores indicate strong oligodendrocyte lineage dominance and a higher likelihood that a cell belongs to the OPC population. Accordingly, the cell with the maximum V5 score and the cluster with the highest mean V5 score represent the most confidently identifiable OPC candidates.

• Scores near zero indicate transcriptionally mixed lineage states, where oligodendrocyte-associated genes and competing lineage genes are expressed at similar levels within the same cell or cluster.

• Negative scores indicate non-OPC identity, as genes associated with other cell types dominate the transcriptional profile.

The following statistics summarize the distribution of V5 scores across all cells in the dataset."""

    page2_text = """This histogram shows the distribution of V5 scores across all cells in the dataset.

Cells on the right side of the distribution exhibit strong oligodendrocyte lineage dominance and therefore represent candidate OPC cells.

Cells with scores near zero represent transcriptionally mixed identities, where OPC and non-OPC gene expression programs are present simultaneously.

Cells with negative scores show dominant expression of non-OPC genes and are therefore unlikely to represent oligodendrocyte progenitor cells.

The rightmost tail of the distribution corresponds to the highest-confidence OPC cells identified by the algorithm. These cells are typically concentrated within the cluster that has the highest average V5 score."""

    page3_text = """Clusters with higher mean V5 scores contain cells with stronger oligodendrocyte lineage dominance and therefore indicate a greater presence of OPC cells.

The cluster with the highest mean V5 score represents the predicted OPC population within the dataset.

Clusters with lower scores represent neuronal, immune, or transcriptionally mixed cell populations."""

    page4_text = """Each point represents a single cell projected into transcriptional similarity space using UMAP dimensionality reduction.

Cells are grouped into clusters based on gene expression similarity.

When a dataset contains many clusters, the top 10 clusters by mean V5 score are highlighted in color, while all remaining clusters are shown in light gray to preserve global structure without overloading the figure.

The OPCODE scoring algorithm evaluates these clusters to determine which population exhibits the strongest oligodendrocyte lineage dominance using the OPC transcriptional signature previously developed from labeled datasets containing the top 60 genes most consistently expressed in OPC cells."""

    page5_text = """Cells are colored according to their V5 score.

Higher values — shown as teal-green colors rather than blue-purple — indicate stronger oligodendrocyte lineage identity.

Clusters enriched in high scores represent the predicted OPC population, while clusters with lower scores correspond to competing neuronal or immune transcriptional programs."""

    page6_text = """This figure highlights the cluster identified by OPCODE as the predicted OPC population because it contains the highest average V5 score among all clusters.

Cells within this cluster are transcriptionally similar and show strong enrichment of genes biologically associated with OPC identity while simultaneously exhibiting minimal expression of competing neuronal and immune gene programs.

This cluster therefore represents the primary pool of OPC candidate cells, which are subsequently analyzed at the single-cell level to identify the highest-confidence OPC cells and derive marker signatures used for experimental purification strategies such as FACS and MACS."""

    page7_text = page7_text = """Red points represent the highest-confidence OPC candidate cells identified by the OPCODE scoring algorithm.

Each point corresponds to a single cell in the UMAP embedding.

These cells exhibit the strongest oligodendrocyte lineage dominance while minimizing competing neuronal and immune gene expression programs.

A sparse red pattern indicates that only a small number of individual cells met the stringent high-confidence cutoff used for downstream marker discovery.

When these red cells localize to a coherent transcriptional region, this supports the interpretation that OPCODE has identified a biologically consistent high-confidence OPC subset rather than diffuse false-positive signal.

These high-confidence OPC cells form the basis for marker discovery and downstream experimental purification design."""

    page8_text = """This table summarizes genes that are significantly enriched in the high-confidence OPC cells highlighted in the previous UMAP.

The metrics used to evaluate marker suitability include:

• Mean expression within high-confidence OPC cells
• Detection rate (fraction of cells expressing the gene)
• Log fold-change measuring enrichment relative to the rest of the dataset

Genes that display high expression, strong detection frequency, and positive fold-change represent reliable transcriptional signals associated with OPC identity.

These genes serve as candidates for identifying cell-surface proteins that may be experimentally targeted using antibody-based purification approaches such as fluorescence-activated cell sorting (FACS) or magnetic-activated cell sorting (MACS)."""

    page9_text = """This ranking prioritizes genes based on their ability to distinguish high-confidence OPC cells from all other cells in the dataset.

Genes ranked at the top of the list typically exhibit:

• Strong enrichment in OPC cells
• High detection frequency among confident OPC candidates
• Minimal expression in non-OPC populations

These characteristics make the highest-ranked genes strong candidates for biological validation and experimental targeting during OPC purification."""

    page10_text = """This chart presents the recommended marker configuration for experimental OPC purification using FACS or MACS.

Markers are categorized into three functional groups:

• Positive surface markers: Proteins predicted to be expressed on the surface of OPC cells and suitable for direct antibody-based gating during sorting.

• RNA validation markers: Genes used to confirm OPC identity after sorting through downstream validation assays such as qPCR, immunocytochemistry, or single-cell sequencing.

• Negative lineage markers: Genes associated with competing cell types such as neurons or immune cells, which can be used to exclude contaminating populations during sorting.

The combination of these markers allows researchers to define a sorting phenotype that enriches OPC cells while minimizing contamination from competing cell types.

This marker panel therefore translates the computational predictions of OPCODE into a practical experimental framework for isolating high-purity OPC populations from heterogeneous mouse brain datasets.

In datasets where no surface marker passes evidence thresholds, experimentally established OPC anchors such as PDGFRα or NG2 remain reasonable candidates to test in follow-up purification experiments, even if RNA support is limited in that specific run due to transcript dropout. Markers such as Sox10, Olig1, and Olig2 should be interpreted as transcriptional lineage markers and RNA validation markers used to confirm OPC identity after purification, rather than as direct surface sorting targets."""

    story.append(Paragraph(f"<b>{escape(project_name)} {escape(opcode_version)} OPC Detection Report</b>", styles["CenterTitle"]))
    story.append(Paragraph(f"<b>Dataset:</b> {escape(dataset_name)}", styles["Small"]))
    story.append(Paragraph(f"<b>Author:</b> {escape(author)}", styles["Small"]))
    story.append(Paragraph(f"<b>Institution:</b> {escape(institution)}", styles["Small"]))
    story.append(Paragraph(f"<b>Year:</b> {escape(str(year))}", styles["Small"]))
    story.append(Spacer(1, 10))
    story.append(Paragraph("<b>OPCODE V5 Score Statistics (Single-Cell Level)</b>", styles["SectionTitle"]))
    story.append(_summary_stats_table(adata, validation_metrics, vis_score))
    story.append(Spacer(1, 12))
    story.append(_interpretation_box(page1_text, styles))
    story.append(PageBreak())
    _report_progress(progress_callback, 12, "Added summary statistics page...")

    _add_page_with_image(
        story, styles,
        "Distribution of V5 Scores",
        page2_text,
        [
            os.path.join("outputs", f"{dataset_name}_score_distribution.png"),
            os.path.join("outputs", "score_distribution.png"),
        ],
        caption="Score distribution across all cells in the dataset.",
    )
    _report_progress(progress_callback, 20, "Added V5 score distribution page...")

    _add_page_with_image(
        story, styles,
        "Mean V5 Score for Each Transcriptional Cluster",
        page3_text,
        [
            os.path.join("outputs", f"{dataset_name}_cluster_means.png"),
            os.path.join("outputs", "cluster_means.png"),
        ],
        caption="Cluster-level average V5 scores used to identify the predicted OPC population.",
    )
    _report_progress(progress_callback, 28, "Added cluster ranking page...")

    _add_page_with_image(
        story, styles,
        "UMAP Visualization Showing Transcriptional Clustering of Cells",
        page4_text,
        [
            os.path.join("outputs", f"{dataset_name}_umap_clusters.png"),
            os.path.join("outputs", "umap_clusters.png"),
        ],
        caption="UMAP projection of all cells, with the top V5-ranked clusters highlighted for readability when cluster counts are large.",
    )
    _report_progress(progress_callback, 36, "Added cluster UMAP page...")

    _add_page_with_image(
        story, styles,
        "UMAP Colored by V5 Score",
        page5_text,
        [
            os.path.join("outputs", f"{dataset_name}_umap_v5.png"),
            os.path.join("outputs", "umap_v5.png"),
        ],
        caption="Cells colored by single-cell V5 score.",
    )
    _report_progress(progress_callback, 44, "Added V5 UMAP page...")

    _add_page_with_image(
        story, styles,
        "UMAP Highlighting the Predicted OPC Cluster",
        page6_text,
        [
            os.path.join("outputs", f"{dataset_name}_umap_opc_cluster.png"),
        ],
        caption="Cluster-level localization of the predicted OPC-enriched population.",
    )
    _report_progress(progress_callback, 52, "Added predicted OPC cluster page...")

    _add_page_with_image(
        story, styles,
        "UMAP Highlighting High-Confidence OPC Candidate Cells",
        page7_text,
        [
            os.path.join("outputs", f"{dataset_name}_umap_highconf.png"),
            os.path.join("outputs", f"{dataset_name}_umap_highconf_cells.png"),
            os.path.join("outputs", "UMAP Highlighting High-Confidence OPC Candidate Cells.png"),
        ],
        caption="Single-cell refinement of the predicted OPC population.",
    )
    _report_progress(progress_callback, 60, "Added high-confidence OPC page...")

    story.append(Paragraph("<b>Marker Discovery</b>", styles["CenterTitle"]))
    marker_csv = os.path.join("outputs", "surface_marker_rankings.csv")
    marker_table = _marker_discovery_table(marker_csv)
    if marker_table is not None:
        story.append(marker_table)
    else:
        story.append(Paragraph("Marker discovery table could not be generated because the surface marker ranking file was not found.", styles["Normal"]))
    story.append(Spacer(1, 10))
    story.append(_interpretation_box(page8_text, styles))
    story.append(PageBreak())
    _report_progress(progress_callback, 70, "Added marker discovery page...")

    story.append(Paragraph("<b>Top 10 Ranked Marker Candidates</b>", styles["CenterTitle"]))
    rank_chart_path = _make_top_marker_chart(marker_csv)
    if rank_chart_path:
        story.append(_scaled_image(rank_chart_path))
        story.append(Spacer(1, 5))
        story.append(Paragraph("<i>Markers are labeled from 1–10 in descending order of ranking strength.</i>", styles["Caption"]))
    else:
        story.append(Paragraph("Top-marker ranking chart could not be generated because the marker ranking file was not found.", styles["Normal"]))
    story.append(Spacer(1, 10))
    story.append(_interpretation_box(page9_text, styles))
    story.append(PageBreak())
    _report_progress(progress_callback, 80, "Added marker ranking page...")

    story.append(Paragraph("<b>Experimental Marker Panel</b>", styles["CenterTitle"]))
    panel = _parse_sorting_panel(os.path.join("outputs", "sorting_panel.yaml"))

    if panel:
        species = panel.get("species", "unknown")
        dataset = panel.get("dataset", dataset_name)
        story.append(Paragraph(f"<b>Dataset:</b> {escape(str(dataset))} &nbsp;&nbsp;&nbsp; <b>Species:</b> {escape(str(species))}", styles["Small"]))
        story.append(Spacer(1, 6))

        positive_markers = []
        for item in _safe_list(panel.get("positive_surface_markers")):
            item = _safe_dict(item)
            gene = item.get("gene", "")
            protein = item.get("protein", "")
            notes = item.get("notes", "")
            txt = gene
            if protein:
                txt += f" — {protein}"
            if notes:
                txt += f" ({notes})"
            if txt.strip():
                positive_markers.append(txt)

        anchors = []
        for item in _safe_list(panel.get("anchors_to_test")):
            item = _safe_dict(item)
            gene = item.get("gene", "")
            protein = item.get("protein", "")
            supported = item.get("supported_by_rna_in_this_run", "")
            txt = gene
            if protein:
                txt += f" — {protein}"
            if supported != "":
                txt += f" | RNA-supported: {supported}"
            if txt.strip():
                anchors.append(txt)

        rna_markers = []
        for item in _safe_list(panel.get("rna_validation_markers")):
            item = _safe_dict(item)
            gene = item.get("gene", "")
            notes = item.get("notes", "")
            txt = gene
            if notes:
                txt += f" — {notes}"
            if txt.strip():
                rna_markers.append(txt)

        negative_markers = [str(x) for x in _safe_list(panel.get("negative_markers")) if str(x).strip()]
        gating_steps = [str(x) for x in _safe_list(panel.get("gating_steps")) if str(x).strip()]
        gating_steps = _clean_gating_steps(gating_steps)

        sections = [
            _simple_list_table("Positive Surface Markers", positive_markers, styles),
            _simple_list_table("Anchors to Test", anchors, styles),
            _simple_list_table("RNA Validation Markers", rna_markers, styles),
            _simple_list_table("Negative Lineage Markers", negative_markers, styles),
            _simple_list_table("Suggested Gating Steps", gating_steps, styles),
        ]

        for tbl in sections:
            if tbl is not None:
                story.append(tbl)
                story.append(Spacer(1, 3))

    else:
        story.append(Paragraph("Experimental marker panel could not be parsed because sorting_panel.yaml was not found or could not be read.", styles["Normal"]))

    story.append(Spacer(1, 6))
    story.append(_interpretation_box(page10_text, styles))
    story.append(PageBreak())
    _report_progress(progress_callback, 90, "Added experimental marker panel page...")

    story.append(Paragraph("<b>References</b>", styles["CenterTitle"]))
    story.append(
        Paragraph(
            "Placeholder for manuscript citations and supporting references.",
            styles["InterpretationText"],
        )
    )

    doc = SimpleDocTemplate(
        output_path,
        pagesize=letter,
        rightMargin=0.55 * inch,
        leftMargin=0.55 * inch,
        topMargin=0.55 * inch,
        bottomMargin=0.55 * inch,
    )

    _report_progress(progress_callback, 96, "Finalizing PDF file...")
    doc.build(story)
    _report_progress(progress_callback, 100, "PDF generation complete.")
    return output_path