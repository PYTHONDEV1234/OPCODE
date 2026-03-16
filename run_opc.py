import argparse
import os
import shutil
import importlib.metadata as imd

import scanpy as sc
import numpy as np
import pandas as pd

from scoring_engine.v5_opc_scoring_engine import score_v5_opc
from reporting.pdf_report import generate_pdf_report

from analysis.visualization import (
    plot_umap_by_cluster,
    plot_umap_by_v5_score,
    plot_umap_highlight_top_clusters
)

from analysis.purification import (
    add_v5_normalized_scores,
    build_purification_outputs
)

from analysis.surface_panel import build_surface_panel_outputs

# Try preferred validation API (matches Streamlit)
try:
    from analysis.validation_metrics import compute_validation_report as _compute_validation
except Exception:
    # Fallback if older API exists
    try:
        from analysis.validation_metrics import compute_validation_metrics as _compute_validation
    except Exception:
        _compute_validation = None


# -------------------------------------------------------
# HELPERS
# -------------------------------------------------------

def should_use_backed(input_path, size_threshold_gb=2):
    file_size_gb = os.path.getsize(input_path) / (1024 ** 3)
    return file_size_gb >= size_threshold_gb


def detect_cluster_column(adata, user_cluster_col: str | None):
    if user_cluster_col:
        if user_cluster_col not in adata.obs.columns:
            raise ValueError(f"--cluster-col '{user_cluster_col}' not found in adata.obs")
        return user_cluster_col
    if "cluster_alias" in adata.obs.columns:
        return "cluster_alias"
    if "clusters" in adata.obs.columns:
        return "clusters"
    return None


def _version(pkg: str) -> str:
    try:
        return imd.version(pkg)
    except Exception:
        return "unknown"


def _read_cell_ids_from_csv(path: str) -> list[str]:
    if not os.path.exists(path):
        return []
    df = pd.read_csv(path)
    if "cell_id" not in df.columns:
        return []
    return df["cell_id"].astype(str).tolist()


def _safe_remove(path: str):
    """Try to remove a file if it exists (useful when regenerating PDFs on Windows)."""
    try:
        if os.path.exists(path):
            os.remove(path)
    except Exception:
        pass


def _tagged_filename(filename: str, tag: str | None) -> str:
    """
    Insert _<tag> before file extension.
    example: foo.csv + run1 -> foo_run1.csv
    """
    if not tag:
        return filename
    base, ext = os.path.splitext(filename)
    return f"{base}_{tag}{ext}"


def _tagged_path(output_dir: str, filename: str, tag: str | None) -> str:
    return os.path.join(output_dir, _tagged_filename(filename, tag))


def _safe_move(src: str, dst: str, overwrite: bool = True) -> str:
    """
    Move/rename a file robustly on Windows.
    Returns dst if moved, otherwise returns src (if src missing).
    """
    if not src or not os.path.exists(src):
        return src
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if overwrite and os.path.exists(dst):
        _safe_remove(dst)
    try:
        shutil.move(src, dst)
        return dst
    except Exception:
        # fallback: copy + remove
        try:
            shutil.copy2(src, dst)
            _safe_remove(src)
            return dst
        except Exception:
            return src


# -------------------------------------------------------
# MAIN
# -------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="OPCODE — V5 CLI")

    parser.add_argument("--input", required=True, help="Path to .h5ad file")
    parser.add_argument("--output", required=True, help="Output directory")

    parser.add_argument("--large", action="store_true", help="Force backed mode")
    parser.add_argument("--cluster-col", default=None, help="Cluster column in adata.obs (default: auto-detect)")

    parser.add_argument("--umap", action="store_true", help="Generate UMAP figures")
    parser.add_argument("--validate", action="store_true", help="Compute validation metrics (if possible)")
    parser.add_argument("--report", action="store_true", help="Generate scientific PDF report")

    # Run tagging / suffixing
    parser.add_argument(
        "--run-tag",
        type=str,
        default=None,
        help="Tag appended to ALL output files (recommended for reproducible runs). Example: run1"
    )
    # Back-compat with your current flag naming
    parser.add_argument(
        "--pdf-suffix",
        type=str,
        default=None,
        help="Alias of --run-tag (kept for backwards compatibility)."
    )
    parser.add_argument(
        "--overwrite-pdf",
        action="store_true",
        help="If set, will attempt to overwrite existing tagged PDF (otherwise it will auto-rename)."
    )

    # Purification Mode
    parser.add_argument(
        "--purification-mode",
        action="store_true",
        help="Export OPC candidate + high-confidence cell lists (for purification workflows)"
    )
    parser.add_argument(
        "--top-k-clusters",
        type=int,
        default=1,
        help="How many top-ranked clusters to include as OPC candidates (default: 1)"
    )

    # High-confidence selection controls
    parser.add_argument(
        "--top-pct",
        type=float,
        default=10.0,
        help="If --v5-threshold is not set: keep the top X%% of candidate cells by V5_0_1 (default: 10)"
    )
    parser.add_argument(
        "--v5-threshold",
        type=float,
        default=None,
        help="Absolute threshold on V5_0_1 for high-confidence gating (overrides --top-pct). Example: 0.85"
    )
    parser.add_argument(
        "--marker-gate-quantile",
        type=float,
        default=0.60,
        help="Marker-score gate quantile within V5-filtered candidates (default: 0.60 keeps top 40%% by marker score)"
    )
    parser.add_argument(
        "--marker-gate-off",
        action="store_true",
        help="Disable marker-score gating (NOT recommended unless markers are missing)"
    )

    # QC controls
    parser.add_argument("--qc-off", action="store_true", help="Disable QC filtering in purification mode")
    parser.add_argument("--min-genes", type=int, default=200, help="QC: minimum genes per cell (if available)")
    parser.add_argument("--min-counts", type=int, default=None, help="QC: minimum total counts per cell (if available)")
    parser.add_argument("--max-mito", type=float, default=20.0, help="QC: maximum mito%% per cell (if available)")
    parser.add_argument("--mito-col", type=str, default=None, help="QC: mito%% column name (default: auto-detect)")

    # Optional strict exclusion
    parser.add_argument(
        "--strict-mature-exclusion",
        action="store_true",
        help="Exclude cells expressing Mbp/Plp1 (>0). More strict, may reduce yield."
    )

    # Sorting-panel outputs
    parser.add_argument(
        "--species",
        type=str,
        default="mouse",
        choices=["mouse", "human"],
        help="Species for surface marker DB defaults (default: mouse)"
    )
    parser.add_argument(
        "--surface-db",
        type=str,
        default=None,
        help="Path to surface marker DB CSV (default: config/surface_markers_mouse.csv for mouse)"
    )
    parser.add_argument(
        "--n-positive-markers",
        type=int,
        default=4,
        help="How many surface-positive markers to propose in sorting_panel.yaml (default: 4)"
    )

    args = parser.parse_args()

    input_path = args.input
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)

    # unify run tagging: --run-tag wins, else --pdf-suffix
    run_tag = args.run_tag or args.pdf_suffix

    dataset_name = os.path.splitext(os.path.basename(input_path))[0]

    print("\n====================================")
    print("OPCODE — V5 CLI")
    print("====================================\n")

    print("Versions:")
    print(" - scanpy:", _version("scanpy"))
    print(" - numpy :", _version("numpy"))
    print(" - pandas:", _version("pandas"))
    if run_tag:
        print(" - run tag:", run_tag)
    print("")

    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")

    # ---------------------------------------------------
    # Load dataset
    # ---------------------------------------------------
    use_backed = args.large or should_use_backed(input_path)

    if use_backed:
        print("Large dataset detected → Using backed mode (read-only)")
        adata = sc.read_h5ad(input_path, backed="r")
    else:
        print("Using in-memory mode")
        adata = sc.read_h5ad(input_path)

    print("\nDataset summary:")
    print(adata)

    cluster_col = detect_cluster_column(adata, args.cluster_col)
    if cluster_col is None:
        print("\n[WARN] No cluster column detected (expected 'clusters' or 'cluster_alias').")
        print("       Cluster ranking + highlighting will be limited unless you supply --cluster-col.")

    # ---------------------------------------------------
    # Scoring
    # ---------------------------------------------------
    print("\nRunning V5 scoring...")
    adata = score_v5_opc(adata)

    # Always add normalized scores for cross-dataset comparisons
    add_v5_normalized_scores(adata)

    # ---------------------------------------------------
    # Save CSV (scored obs) — TAGGED
    # ---------------------------------------------------
    scored_csv_name = f"{dataset_name}_v5_scored_output.csv"
    csv_path = _tagged_path(output_dir, scored_csv_name, run_tag)

    obs_out = adata.obs.copy()
    obs_out.insert(0, "cell_id", obs_out.index.astype(str))
    obs_out.to_csv(csv_path, index=False)
    print("\nSaved scored dataset:", csv_path)

    # ---------------------------------------------------
    # Purification mode outputs + surface panel package
    # ---------------------------------------------------
    purification_summary = None
    surface_summary = None

    if args.purification_mode:
        print("\nPurification mode: generating OPC candidate + high-confidence exports...")

        purification_summary = build_purification_outputs(
            adata=adata,
            output_dir=output_dir,
            cluster_col=cluster_col,
            top_k_clusters=args.top_k_clusters,
            qc_on=(not args.qc_off),
            min_genes=args.min_genes,
            min_counts=args.min_counts,
            max_mito=args.max_mito,
            mito_col=args.mito_col,
            v5_threshold=args.v5_threshold,
            top_pct=args.top_pct,
            marker_gate_quantile=args.marker_gate_quantile,
            use_marker_gate=(not args.marker_gate_off),
            strict_mature_exclusion=args.strict_mature_exclusion,
        )

        # TAG/RENAME purification outputs to prevent overwrites
        if run_tag:
            cand_src = purification_summary.get("opc_candidates_path")
            high_src = purification_summary.get("opc_highconf_path")

            cand_dst = _tagged_path(output_dir, "opc_candidates.csv", run_tag)
            high_dst = _tagged_path(output_dir, "opc_highconf.csv", run_tag)

            purification_summary["opc_candidates_path"] = _safe_move(cand_src, cand_dst, overwrite=True)
            purification_summary["opc_highconf_path"] = _safe_move(high_src, high_dst, overwrite=True)

        print("Purification exports written:")
        print(" -", purification_summary["opc_candidates_path"])
        print(" -", purification_summary["opc_highconf_path"])
        print("High-confidence fraction (post-QC candidates):",
              f"{purification_summary['highconf_fraction_of_candidates']:.3f}")

        if purification_summary.get("used_v5_threshold") is not None:
            print("Highconf selection: V5_0_1 >= ", purification_summary["used_v5_threshold"])
        else:
            print("Highconf selection: top", purification_summary.get("used_top_pct", args.top_pct), "% by V5_0_1")

        print("Marker gate enabled:", purification_summary.get("used_marker_gate", True))
        if purification_summary.get("used_marker_gate", True):
            print("Marker gate quantile:", purification_summary.get("marker_gate_quantile", args.marker_gate_quantile))
        print("Strict mature exclusion:", purification_summary.get("strict_mature_exclusion", False))
        print("Mito column used:", purification_summary.get("mito_col_used", None))

        # --- Surface panel package ---
        highconf_ids = _read_cell_ids_from_csv(purification_summary["opc_highconf_path"])

        if len(highconf_ids) == 0:
            print("\n[WARN] opc_highconf.csv is empty; skipping surface marker ranking + sorting_panel.yaml.")
        else:
            surface_db = args.surface_db
            if surface_db is None:
                if args.species == "mouse":
                    surface_db = os.path.join("config", "surface_markers_mouse.csv")
                else:
                    surface_db = None

            print("\nBuilding surface marker rankings + sorting panel proposal...")
            surface_summary = build_surface_panel_outputs(
                adata=adata,
                output_dir=output_dir,
                dataset_name=dataset_name,
                highconf_cell_ids=highconf_ids,
                species=args.species,
                surface_db_path=surface_db,
                n_positive=args.n_positive_markers,
            )

            # TAG/RENAME surface outputs too
            if run_tag:
                sm_src = surface_summary.get("surface_marker_rankings_path")
                sp_src = surface_summary.get("sorting_panel_path")

                # preserve original filenames but tag them
                sm_dst = _tagged_path(output_dir, os.path.basename(sm_src), run_tag) if sm_src else None
                sp_dst = _tagged_path(output_dir, os.path.basename(sp_src), run_tag) if sp_src else None

                if sm_src and sm_dst:
                    surface_summary["surface_marker_rankings_path"] = _safe_move(sm_src, sm_dst, overwrite=True)
                if sp_src and sp_dst:
                    surface_summary["sorting_panel_path"] = _safe_move(sp_src, sp_dst, overwrite=True)

            print("Surface outputs written:")
            print(" -", surface_summary["surface_marker_rankings_path"])
            print(" -", surface_summary["sorting_panel_path"])

            print("Top proposed surface-positive genes:",
                  ", ".join(surface_summary.get("positive_surface_genes", [])))
            print("RNA validation markers:",
                  ", ".join(surface_summary.get("rna_validation_markers", [])))
            print("Anchors to test:",
                  ", ".join(surface_summary.get("anchors_to_test", [])))
            if surface_summary.get("negative_markers"):
                print("Negative depletion markers present:",
                      ", ".join(surface_summary.get("negative_markers", [])))

    # ---------------------------------------------------
    # UMAP (Optional) — TAGGED
    # ---------------------------------------------------
    if args.umap:
        print("\nGenerating UMAP visualizations...")

        if "X_umap" not in adata.obsm:
            if "X_pca" not in adata.obsm:
                sc.pp.pca(adata)
            if "neighbors" not in adata.uns:
                sc.pp.neighbors(adata)
            sc.tl.umap(adata, random_state=42)

        if cluster_col is not None:
            fig_cluster = plot_umap_by_cluster(adata, cluster_col)
            cluster_name = f"{dataset_name}_umap_clusters.png"
            cluster_path = _tagged_path(output_dir, cluster_name, run_tag)
            fig_cluster.savefig(cluster_path, dpi=300, bbox_inches="tight")
            print("Saved:", cluster_path)

            fig_highlight = plot_umap_highlight_top_clusters(adata, cluster_col)
            highlight_name = f"{dataset_name}_umap_highlight.png"
            highlight_path = _tagged_path(output_dir, highlight_name, run_tag)
            fig_highlight.savefig(highlight_path, dpi=300, bbox_inches="tight")
            print("Saved:", highlight_path)

        fig_score = plot_umap_by_v5_score(adata, "V5_OPC_score")
        score_name = f"{dataset_name}_umap_v5.png"
        score_path = _tagged_path(output_dir, score_name, run_tag)
        fig_score.savefig(score_path, dpi=300, bbox_inches="tight")
        print("Saved:", score_path)

        print("UMAP figures saved.")

    # ---------------------------------------------------
    # Validation (Optional) — TAGGED
    # ---------------------------------------------------
    validation_metrics = None
    vis_score = None

    if args.validate:
        if _compute_validation is None:
            print("\n[WARN] Validation module not available (no compute_validation_report/metrics found). Skipping.")
        else:
            print("\nComputing validation metrics...")
            try:
                if _compute_validation.__name__ == "compute_validation_report":
                    validation_metrics = _compute_validation(
                        adata,
                        cluster_column=cluster_col or "clusters",
                        class_column="class"
                    )
                    vis_score = validation_metrics.get("VIS", None) if isinstance(validation_metrics, dict) else None
                else:
                    validation_metrics = _compute_validation(adata)

                metrics_name = f"{dataset_name}_validation_metrics.csv"
                metrics_path = _tagged_path(output_dir, metrics_name, run_tag)
                pd.DataFrame([validation_metrics]).to_csv(metrics_path, index=False)
                print("Validation metrics saved:", metrics_path)
            except Exception as e:
                print(f"[WARN] Validation failed: {e}")

    # ---------------------------------------------------
    # PDF Report (Optional) — TAGGED
    # ---------------------------------------------------
    if args.report:
        print("\nGenerating scientific report...")
        pdf_name = f"{dataset_name}_V5_OPC_Report.pdf"
        pdf_path = _tagged_path(output_dir, pdf_name, run_tag)

        if args.overwrite_pdf:
            _safe_remove(pdf_path)
        else:
            # if exists, auto-bump name
            if os.path.exists(pdf_path):
                base, ext = os.path.splitext(pdf_path)
                i = 2
                while os.path.exists(f"{base}_v{i}{ext}"):
                    i += 1
                pdf_path = f"{base}_v{i}{ext}"

        try:
            generate_pdf_report(
                adata=adata,
                project_name="OPCODE",
                author="Ansel Belani",
                institution="Delaware County Community College",
                year="2026",
                dataset_name=dataset_name,
                output_path=pdf_path,
                validation_metrics=validation_metrics,
                vis_score=vis_score,
                opcode_version="1.0.1",
            )
            print("Report saved:", pdf_path)
        except Exception as e:
            print(f"[WARN] PDF report generation failed: {e}")
            print("       If you're on Windows, close the PDF in any viewer and rerun.")
            print("       If upload is failing, try copying the PDF to a new name before uploading.")

    # ---------------------------------------------------
    print("\n====================================")
    print("Process Complete")
    print("====================================")
    print("Total cells:", adata.n_obs)
    print("Mean V5 score:", float(np.mean(adata.obs["V5_OPC_score"])))
    print("Max V5 score:", float(np.max(adata.obs["V5_OPC_score"])))
    print("Output directory:", output_dir)
    if run_tag:
        print("Run tag:", run_tag)

    if purification_summary is not None:
        print("\nPurification summary:")
        print("Candidates:", purification_summary["n_candidates"])
        print("Candidates after QC:", purification_summary["n_candidates_after_qc"])
        print("High-confidence OPCs:", purification_summary["n_highconf"])

    if surface_summary is not None:
        print("\nSurface panel summary:")
        print("Rankings:", surface_summary["surface_marker_rankings_path"])
        print("Panel:", surface_summary["sorting_panel_path"])
        print("Top surface-positive genes:", ", ".join(surface_summary.get("positive_surface_genes", [])))
        print("RNA validation markers:", ", ".join(surface_summary.get("rna_validation_markers", [])))
        print("Anchors to test:", ", ".join(surface_summary.get("anchors_to_test", [])))
        if surface_summary.get("negative_markers"):
            print("Negative depletion markers present:", ", ".join(surface_summary.get("negative_markers", [])))


if __name__ == "__main__":
    main()