import os
from typing import Optional, Dict, Any, List, Tuple

import numpy as np
import pandas as pd


# ============================================================
# Sparse-safe helpers
# ============================================================

def _is_sparse(x) -> bool:
    try:
        import scipy.sparse as sp
        return sp.issparse(x)
    except Exception:
        return False


def _to_1d(a) -> np.ndarray:
    try:
        return np.asarray(a).ravel()
    except Exception:
        return np.array(a).ravel()


def _mean_and_detect(x) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns (mean_expr, detect_rate) for each gene in matrix x (cells x genes).
    detect_rate: fraction of cells with expr > 0.
    """
    n = x.shape[0]
    if n == 0:
        return np.zeros(x.shape[1], dtype=float), np.zeros(x.shape[1], dtype=float)

    if _is_sparse(x):
        mean = _to_1d(x.mean(axis=0)).astype(float)
        nnz = _to_1d((x > 0).sum(axis=0)).astype(float)
        detect = nnz / float(n)
        return mean, detect

    x = np.asarray(x)
    mean = x.mean(axis=0).astype(float)
    detect = (x > 0).mean(axis=0).astype(float)
    return mean, detect


def _safe_log2fc(a_mean: np.ndarray, b_mean: np.ndarray, eps: float = 1e-3) -> np.ndarray:
    return np.log2((a_mean + eps) / (b_mean + eps))


def _read_surface_db(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    if "gene" not in cols:
        raise ValueError("surface marker DB must contain a 'gene' column.")
    gene_col = cols["gene"]
    protein_col = cols.get("protein", None)
    notes_col = cols.get("notes", None)

    out = pd.DataFrame()
    out["gene"] = df[gene_col].astype(str)
    out["protein"] = df[protein_col].astype(str) if protein_col else ""
    out["notes"] = df[notes_col].astype(str) if notes_col else ""
    out = out.dropna(subset=["gene"]).drop_duplicates(subset=["gene"])
    return out


def _default_negative_markers(species: str = "mouse") -> List[str]:
    if species.lower() == "human":
        return ["PTPRC", "PECAM1", "EPCAM"]
    return ["Ptprc", "Pecam1", "Epcam"]


def _write_yaml(path: str, payload: Dict[str, Any]):
    """
    Minimal YAML writer (no external deps). Handles nested dict/list/scalars.
    """
    def fmt_scalar(v):
        if v is None:
            return "null"
        if isinstance(v, bool):
            return "true" if v else "false"
        if isinstance(v, (int, float)):
            if isinstance(v, float):
                return f"{v:.6g}"
            return str(v)
        s = str(v)
        if any(ch in s for ch in [":", "{", "}", "[", "]", "#", "\n", "\r", "\t"]):
            return '"' + s.replace('"', '\\"') + '"'
        return s

    lines: List[str] = []

    def emit(obj, indent=0, key=None):
        pre = "  " * indent
        if isinstance(obj, dict):
            if key is not None:
                lines.append(f"{pre}{key}:")
                pre = "  " * (indent + 1)
            for k, v in obj.items():
                if isinstance(v, (dict, list)):
                    emit(v, indent + (1 if key is not None else 0), k)
                else:
                    lines.append(f"{pre}{k}: {fmt_scalar(v)}")
        elif isinstance(obj, list):
            if key is not None:
                lines.append(f"{pre}{key}:")
                pre = "  " * (indent + 1)
            for item in obj:
                if isinstance(item, dict):
                    lines.append(f"{pre}-")
                    for k, v in item.items():
                        if isinstance(v, (dict, list)):
                            emit(v, indent + 2, k)
                        else:
                            lines.append(f"{'  ' * (indent + 2)}{k}: {fmt_scalar(v)}")
                else:
                    lines.append(f"{pre}- {fmt_scalar(item)}")
        else:
            if key is None:
                lines.append(f"{pre}{fmt_scalar(obj)}")
            else:
                lines.append(f"{pre}{key}: {fmt_scalar(obj)}")

    emit(payload, indent=0)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


# ============================================================
# Core scoring / ranking
# ============================================================

def compute_surface_marker_rankings(
    adata,
    highconf_cell_ids: List[str],
    surface_db_path: Optional[str] = None,
    max_genes: int = 5000,
) -> pd.DataFrame:
    """
    Compare highconf vs rest and compute ranking metrics.

    Columns:
      mean_highconf, mean_rest
      detect_highconf, detect_rest
      log2fc_highconf_vs_rest
      specificity_detect_delta
      ranking_score (general)
    """
    if len(highconf_cell_ids) == 0:
        raise ValueError("highconf_cell_ids is empty; cannot build surface marker rankings.")

    highconf_cell_ids = [cid for cid in highconf_cell_ids if cid in adata.obs_names]
    if len(highconf_cell_ids) == 0:
        raise ValueError("None of the provided highconf_cell_ids exist in adata.obs_names.")

    is_high = adata.obs_names.isin(highconf_cell_ids)
    if is_high.sum() == 0:
        raise ValueError("highconf mask resulted in 0 cells.")
    if (~is_high).sum() == 0:
        raise ValueError("No 'rest' cells available (dataset only contains highconf cells).")

    X_high = adata[is_high].X
    X_rest = adata[~is_high].X

    mean_h, det_h = _mean_and_detect(X_high)
    mean_r, det_r = _mean_and_detect(X_rest)

    log2fc = _safe_log2fc(mean_h, mean_r, eps=1e-3)
    spec = det_h - det_r

    # general-purpose score (not the same as the "panel score")
    log2fc_clip = np.clip(log2fc, 0, 8)
    spec_clip = np.clip(spec, 0, 1)
    score = log2fc_clip * spec_clip * np.clip(det_h, 0, 1)

    df = pd.DataFrame({
        "gene": adata.var_names.astype(str),
        "mean_highconf": mean_h,
        "mean_rest": mean_r,
        "detect_highconf": det_h,
        "detect_rest": det_r,
        "log2fc_highconf_vs_rest": log2fc,
        "specificity_detect_delta": spec,
        "ranking_score": score,
    }).sort_values("ranking_score", ascending=False)

    if surface_db_path and os.path.exists(surface_db_path):
        surface = _read_surface_db(surface_db_path)
        df = df.merge(surface, on="gene", how="left")
        df["is_surface_db"] = df["protein"].notna() & (df["protein"].astype(str) != "")
    else:
        df["protein"] = ""
        df["notes"] = ""
        df["is_surface_db"] = False

    if max_genes and len(df) > int(max_genes):
        df = df.head(int(max_genes)).copy()

    return df


def _panel_score_row(r: pd.Series) -> float:
    """
    Sorting-oriented score:
      prioritize markers that are:
        - detectable in highconf
        - specific (detect_highconf - detect_rest)
        - enriched (log2fc > 0)
    """
    det_h = float(r["detect_highconf"])
    det_r = float(r["detect_rest"])
    spec = max(0.0, det_h - det_r)
    l2 = max(0.0, float(r["log2fc_highconf_vs_rest"]))
    return det_h * spec * l2


# ============================================================
# Panel building
# ============================================================

def build_sorting_panel(
    rankings: pd.DataFrame,
    adata,
    dataset_name: str,
    species: str = "mouse",
    n_positive: int = 4,
    # evidence thresholds (wet-lab realism)
    min_detect_highconf: float = 0.05,
    min_specificity: float = 0.05,
    require_positive_log2fc: bool = True,
    # anchors: always reported (even if fail evidence)
    anchor_markers: Optional[List[str]] = None,
    # RNA-only validation markers (NOT sortable)
    rna_validation_markers: Optional[List[str]] = None,
    # allow positive surface markers ONLY from DB
    require_surface_db_for_positive: bool = True,
) -> Dict[str, Any]:
    """
    Build a panel spec from rankings.

    Outputs:
      - positive_surface_markers: evidence-supported AND surface_db-listed (sortable candidates)
      - anchors_to_test: classic OPC anchors (PDGFRA/CSPG4) with support flags
      - rna_validation_markers: lineage validators (SOX10/OLIG1/OLIG2) for post-sort confirmation
      - negative_markers: depletion markers included only if detected in rest
      - gating_steps: realistic boolean expressions (AND/OR)
    """
    if anchor_markers is None:
        anchor_markers = ["Pdgfra", "Cspg4"] if species.lower() == "mouse" else ["PDGFRA", "CSPG4"]

    if rna_validation_markers is None:
        if species.lower() == "human":
            rna_validation_markers = ["SOX10", "OLIG1", "OLIG2", "NKX2-2", "ASCL1"]
        else:
            rna_validation_markers = ["Sox10", "Olig1", "Olig2", "Nkx2-2", "Ascl1"]

    df = rankings.copy()
    df["panel_score"] = df.apply(_panel_score_row, axis=1)

    # evidence filters
    keep = (df["detect_highconf"] >= float(min_detect_highconf)) & \
           ((df["detect_highconf"] - df["detect_rest"]) >= float(min_specificity))

    if require_positive_log2fc:
        keep &= (df["log2fc_highconf_vs_rest"] > 0)

    df_pass = df[keep].copy()

    # enforce "surface only" for positive gating
    if require_surface_db_for_positive:
        df_pass = df_pass[df_pass["is_surface_db"] == True].copy()

    # rank within pass set
    df_pass = df_pass.sort_values(
        by=["panel_score", "ranking_score"],
        ascending=[False, False]
    )

    top_pos = df_pass.head(int(n_positive)).copy()

    positive_surface = []
    for _, r in top_pos.iterrows():
        positive_surface.append({
            "gene": str(r["gene"]),
            "protein": str(r.get("protein", "")) if pd.notna(r.get("protein", "")) else "",
            "detect_highconf": float(r["detect_highconf"]),
            "detect_rest": float(r["detect_rest"]),
            "log2fc": float(r["log2fc_highconf_vs_rest"]),
            "specificity": float(r["detect_highconf"] - r["detect_rest"]),
            "panel_score": float(r["panel_score"]),
            "notes": str(r.get("notes", "")) if pd.notna(r.get("notes", "")) else "",
        })

    # anchors-to-test (reported even if they fail)
    anchors = []
    for g in anchor_markers:
        if g in df["gene"].values:
            rr = df.loc[df["gene"] == g].iloc[0]
            supported = (
                (float(rr["detect_highconf"]) >= float(min_detect_highconf)) and
                ((float(rr["detect_highconf"]) - float(rr["detect_rest"])) >= float(min_specificity)) and
                ((float(rr["log2fc_highconf_vs_rest"]) > 0) if require_positive_log2fc else True)
            )
            anchors.append({
                "gene": str(rr["gene"]),
                "protein": str(rr.get("protein", "")) if pd.notna(rr.get("protein", "")) else "",
                "detect_highconf": float(rr["detect_highconf"]),
                "detect_rest": float(rr["detect_rest"]),
                "log2fc": float(rr["log2fc_highconf_vs_rest"]),
                "specificity": float(rr["detect_highconf"] - rr["detect_rest"]),
                "supported_by_rna_in_this_run": bool(supported),
                "notes": str(rr.get("notes", "")) if pd.notna(rr.get("notes", "")) else "",
            })

    # RNA validation markers (lineage confirmation; not sorting)
    rna_valid = []
    for g in rna_validation_markers:
        if g in df["gene"].values:
            rr = df.loc[df["gene"] == g].iloc[0]
            rna_valid.append({
                "gene": str(rr["gene"]),
                "detect_highconf": float(rr["detect_highconf"]),
                "detect_rest": float(rr["detect_rest"]),
                "log2fc": float(rr["log2fc_highconf_vs_rest"]),
                "notes": "RNA validation marker (not a surface sort target)."
            })

    # Negative markers: include only if present AND meaningfully detected in "rest"
    neg_candidates = _default_negative_markers(species=species)
    negative = []
    for g in neg_candidates:
        if g in df["gene"].values:
            rr = df.loc[df["gene"] == g].iloc[0]
            if float(rr["detect_rest"]) >= 0.01:
                negative.append(g)

    # Build realistic gating strings
    gating_steps = [
        "1) FSC/SSC gate to remove debris.",
        "2) Singlet gate (FSC-A vs FSC-H) to remove doublets.",
        "3) Live/Dead stain gate to keep viable cells.",
    ]
    if negative:
        gating_steps.append(f"4) Depletion gate (exclude): {' OR '.join([m + '+' for m in negative])}")

    if positive_surface:
        if len(positive_surface) >= 2:
            base = f"{positive_surface[0]['gene']}+ AND {positive_surface[1]['gene']}+"
            extras = [p["gene"] + "+" for p in positive_surface[2:]]
            if extras:
                gating_steps.append(
                    "5) Positive gate (keep): "
                    + base
                    + "  (optional refine: AND "
                    + " AND ".join(extras)
                    + ")"
                )
            else:
                gating_steps.append(f"5) Positive gate (keep): {base}")
        else:
            gating_steps.append(f"5) Positive gate (keep): {positive_surface[0]['gene']}+")
    else:
        # fall back suggestion based on anchors_to_test
        anchors_supported = [a["gene"] for a in anchors if a.get("supported_by_rna_in_this_run", False)]
        if anchors_supported:
            gating_steps.append(
                "5) Positive gate (keep): "
                + anchors_supported[0]
                + "+  (anchor supported by RNA; consider adding another surface marker if available)"
            )
        else:
            gating_steps.append(
                "5) Positive gate: No surface markers passed evidence thresholds. "
                "Try anchors_to_test experimentally (PDGFRα/NG2) and/or expand surface DB."
            )

    panel = {
        "dataset": dataset_name,
        "species": species,
        "evidence_thresholds": {
            "min_detect_highconf": float(min_detect_highconf),
            "min_specificity": float(min_specificity),
            "require_positive_log2fc": bool(require_positive_log2fc),
            "require_surface_db_for_positive": bool(require_surface_db_for_positive),
        },
        "positive_surface_markers": positive_surface,
        "anchors_to_test": anchors,
        "rna_validation_markers": rna_valid,
        "negative_markers": negative,
        "gating_steps": gating_steps,
        "notes": [
            "positive_surface_markers are the ONLY markers proposed for direct FACS/MACS positive gating.",
            "anchors_to_test are classic OPC anchors; RNA dropout is common, so lack of RNA support does not disqualify the protein experimentally.",
            "rna_validation_markers are for confirming OPC identity after sorting (qPCR/ICC/scRNA), not for surface gating.",
            "If no positive_surface_markers appear, expand config/surface_markers_mouse.csv with additional known surface candidates.",
        ],
    }
    return panel


def build_surface_panel_outputs(
    adata,
    output_dir: str,
    dataset_name: str,
    highconf_cell_ids: List[str],
    species: str = "mouse",
    surface_db_path: Optional[str] = None,
    n_positive: int = 4,
) -> Dict[str, Any]:
    os.makedirs(output_dir, exist_ok=True)

    rankings = compute_surface_marker_rankings(
        adata=adata,
        highconf_cell_ids=highconf_cell_ids,
        surface_db_path=surface_db_path,
        max_genes=5000,
    )
    rankings_path = os.path.join(output_dir, "surface_marker_rankings.csv")
    rankings.to_csv(rankings_path, index=False)

    panel = build_sorting_panel(
        rankings=rankings,
        adata=adata,
        dataset_name=dataset_name,
        species=species,
        n_positive=n_positive,
    )
    panel_path = os.path.join(output_dir, "sorting_panel.yaml")
    _write_yaml(panel_path, panel)

    return {
        "surface_marker_rankings_path": rankings_path,
        "sorting_panel_path": panel_path,
        "n_ranked_genes": int(len(rankings)),
        "positive_surface_genes": [p["gene"] for p in panel.get("positive_surface_markers", [])],
        "anchors_to_test": [a["gene"] for a in panel.get("anchors_to_test", [])],
        "negative_markers": panel.get("negative_markers", []),
        "rna_validation_markers": [v["gene"] for v in panel.get("rna_validation_markers", [])],
    }