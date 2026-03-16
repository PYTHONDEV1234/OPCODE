import numpy as np
import pandas as pd


# =====================================================
# Small helpers
# =====================================================

def _to_float_or_nan(x):
    try:
        if x is None:
            return np.nan
        fx = float(x)
        if np.isnan(fx):
            return np.nan
        return fx
    except Exception:
        return np.nan


def _has_clusters(adata, cluster_column: str) -> bool:
    return hasattr(adata, "obs") and cluster_column in adata.obs.columns


def _has_usable_class_labels(adata, class_column: str) -> bool:
    """
    "Usable" means: column exists AND has at least 2 non-empty classes.
    We do not assume any particular taxonomy beyond that.
    """
    if not hasattr(adata, "obs") or class_column not in adata.obs.columns:
        return False
    s = adata.obs[class_column].astype(str)
    s = s.replace(["nan", "None", ""], np.nan).dropna()
    return s.nunique() >= 2


# =====================================================
# Core metrics (unlabeled-friendly)
# =====================================================

def compute_dynamic_range(adata, cluster_column="clusters", score_column="V5_OPC_score"):
    """
    Dynamic range of cluster mean scores: max(mean_cluster) - min(mean_cluster)
    """
    if not _has_clusters(adata, cluster_column) or score_column not in adata.obs.columns:
        return np.nan

    cluster_means = (
        adata.obs.groupby(cluster_column, observed=True)[score_column]
        .mean()
        .sort_values(ascending=False)
    )
    if len(cluster_means) < 2:
        return np.nan
    return float(cluster_means.max() - cluster_means.min())


def compute_top_cluster_gap(adata, cluster_column="clusters", score_column="V5_OPC_score"):
    """
    Gap between top and second cluster mean.
    """
    if not _has_clusters(adata, cluster_column) or score_column not in adata.obs.columns:
        return np.nan

    cluster_means = (
        adata.obs.groupby(cluster_column, observed=True)[score_column]
        .mean()
        .sort_values(ascending=False)
    )
    if len(cluster_means) < 2:
        return np.nan
    return float(cluster_means.iloc[0] - cluster_means.iloc[1])


def compute_positive_fraction(adata, score_column="V5_OPC_score"):
    """
    Fraction of cells with score > 0.
    Reported as evidence; not treated as inherently "good."
    """
    if score_column not in adata.obs.columns:
        return np.nan
    v = adata.obs[score_column].astype(float).values
    if len(v) == 0:
        return np.nan
    return float(np.mean(v > 0))


# =====================================================
# Label-based metrics (only meaningful if labels exist)
# =====================================================

def cohens_d(a: np.ndarray, b: np.ndarray):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size < 2 or b.size < 2:
        return np.nan

    sa = np.nanstd(a, ddof=1)
    sb = np.nanstd(b, ddof=1)
    pooled = np.sqrt((sa ** 2 + sb ** 2) / 2.0)
    if pooled <= 1e-12 or np.isnan(pooled):
        return np.nan
    return float((np.nanmean(a) - np.nanmean(b)) / pooled)


def compute_effect_sizes(adata, class_column="class", score_column="V5_OPC_score"):
    """
    Conservative effect sizes focused on oligodendrocyte-lineage vs excitatory/inhibitory neurons
    IF the class labels contain those substrings. Otherwise returns NaNs.
    """
    out = {"Oligo_vs_Glut_d": np.nan, "Oligo_vs_GABA_d": np.nan}

    if not _has_usable_class_labels(adata, class_column) or score_column not in adata.obs.columns:
        return out

    cls = adata.obs[class_column].astype(str)
    scores = adata.obs[score_column].astype(float)

    oligo_mask = cls.str.contains("Oligo", case=False, na=False)
    glut_mask = cls.str.contains("Glut", case=False, na=False)
    gaba_mask = cls.str.contains("GABA", case=False, na=False)

    oligo = scores[oligo_mask].values
    glut = scores[glut_mask].values
    gaba = scores[gaba_mask].values

    out["Oligo_vs_Glut_d"] = cohens_d(oligo, glut)
    out["Oligo_vs_GABA_d"] = cohens_d(oligo, gaba)

    return out


def compute_cluster_rank_integrity(
    adata,
    cluster_column="clusters",
    class_column="class",
    score_column="V5_OPC_score",
):
    """
    CRIS: Cluster Rank Integrity Score.
    Interprets "Oligo*" as the target class (if present).
    If labels don't contain Oligo*, returns NaNs.

    Idea:
      - rank clusters by mean score (descending)
      - find clusters containing Oligo-labeled cells
      - compute average rank of those clusters
      - CRIS = 1 - (avg_rank - 1)/(total_clusters - 1)
    """
    if not _has_clusters(adata, cluster_column) or score_column not in adata.obs.columns:
        return {"Average_Oligo_Rank": np.nan, "Total_Clusters": np.nan, "CRIS": np.nan}

    cluster_means = (
        adata.obs.groupby(cluster_column, observed=True)[score_column]
        .mean()
        .sort_values(ascending=False)
    )
    total_clusters = int(len(cluster_means))
    if total_clusters < 2:
        return {"Average_Oligo_Rank": np.nan, "Total_Clusters": total_clusters, "CRIS": np.nan}

    if not _has_usable_class_labels(adata, class_column):
        return {"Average_Oligo_Rank": np.nan, "Total_Clusters": total_clusters, "CRIS": np.nan}

    cls = adata.obs[class_column].astype(str)
    oligo_clusters = (
        adata.obs.loc[cls.str.contains("Oligo", case=False, na=False), cluster_column]
        .astype(str)
        .unique()
    )
    if len(oligo_clusters) == 0:
        return {"Average_Oligo_Rank": np.nan, "Total_Clusters": total_clusters, "CRIS": np.nan}

    idx_as_str = cluster_means.index.astype(str)
    ranks = []
    for c in oligo_clusters:
        pos = np.where(idx_as_str == str(c))[0]
        if pos.size:
            ranks.append(int(pos[0]) + 1)

    if len(ranks) == 0:
        return {"Average_Oligo_Rank": np.nan, "Total_Clusters": total_clusters, "CRIS": np.nan}

    avg_rank = float(np.mean(ranks))
    cris = 1.0 - ((avg_rank - 1.0) / (total_clusters - 1.0))
    return {"Average_Oligo_Rank": avg_rank, "Total_Clusters": total_clusters, "CRIS": float(cris)}


# =====================================================
# Legacy: single dict (kept for compatibility)
# =====================================================

def compute_validation_metrics(
    adata,
    cluster_column="clusters",
    class_column="class",
    score_column="V5_OPC_score",
):
    """
    Legacy output: one flat dict.
    """
    dyn = compute_dynamic_range(adata, cluster_column=cluster_column, score_column=score_column)
    gap = compute_top_cluster_gap(adata, cluster_column=cluster_column, score_column=score_column)
    pos = compute_positive_fraction(adata, score_column=score_column)

    eff = compute_effect_sizes(adata, class_column=class_column, score_column=score_column)
    cris = compute_cluster_rank_integrity(
        adata,
        cluster_column=cluster_column,
        class_column=class_column,
        score_column=score_column,
    )

    return {
        "Dynamic_Range": _to_float_or_nan(dyn),
        "Top_Cluster_Gap": _to_float_or_nan(gap),
        "Positive_Fraction": _to_float_or_nan(pos),
        "Oligo_vs_Glut_d": _to_float_or_nan(eff.get("Oligo_vs_Glut_d")),
        "Oligo_vs_GABA_d": _to_float_or_nan(eff.get("Oligo_vs_GABA_d")),
        "CRIS": _to_float_or_nan(cris.get("CRIS")),
        "Average_Oligo_Rank": _to_float_or_nan(cris.get("Average_Oligo_Rank")),
        "Total_Clusters": _to_float_or_nan(cris.get("Total_Clusters")),
    }


def compute_validation_integrity_score(metrics: dict):
    """
    Legacy VIS: robust mean of available numeric metrics (ignores NaN/non-numeric).
    """
    if not isinstance(metrics, dict):
        return np.nan

    vals = []
    for v in metrics.values():
        fv = _to_float_or_nan(v)
        if np.isfinite(fv):
            vals.append(float(fv))
    return float(np.mean(vals)) if len(vals) else np.nan


# =====================================================
# Preferred: structured report (unlabeled vs labeled)
# =====================================================

def compute_validation_report(
    adata,
    cluster_column="clusters",
    class_column="class",
    score_column="V5_OPC_score",
):
    """
    Returns a structured report:
      - signal_metrics: unlabeled-friendly evidence
      - label_metrics: label-based evidence (NaN if not usable)
      - Signal_Score: simple transparent aggregate of signal_metrics
      - Label_Validation_Score: aggregate of label_metrics when usable
      - Combined_Score: average of the two when available

    Also includes flags about whether clusters/labels were usable.
    """

    # Handle cluster column fallback
    if cluster_column not in getattr(adata, "obs", {}):
        if "cluster_alias" in getattr(adata, "obs", {}):
            cluster_column = "cluster_alias"

    clusters_ok = _has_clusters(adata, cluster_column)
    labeled_ok = _has_usable_class_labels(adata, class_column)

    # Signal metrics (unlabeled OK)
    dyn = compute_dynamic_range(adata, cluster_column=cluster_column, score_column=score_column) if clusters_ok else np.nan
    gap = compute_top_cluster_gap(adata, cluster_column=cluster_column, score_column=score_column) if clusters_ok else np.nan
    pos = compute_positive_fraction(adata, score_column=score_column)

    signal_metrics = {
        "Dynamic_Range": _to_float_or_nan(dyn),
        "Top_Cluster_Gap": _to_float_or_nan(gap),
        "Positive_Fraction": _to_float_or_nan(pos),
    }

    # Transparent heuristic for a signal score
    sig_parts = []
    if np.isfinite(signal_metrics["Dynamic_Range"]):
        sig_parts.append(float(np.log1p(max(signal_metrics["Dynamic_Range"], 0.0))))
    if np.isfinite(signal_metrics["Top_Cluster_Gap"]):
        sig_parts.append(float(max(signal_metrics["Top_Cluster_Gap"], 0.0)))
    if np.isfinite(signal_metrics["Positive_Fraction"]):
        sig_parts.append(float(signal_metrics["Positive_Fraction"]))

    signal_score = float(np.mean(sig_parts)) if len(sig_parts) else np.nan

    # Label metrics (only if labeled)
    label_metrics = {"Oligo_vs_Glut_d": np.nan, "Oligo_vs_GABA_d": np.nan, "CRIS": np.nan}
    label_score = np.nan

    if labeled_ok:
        eff = compute_effect_sizes(adata, class_column=class_column, score_column=score_column)
        cris = compute_cluster_rank_integrity(
            adata,
            cluster_column=cluster_column,
            class_column=class_column,
            score_column=score_column,
        )
        label_metrics.update({
            "Oligo_vs_Glut_d": _to_float_or_nan(eff.get("Oligo_vs_Glut_d")),
            "Oligo_vs_GABA_d": _to_float_or_nan(eff.get("Oligo_vs_GABA_d")),
            "CRIS": _to_float_or_nan(cris.get("CRIS")),
        })

        lab_parts = [v for v in label_metrics.values() if np.isfinite(_to_float_or_nan(v))]
        label_score = float(np.mean(lab_parts)) if len(lab_parts) else np.nan

    # Combined score
    combined_parts = [v for v in [signal_score, label_score] if np.isfinite(_to_float_or_nan(v))]
    combined_score = float(np.mean(combined_parts)) if len(combined_parts) else np.nan

    return {
        "signal_metrics": signal_metrics,
        "label_metrics": label_metrics,
        "Signal_Score": _to_float_or_nan(signal_score),
        "Label_Validation_Score": _to_float_or_nan(label_score),
        "Combined_Score": _to_float_or_nan(combined_score),
        # flags / metadata
        "Has_Usable_Class_Labels": bool(labeled_ok),
        "Has_Clusters": bool(clusters_ok),
        "Class_Column": class_column if class_column in getattr(adata, "obs", {}) else None,
        "Cluster_Column": cluster_column if cluster_column in getattr(adata, "obs", {}) else None,
    }