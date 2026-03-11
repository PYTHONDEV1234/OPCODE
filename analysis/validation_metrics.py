import numpy as np
import pandas as pd


# =========================================================
# EFFECT SIZE (Cohen's d)
# =========================================================

def cohens_d(a, b):
    a = np.asarray(a)
    b = np.asarray(b)

    if len(a) == 0 or len(b) == 0:
        return np.nan

    pooled_std = np.sqrt(((a.std() ** 2) + (b.std() ** 2)) / 2)

    if pooled_std == 0:
        return np.nan

    return (a.mean() - b.mean()) / pooled_std


def compute_effect_sizes(adata, score_column="V5_OPC_score", class_column="class"):

    if class_column not in adata.obs.columns:
        return {"Oligo_vs_Glut_d": np.nan,
                "Oligo_vs_GABA_d": np.nan}

    oligo_mask = adata.obs[class_column].str.contains("Oligo", case=False, na=False)
    glut_mask = adata.obs[class_column].str.contains("Glut", case=False, na=False)
    gaba_mask = adata.obs[class_column].str.contains("GABA", case=False, na=False)

    oligo_scores = adata.obs.loc[oligo_mask, score_column]
    glut_scores = adata.obs.loc[glut_mask, score_column]
    gaba_scores = adata.obs.loc[gaba_mask, score_column]

    return {
        "Oligo_vs_Glut_d": cohens_d(oligo_scores, glut_scores),
        "Oligo_vs_GABA_d": cohens_d(oligo_scores, gaba_scores),
    }


# =========================================================
# DYNAMIC RANGE
# =========================================================

def compute_dynamic_range(adata, score_column="V5_OPC_score"):
    scores = adata.obs[score_column]
    return float(scores.max() - scores.min())


# =========================================================
# CLUSTER RANK INTEGRITY SCORE (CRIS)
# =========================================================

def compute_cluster_rank_integrity(
    adata,
    cluster_column,
    score_column="V5_OPC_score"
):
    if cluster_column not in adata.obs.columns:
        return {
            "Average_Oligo_Rank": np.nan,
            "Total_Clusters": np.nan,
            "CRIS": np.nan
        }

    cluster_means = (
        adata.obs
        .groupby(cluster_column)[score_column]
        .mean()
        .sort_values(ascending=False)
    )

    total_clusters = len(cluster_means)

    oligo_clusters = cluster_means[
        cluster_means.index.astype(str).str.contains("Oligo|OPC", case=False)
    ]

    if len(oligo_clusters) == 0:
        return {
            "Average_Oligo_Rank": np.nan,
            "Total_Clusters": total_clusters,
            "CRIS": np.nan
        }

    ranks = [
        cluster_means.index.get_loc(idx) + 1
        for idx in oligo_clusters.index
    ]

    avg_rank = np.mean(ranks)

    cris = 1 - (avg_rank / total_clusters)

    return {
        "Average_Oligo_Rank": float(avg_rank),
        "Total_Clusters": total_clusters,
        "CRIS": float(cris)
    }


# =========================================================
# VALIDATION INTEGRITY SCORE
# =========================================================

def compute_validation_integrity_score(metrics_dict):
    """
    Combines dynamic range + effect sizes into one dataset-level score.
    """

    dynamic_range = metrics_dict.get("Dynamic_Range", 0)

    d1 = metrics_dict.get("Oligo_vs_Glut_d", 0)
    d2 = metrics_dict.get("Oligo_vs_GABA_d", 0)

    effect_component = 0

    if not np.isnan(d1):
        effect_component += abs(d1)

    if not np.isnan(d2):
        effect_component += abs(d2)

    return float(dynamic_range + effect_component)