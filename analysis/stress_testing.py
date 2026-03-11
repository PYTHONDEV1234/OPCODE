import pandas as pd
import numpy as np

from analysis.validation_metrics import (
    compute_effect_sizes,
    compute_dynamic_range,
    compute_cluster_rank_integrity
)


def run_region_stress_test(
    adata,
    dataset_name,
    region_name,
    cluster_column="cluster_alias",
    class_column="class",
    score_column="V5_OPC_score"
):
    """
    Runs structured stress metrics on a single region dataset.
    """

    print(f"\nRunning stress test for {dataset_name} - {region_name}")

    # ---------------------------------
    # Effect sizes
    # ---------------------------------
    effect_metrics = compute_effect_sizes(
        adata,
        class_column=class_column,
        score_column=score_column
    )

    # ---------------------------------
    # Dynamic range
    # ---------------------------------
    dynamic_range = compute_dynamic_range(
        adata,
        score_column=score_column
    )

    # ---------------------------------
    # CRIS (if clusters available)
    # ---------------------------------
    if cluster_column in adata.obs.columns:
        cris_metrics = compute_cluster_rank_integrity(
            adata,
            cluster_column=cluster_column,
            score_column=score_column
        )
    else:
        cris_metrics = {
            "Average_Oligo_Rank": np.nan,
            "Total_Clusters": np.nan,
            "CRIS": np.nan
        }

    # ---------------------------------
    # Combine metrics
    # ---------------------------------
    result = {
        "Dataset": dataset_name,
        "Region": region_name,
        "Oligo_vs_Glut_d": effect_metrics.get("Oligo_vs_Glut_d"),
        "Oligo_vs_GABA_d": effect_metrics.get("Oligo_vs_GABA_d"),
        "Dynamic_Range": dynamic_range,
        "Average_Oligo_Rank": cris_metrics.get("Average_Oligo_Rank"),
        "Total_Clusters": cris_metrics.get("Total_Clusters"),
        "CRIS": cris_metrics.get("CRIS")
    }

    return result


def save_stress_results(results_list, output_path):
    df = pd.DataFrame(results_list)
    df.to_csv(output_path, index=False)
    print(f"\nSaved stress results to: {output_path}")