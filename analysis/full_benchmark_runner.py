import pandas as pd
from analysis.validation_metrics import (
    compute_effect_sizes,
    compute_dynamic_range,
    compute_validation_integrity_score,
    compute_cluster_rank_integrity
)

def run_full_benchmark(
    adata,
    dataset_name,
    score_column="V5_OPC_score",
    class_column="class",
    cluster_column=None
):
    print(f"\nRunning FULL benchmark for {dataset_name}")

    results = {}

    # ----------------------------------------
    # EFFECT SIZES
    # ----------------------------------------

    effect_results = compute_effect_sizes(
        adata,
        score_column=score_column,
        class_column=class_column
    )

    results.update(effect_results)

    # ----------------------------------------
    # DYNAMIC RANGE
    # ----------------------------------------

    dynamic_range = compute_dynamic_range(
        adata,
        score_column=score_column
    )

    results["Dynamic_Range"] = dynamic_range

    # ----------------------------------------
    # VALIDATION INTEGRITY SCORE
    # ----------------------------------------

    vis = compute_validation_integrity_score(
        results
    )

    results["Validation_Integrity_Score"] = vis

    # ----------------------------------------
    # CLUSTER RANK INTEGRITY (if available)
    # ----------------------------------------

    if cluster_column and cluster_column in adata.obs.columns:
        cris = compute_cluster_rank_integrity(
            adata,
            cluster_column=cluster_column,
            score_column=score_column
        )
        results.update(cris)

    results["Dataset"] = dataset_name

    return pd.DataFrame([results])