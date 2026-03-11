import pandas as pd
import numpy as np

from analysis.validation_metrics import (
    compute_effect_sizes,
    compute_dynamic_range,
    compute_validation_integrity_score,
)

def run_benchmark(adata, dataset_name):

    print(f"\nRunning benchmark for {dataset_name}")

    results = []

    score_column = "V5_OPC_score"

    if score_column not in adata.obs.columns:
        raise ValueError("V5 score not found in AnnData object.")

    # -------------------------
    # Effect sizes
    # -------------------------

    effect_metrics = compute_effect_sizes(adata, score_column)

    if effect_metrics is None:
        effect_metrics = {
            "Oligo_vs_Glut_d": np.nan,
            "Oligo_vs_GABA_d": np.nan,
        }

    # -------------------------
    # Dynamic range
    # -------------------------

    dynamic_range = compute_dynamic_range(adata, score_column)

    if dynamic_range is None:
        dynamic_range = np.nan

    # -------------------------
    # Integrity score
    # -------------------------

    integrity_score = compute_validation_integrity_score(
        effect_metrics["Oligo_vs_Glut_d"],
        effect_metrics["Oligo_vs_GABA_d"],
        dynamic_range,
    )

    if integrity_score is None:
        integrity_score = np.nan

    # -------------------------
    # Assemble row
    # -------------------------

    results.append({
        "Dataset": dataset_name,
        "Oligo_vs_Glut_d": effect_metrics["Oligo_vs_Glut_d"],
        "Oligo_vs_GABA_d": effect_metrics["Oligo_vs_GABA_d"],
        "Dynamic_Range": dynamic_range,
        "Validation_Integrity_Score": integrity_score,
    })

    df = pd.DataFrame(results)

    return df