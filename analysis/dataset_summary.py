import pandas as pd
import numpy as np
from analysis.effect_size import cohens_d

def summarize_dataset(adata):
    """
    Compute high-level summary statistics.
    """
    scores = adata.obs["V5_OPC_score"].values

    summary = {
        "Total Cells": len(scores),
        "Mean Score": float(np.mean(scores)),
        "Max Score": float(np.max(scores)),
        "Min Score": float(np.min(scores)),
        "Dynamic Range": float(np.max(scores) - np.min(scores))
    }

    if "class" in adata.obs.columns:
        oligo_mask = adata.obs["class"].str.contains("Oligo", case=False, na=False)
        glut_mask = adata.obs["class"].str.contains("Glut", case=False, na=False)
        gaba_mask = adata.obs["class"].str.contains("GABA", case=False, na=False)

        if oligo_mask.sum() > 0 and glut_mask.sum() > 0:
            summary["Effect Oligo vs Glut"] = float(
                cohens_d(
                    adata.obs.loc[oligo_mask, "V5_OPC_score"],
                    adata.obs.loc[glut_mask, "V5_OPC_score"]
                )
            )

        if oligo_mask.sum() > 0 and gaba_mask.sum() > 0:
            summary["Effect Oligo vs GABA"] = float(
                cohens_d(
                    adata.obs.loc[oligo_mask, "V5_OPC_score"],
                    adata.obs.loc[gaba_mask, "V5_OPC_score"]
                )
            )

    return summary