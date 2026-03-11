import numpy as np
import scanpy as sc
from analysis.validation_metrics import (
    compute_effect_sizes,
    compute_dynamic_range,
    compute_validation_integrity_score,
    compute_cluster_rank_integrity
)

# ============================================================
# RAW OPC AVERAGE (Naive Baseline)
# ============================================================

def score_raw_opc_average(adata):

    opc_genes = [
        "Pdgfra", "Cspg4", "Ptprz1",
        "Sox10", "Olig1", "Olig2"
    ]

    # Resolve gene naming automatically
    if "gene_symbol" in adata.var.columns:
        gene_names = adata.var["gene_symbol"].values
    else:
        gene_names = adata.var_names

    valid = [g for g in opc_genes if g in gene_names]

    if len(valid) == 0:
        adata.obs["RAW_OPC_score"] = 0
        return adata

    indices = [np.where(gene_names == g)[0][0] for g in valid]

    X = adata[:, indices].X

    if not isinstance(X, np.ndarray):
        X = X.toarray()

    adata.obs["RAW_OPC_score"] = X.mean(axis=1)

    return adata


# ============================================================
# SCANPY MODULE SCORE (Common Literature Approach)
# ============================================================

def score_scanpy_module(adata):

    opc_genes = [
        "Pdgfra", "Cspg4", "Ptprz1",
        "Sox10", "Olig1", "Olig2"
    ]

    if "gene_symbol" in adata.var.columns:
        gene_names = adata.var["gene_symbol"].values
    else:
        gene_names = adata.var_names

    valid = [g for g in opc_genes if g in gene_names]

    if len(valid) == 0:
        adata.obs["MODULE_OPC_score"] = 0
        return adata

    sc.tl.score_genes(
        adata,
        gene_list=valid,
        score_name="MODULE_OPC_score",
        use_raw=False
    )

    return adata


# ============================================================
# METHOD COMPARISON RUNNER
# ============================================================

def run_method_comparison(
    adata,
    dataset_name,
    cluster_column
):

    results = []

    methods = {
        "RAW": "RAW_OPC_score",
        "MODULE": "MODULE_OPC_score",
        "V5": "V5_OPC_score"
    }

    for method_name, score_column in methods.items():

        if score_column not in adata.obs.columns:
            continue

        print(f"\nEvaluating {method_name} on {dataset_name}")

        effect_metrics = compute_effect_sizes(
            adata,
            score_column=score_column
        )

        dynamic_range = compute_dynamic_range(
            adata,
            score_column=score_column
        )

        integrity_score = compute_validation_integrity_score(
            effect_metrics,
            dynamic_range
        )

        cris_data = compute_cluster_rank_integrity(
            adata,
            cluster_column=cluster_column,
            score_column=score_column
        )

        results.append({
            "Dataset": dataset_name,
            "Method": method_name,
            "Oligo_vs_Glut_d": effect_metrics["Oligo_vs_Glut_d"],
            "Oligo_vs_GABA_d": effect_metrics["Oligo_vs_GABA_d"],
            "Dynamic_Range": dynamic_range,
            "Validation_Integrity_Score": integrity_score,
            "CRIS": cris_data["CRIS"]
        })

    return results