import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import entropy

from scoring_engine.v5_opc_scoring_engine import score_v5_opc
from analysis.purification import add_v5_normalized_scores


# -----------------------------
# CONFIG
# -----------------------------

OPC_GENES = ["Pdgfra", "Cspg4", "Sox10", "Olig1", "Olig2"]
NEURON_GENES = ["Snap25", "Syt1", "Rbfox3", "Tubb3"]
TOP_PERCENTILE = 0.95


# -----------------------------
# UTILITIES
# -----------------------------

def zscore(x):
    x = np.array(x)
    if np.std(x) == 0:
        return np.zeros_like(x)
    return (x - np.mean(x)) / np.std(x)


def canonical_average_score(adata, genes):
    valid = [g for g in genes if g in adata.var_names]
    if not valid:
        return np.zeros(adata.n_obs)
    X = adata[:, valid].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    return X.mean(axis=1)


def neuronal_contamination(adata, mask):
    valid = [g for g in NEURON_GENES if g in adata.var_names]
    if not valid:
        return 0
    X = adata[:, valid].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    return X[mask].mean()


def ensure_clusters(adata):
    """
    Ensures adata.obs contains 'clusters'.
    If not present, compute Leiden clusters.
    """
    if "clusters" in adata.obs.columns:
        return

    print("No 'clusters' column found. Computing Leiden clusters...")

    if "X_pca" not in adata.obsm:
        sc.pp.pca(adata)

    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata)

    sc.tl.leiden(adata, key_added="clusters")


def cluster_separation(adata, scores, cluster_col="clusters"):
    df = pd.DataFrame({
        "cluster": adata.obs[cluster_col],
        "score": scores
    })
    means = df.groupby("cluster", observed=True)["score"].mean().sort_values(ascending=False)
    if len(means) < 2:
        return 0
    return means.iloc[0] - means.iloc[1]


def cluster_entropy(adata, mask, cluster_col="clusters"):
    clusters = adata.obs.loc[mask, cluster_col]
    probs = clusters.value_counts(normalize=True)
    return entropy(probs)


# -----------------------------
# MAIN BENCHMARK FUNCTION
# -----------------------------

def benchmark_dataset(path):
    print(f"\nProcessing: {path}")
    adata = sc.read_h5ad(path)

    # Ensure clustering exists
    ensure_clusters(adata)

    # OPCODE
    adata = score_v5_opc(adata)
    add_v5_normalized_scores(adata)
    opcode_score = adata.obs["V5_OPC_score"]

    # Scanpy additive scoring
    sc.tl.score_genes(adata, OPC_GENES, score_name="scanpy_score")
    scanpy_score = adata.obs["scanpy_score"]

    # Canonical average
    canonical_score = canonical_average_score(adata, OPC_GENES)

    methods = {
        "OPCODE": opcode_score,
        "Scanpy": scanpy_score,
        "Canonical": canonical_score
    }

    results = []

    for name, score in methods.items():
        score_z = zscore(score)
        threshold = np.quantile(score_z, TOP_PERCENTILE)
        mask = score_z >= threshold

        sep = cluster_separation(adata, score_z)
        contam = neuronal_contamination(adata, mask)
        ent = cluster_entropy(adata, mask)

        results.append({
            "dataset": os.path.basename(path),
            "method": name,
            "delta_sep": sep,
            "neuronal_contamination": contam,
            "cluster_entropy": ent
        })

    return results


# -----------------------------
# RUN
# -----------------------------

if __name__ == "__main__":

    DATASETS = [
        r"C:\Users\ansel\Desktop\OPC Project\raw_datasets\GSM2906405_Brain1_processed.h5ad",
        r"C:\Users\ansel\Desktop\OPC Project\raw_datasets\GSM2906406_Brain2_processed.h5ad",
        r"C:\Users\ansel\Desktop\OPC Project\raw_datasets\GSE60361_processed.h5ad",
        r"C:\Users\ansel\Desktop\OPC Project\raw_datasets\GSE115746_processed.h5ad",
    ]

    all_results = []

    for ds in DATASETS:
        all_results.extend(benchmark_dataset(ds))

    df = pd.DataFrame(all_results)
    df.to_csv("benchmark_results.csv", index=False)

    print("\nBenchmark complete.")
    print(df)