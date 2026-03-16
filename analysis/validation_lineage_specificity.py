import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scoring_engine.v5_opc_scoring_engine import score_v5_opc

DATASETS = [
r"C:\Users\ansel\Desktop\OPC Project\raw_datasets\GSE115746_processed.h5ad"
]

def map_lineage(label):

    label = str(label).lower()

    if "opc" in label:
        return "OPC"

    if "oligodendro" in label:
        return "Oligodendrocyte"

    if "astro" in label:
        return "Astrocyte"

    if "microglia" in label or "immune" in label:
        return "Immune"

    if "gaba" in label or "glutamatergic" in label or "neuron" in label:
        return "Neuron"

    return "Other"


all_scores = []

for path in DATASETS:

    print("\nProcessing:", path)

    adata = sc.read_h5ad(path)

    # normalize counts
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # run OPCODE scoring
    adata = score_v5_opc(adata)

    # determine annotation column
    if "cell_class" in adata.obs.columns:
        celltype_col = "cell_class"

    elif "cell_subclass" in adata.obs.columns:
        celltype_col = "cell_subclass"

    else:
        raise RuntimeError("No recognizable annotation column.")

    adata.obs["lineage"] = adata.obs[celltype_col].apply(map_lineage)

    df = pd.DataFrame({
        "lineage": adata.obs["lineage"],
        "V5_score": adata.obs["V5_OPC_score"]
    })

    all_scores.append(df)

df = pd.concat(all_scores)

df = df[df.lineage != "Other"]

sns.set(style="whitegrid")

plt.figure(figsize=(8,6))

sns.violinplot(
    data=df,
    x="lineage",
    y="V5_score",
    order=["OPC","Oligodendrocyte","Neuron","Astrocyte","Immune"]
)

plt.axhline(0, linestyle="--", color="black")

plt.title("Lineage Specificity of the OPCODE V5 Dominance Score")

plt.xlabel("Cell lineage")

plt.ylabel("OPCODE V5 score")

plt.tight_layout()

plt.savefig(
r"C:\Users\ansel\Desktop\OPCODE\figure_lineage_specificity.png",
dpi=300
)

plt.show()