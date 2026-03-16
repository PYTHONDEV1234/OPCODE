import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# -----------------------------
# CONFIG
# -----------------------------

INPUT_FILE = "benchmark_results.csv"
OUTPUT_FILE = "benchmark_figure_interpretable.png"

DATASET_NAME_MAP = {
    "GSM2906405_Brain1_processed.h5ad": "Brain1",
    "GSM2906406_Brain2_processed.h5ad": "Brain2",
    "GSE60361_processed.h5ad": "GSE60361",
    "GSE115746_processed.h5ad": "GSE115746"
}

METHOD_ORDER = ["OPCODE", "Scanpy", "Canonical"]

sns.set_style("white")

# -----------------------------
# LOAD
# -----------------------------

df = pd.read_csv(INPUT_FILE)
df["dataset_short"] = df["dataset"].map(DATASET_NAME_MAP)
df["method"] = pd.Categorical(df["method"], categories=METHOD_ORDER, ordered=True)

# -----------------------------
# FIGURE
# -----------------------------

fig = plt.figure(figsize=(14, 12))
gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.4])

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])
ax_caption = fig.add_subplot(gs[2, :])
ax_caption.axis("off")

# -----------------------------
# PANEL A — Δ_sep
# -----------------------------

sns.lineplot(
    data=df,
    x="dataset_short",
    y="delta_sep",
    hue="method",
    marker="o",
    ax=ax1
)

ax1.set_title("(A) Cluster Separation (Δ_sep)\nHigher = Stronger Top Cluster", fontsize=12)
ax1.set_xlabel("")
ax1.set_ylabel("Δ_sep")
ax1.tick_params(axis='x', rotation=0)

# -----------------------------
# PANEL B — Contamination
# -----------------------------

sns.lineplot(
    data=df,
    x="dataset_short",
    y="neuronal_contamination",
    hue="method",
    marker="o",
    ax=ax2,
    legend=False
)

ax2.set_title("(B) Neuronal Contamination\nLower = Purer OPC Selection", fontsize=12)
ax2.set_xlabel("")
ax2.set_ylabel("Mean Neuronal Marker Expression")
ax2.tick_params(axis='x', rotation=0)

# -----------------------------
# PANEL C — Entropy
# -----------------------------

sns.lineplot(
    data=df,
    x="dataset_short",
    y="cluster_entropy",
    hue="method",
    marker="o",
    ax=ax3,
    legend=False
)

ax3.set_title("(C) Cluster Entropy\nLower = More Cluster-Specific Selection", fontsize=12)
ax3.set_xlabel("")
ax3.set_ylabel("Entropy")
ax3.tick_params(axis='x', rotation=0)

# -----------------------------
# PANEL D — Summary
# -----------------------------

mean_df = df.groupby("method", observed=True)[
    ["delta_sep", "neuronal_contamination", "cluster_entropy"]
].mean().reset_index()

mean_melt = mean_df.melt(id_vars="method", var_name="metric", value_name="value")

sns.barplot(
    data=mean_melt,
    x="method",
    y="value",
    hue="metric",
    ax=ax4
)

ax4.set_title("(D) Cross-Dataset Mean Metrics", fontsize=12)
ax4.set_xlabel("")
ax4.set_ylabel("Mean Value")

# -----------------------------
# GLOBAL LEGEND
# -----------------------------

handles, labels = ax1.get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False)

ax1.legend_.remove()

# -----------------------------
# CAPTION
# -----------------------------

caption_text = (
    "OPCODE (dominance-based scoring) consistently reduces neuronal contamination and cluster entropy\n"
    "across four independent scRNA-seq datasets compared to additive scoring methods (Scanpy and Canonical).\n"
    "Although additive methods often show higher raw cluster separation (Δ_sep), this is accompanied by\n"
    "substantially increased lineage impurity and cluster dispersion."
)

ax_caption.text(
    0.5, 0.5,
    caption_text,
    ha="center",
    va="center",
    fontsize=11
)

plt.tight_layout(rect=[0, 0.05, 1, 0.95])
plt.savefig(OUTPUT_FILE, dpi=300)
plt.show()

print("Figure saved:", OUTPUT_FILE)