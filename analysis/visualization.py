import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

MAX_CLUSTER_LEGEND = 10
TOP_CLUSTER_LABELS = 10


# =====================================================
# HELPERS
# =====================================================

def _has_umap(adata):
    return "X_umap" in adata.obsm and adata.obsm["X_umap"] is not None


def _has_neighbors(adata):
    return "neighbors" in adata.uns


def _is_logged(adata):
    """
    Heuristic:
    if max expression is modest, data is probably already log-transformed.
    """
    try:
        x_max = adata.X.max()
        if hasattr(x_max, "item"):
            x_max = x_max.item()
        return x_max < 50
    except Exception:
        return True


def ensure_umap(adata):
    """
    Ensures UMAP coordinates exist.
    If missing, computes:
      - normalization (if needed)
      - log1p (if needed)
      - HVGs
      - PCA
      - neighbors
      - UMAP
    """
    if _has_umap(adata):
        return adata

    adata = adata.copy()

    try:
        total_counts = np.array(adata.X.sum(axis=1)).flatten()
        median_counts = float(np.median(total_counts))
    except Exception:
        median_counts = 1e4

    if median_counts > 100:
        sc.pp.normalize_total(adata, target_sum=1e4)

    if not _is_logged(adata):
        sc.pp.log1p(adata)

    if "highly_variable" not in adata.var.columns:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=min(2000, adata.shape[1]),
            flavor="seurat",
            subset=False,
        )

    hvg_mask = adata.var["highly_variable"].values
    if hvg_mask.sum() > 50:
        adata_for_pca = adata[:, hvg_mask].copy()
    else:
        adata_for_pca = adata

    if "X_pca" not in adata_for_pca.obsm:
        n_comps = min(50, max(2, adata_for_pca.shape[1] - 1))
        sc.pp.pca(adata_for_pca, n_comps=n_comps)

    if not _has_neighbors(adata_for_pca):
        sc.pp.neighbors(
            adata_for_pca,
            n_neighbors=15,
            n_pcs=min(30, adata_for_pca.obsm["X_pca"].shape[1]),
        )

    if not _has_umap(adata_for_pca):
        sc.tl.umap(adata_for_pca)

    if "X_pca" in adata_for_pca.obsm:
        adata.obsm["X_pca"] = adata_for_pca.obsm["X_pca"]
    if "neighbors" in adata_for_pca.uns:
        adata.uns["neighbors"] = adata_for_pca.uns["neighbors"]

    adata.obsm["X_umap"] = adata_for_pca.obsm["X_umap"]
    return adata


def _get_top_clusters_by_v5(adata, cluster_column, top_n=TOP_CLUSTER_LABELS):
    cluster_means = (
        adata.obs.groupby(cluster_column)["V5_OPC_score"]
        .mean()
        .sort_values(ascending=False)
    )
    return list(cluster_means.head(top_n).index.astype(str))


def _top_cluster_colors(n: int):
    cmap = plt.get_cmap("tab10")
    return [cmap(i) for i in range(min(n, 10))]


def _style_legend(legend):
    if legend is None:
        return
    handles = getattr(legend, "legend_handles", None)
    if handles is None:
        handles = getattr(legend, "legendHandles", [])
    for handle in handles:
        try:
            handle.set_sizes([60])
        except Exception:
            pass


# =====================================================
# PLOTTING
# =====================================================

def plot_umap_by_cluster(adata, cluster_column):
    adata = ensure_umap(adata)

    fig, ax = plt.subplots(figsize=(9, 6.5))

    coords = adata.obsm["X_umap"]
    labels = adata.obs[cluster_column].astype(str)
    unique_labels = sorted(labels.unique())
    n_clusters = len(unique_labels)

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("")

    if n_clusters <= MAX_CLUSTER_LEGEND:
        colors = _top_cluster_colors(n_clusters)

        for color, label in zip(colors, unique_labels):
            mask = labels == label
            ax.scatter(
                coords[mask, 0],
                coords[mask, 1],
                s=8,
                c=[color],
                label=label,
                linewidths=0,
            )

        legend = ax.legend(
            title=cluster_column,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False,
            markerscale=2.8,
            scatterpoints=1,
        )
        _style_legend(legend)
    else:
        if "V5_OPC_score" in adata.obs.columns:
            top_clusters = _get_top_clusters_by_v5(adata, cluster_column, TOP_CLUSTER_LABELS)
        else:
            top_clusters = unique_labels[:TOP_CLUSTER_LABELS]

        top_set = set(top_clusters)

        other_mask = ~labels.isin(top_set)
        ax.scatter(
            coords[other_mask, 0],
            coords[other_mask, 1],
            s=8,
            c="lightgrey",
            linewidths=0,
            alpha=0.9,
        )

        colors = _top_cluster_colors(len(top_clusters))
        for color, label in zip(colors, top_clusters):
            mask = labels == label
            ax.scatter(
                coords[mask, 0],
                coords[mask, 1],
                s=8,
                c=[color],
                label=label,
                linewidths=0,
            )

        legend = ax.legend(
            title=f"Top {len(top_clusters)} clusters",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            frameon=False,
            markerscale=2.8,
            scatterpoints=1,
        )
        _style_legend(legend)

    plt.tight_layout()
    return fig


def plot_umap_by_v5_score(adata):
    adata = ensure_umap(adata)

    fig, ax = plt.subplots(figsize=(9, 6.5))

    coords = adata.obsm["X_umap"]
    scores = adata.obs["V5_OPC_score"].astype(float).values

    sc_plot = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=scores,
        s=8,
        linewidths=0,
    )

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("")

    cbar = plt.colorbar(sc_plot, ax=ax)
    cbar.set_label("V5 score")

    plt.tight_layout()
    return fig


def plot_umap_highlight_top_clusters(adata, cluster_column):
    adata = ensure_umap(adata)

    cluster_means = (
        adata.obs.groupby(cluster_column)["V5_OPC_score"]
        .mean()
        .sort_values(ascending=False)
    )
    top_cluster = str(cluster_means.index[0])

    labels = adata.obs[cluster_column].astype(str)
    coords = adata.obsm["X_umap"]
    highlight = labels == top_cluster

    fig, ax = plt.subplots(figsize=(9, 6.5))

    ax.scatter(
        coords[~highlight, 0],
        coords[~highlight, 1],
        s=8,
        c="lightgrey",
        linewidths=0,
        label="Other clusters",
    )

    ax.scatter(
        coords[highlight, 0],
        coords[highlight, 1],
        s=12,
        c="red",
        linewidths=0,
        label=f"Cluster {top_cluster}",
    )

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("")
    legend = ax.legend(frameon=False, markerscale=2.8, scatterpoints=1)
    _style_legend(legend)

    plt.tight_layout()
    return fig