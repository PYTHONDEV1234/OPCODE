# analysis/purification.py
"""
OPCODE — Purification utilities

This module provides:
1) add_v5_normalized_scores(adata)
   - Adds cross-dataset comparable normalization columns:
       - V5_z   : z-score of V5_OPC_score within the dataset
       - V5_0_1 : sigmoid(V5_z), bounded (0,1)

2) build_purification_outputs(...)
   - Identifies top OPC-like clusters (by mean V5_OPC_score)
   - Exports:
       - opc_candidates.csv : all cells in top-k clusters
       - opc_highconf.csv   : high-confidence subset based on V5_0_1 + optional marker gate
   - Returns a summary dict used by the CLI and report builder.

Design goals:
- Deterministic, dataset-only (no web), robust to missing columns.
- Avoid pandas groupby FutureWarnings by using observed=True.
- Produce outputs that are easy to map into downstream wet-lab gating workflows.

NOTE on anchors:
- The *surface panel builder* is responsible for writing sorting_panel.yaml,
  but we also attach a standard anchor-refinement hint here so it can be reused
  consistently by other modules/CLI if desired.
"""

from __future__ import annotations

import os
import math
from typing import Dict, Optional, List, Any

import numpy as np
import pandas as pd

try:
    import scipy.sparse as sp
except Exception:
    sp = None


ANCHOR_REFINE_TEXT = (
    "Optional refine/test: AND PDGFRα(CD140a)+ / NG2+ / ITGA6+ (anchors; validate experimentally)"
)


# ----------------------------
# Helpers
# ----------------------------

def _sigmoid(x: np.ndarray) -> np.ndarray:
    # Numerically stable sigmoid
    x = np.clip(x, -50, 50)
    return 1.0 / (1.0 + np.exp(-x))


def _ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def _detect_mito_col(obs: pd.DataFrame, user_col: Optional[str]) -> Optional[str]:
    if user_col and user_col in obs.columns:
        return user_col

    candidates = [
        "pct_counts_mito",
        "percent_mito",
        "pct_mito",
        "mito_pct",
        "mitochondrial_percent",
        "mt_percent",
    ]
    for c in candidates:
        if c in obs.columns:
            return c
    return None


def _to_dense_vector(x) -> np.ndarray:
    # Converts a vector (n,) or (n,1) to a dense 1D float array
    if sp is not None and sp.issparse(x):
        x = x.toarray()
    x = np.asarray(x)
    if x.ndim == 2 and x.shape[1] == 1:
        x = x[:, 0]
    return x.astype(float, copy=False)


def _get_gene_indices(var_names: pd.Index, genes: List[str]) -> List[int]:
    lower_map = {g.lower(): i for i, g in enumerate(var_names)}
    idx = []
    for g in genes:
        i = lower_map.get(g.lower())
        if i is not None:
            idx.append(i)
    return idx


def _mean_expr_score(adata, genes: List[str]) -> np.ndarray:
    """
    Compute a simple mean expression score over a list of genes.
    Uses adata.X (normalized space). Falls back gracefully if genes missing.
    """
    var_names = getattr(adata, "var_names", None)
    if var_names is None:
        return np.zeros(adata.n_obs, dtype=float)

    idx = _get_gene_indices(var_names, genes)
    if len(idx) == 0:
        return np.zeros(adata.n_obs, dtype=float)

    X = adata.X
    if sp is not None and sp.issparse(X):
        sub = X[:, idx]
        # mean over columns
        score = np.asarray(sub.mean(axis=1)).ravel()
    else:
        sub = np.asarray(X)[:, idx]
        score = sub.mean(axis=1)

    score = np.asarray(score).ravel().astype(float, copy=False)
    return score


def _gene_positive_mask(adata, gene: str, threshold: float = 0.0) -> np.ndarray:
    """
    Returns boolean mask where expression of `gene` > threshold.
    Uses adata.X (normalized space). If gene missing, returns all False.
    """
    var_names = getattr(adata, "var_names", None)
    if var_names is None:
        return np.zeros(adata.n_obs, dtype=bool)

    idx = _get_gene_indices(var_names, [gene])
    if not idx:
        return np.zeros(adata.n_obs, dtype=bool)

    X = adata.X
    if sp is not None and sp.issparse(X):
        vec = X[:, idx[0]]
        vec = np.asarray(vec.toarray()).ravel()
    else:
        vec = np.asarray(X)[:, idx[0]].ravel()

    return vec > threshold


# ----------------------------
# Public API
# ----------------------------

def add_v5_normalized_scores(adata) -> None:
    """
    Adds:
      - adata.obs['V5_z']
      - adata.obs['V5_0_1']

    Also stores the anchor refine hint in adata.uns for reuse.
    """
    if "V5_OPC_score" not in adata.obs.columns:
        raise ValueError("Expected adata.obs['V5_OPC_score'] to exist. Run scoring first.")

    scores = pd.to_numeric(adata.obs["V5_OPC_score"], errors="coerce").fillna(0.0).to_numpy(dtype=float)
    mu = float(np.mean(scores))
    sd = float(np.std(scores))

    if sd == 0.0 or not np.isfinite(sd):
        z = np.zeros_like(scores)
    else:
        z = (scores - mu) / sd

    adata.obs["V5_z"] = z
    adata.obs["V5_0_1"] = _sigmoid(z)

    # Store for downstream modules (surface panel builder, report)
    if not hasattr(adata, "uns") or adata.uns is None:
        adata.uns = {}
    adata.uns["opcode_anchor_refine_text"] = ANCHOR_REFINE_TEXT
    adata.uns["opcode_v5_norm_mu"] = mu
    adata.uns["opcode_v5_norm_sd"] = sd


def build_purification_outputs(
    adata,
    output_dir: str,
    cluster_col: Optional[str],
    top_k_clusters: int = 1,
    qc_on: bool = True,
    min_genes: int = 200,
    min_counts: Optional[int] = None,
    max_mito: float = 20.0,
    mito_col: Optional[str] = None,
    v5_threshold: Optional[float] = None,
    top_pct: float = 10.0,
    marker_gate_quantile: float = 0.60,
    use_marker_gate: bool = True,
    strict_mature_exclusion: bool = False,
) -> Dict[str, Any]:
    """
    Returns dict with:
      - opc_candidates_path
      - opc_highconf_path
      - n_candidates
      - n_candidates_after_qc
      - n_highconf
      - highconf_fraction_of_candidates
      - used_v5_threshold (or None)
      - used_top_pct
      - used_marker_gate
      - marker_gate_quantile
      - strict_mature_exclusion
      - mito_col_used
      - selected_clusters (list)
      - anchor_refine_text (string)
    """
    _ensure_dir(output_dir)

    if "V5_OPC_score" not in adata.obs.columns:
        raise ValueError("Expected adata.obs['V5_OPC_score'] to exist. Run scoring first.")
    if "V5_0_1" not in adata.obs.columns or "V5_z" not in adata.obs.columns:
        add_v5_normalized_scores(adata)

    obs = adata.obs.copy()

    # Ensure cell_id column for exports
    obs["cell_id"] = obs.index.astype(str)

    # ----------------------------
    # Cluster ranking / selection
    # ----------------------------
    if cluster_col is None or cluster_col not in obs.columns:
        # No clusters: treat entire dataset as "candidates"
        selected_clusters = []
        candidate_mask = np.ones(adata.n_obs, dtype=bool)
    else:
        # mean V5 per cluster
        cluster_means = (
            obs.groupby(cluster_col, observed=True)["V5_OPC_score"]
            .mean()
            .sort_values(ascending=False)
        )
        selected_clusters = cluster_means.head(max(1, int(top_k_clusters))).index.astype(str).tolist()
        candidate_mask = obs[cluster_col].astype(str).isin(selected_clusters).to_numpy()

    candidates = obs.loc[candidate_mask].copy()

    # ----------------------------
    # QC filtering (optional)
    # ----------------------------
    mito_used = _detect_mito_col(candidates, mito_col)
    qc_mask = np.ones(len(candidates), dtype=bool)

    if qc_on:
        if "n_genes" in candidates.columns:
            qc_mask &= pd.to_numeric(candidates["n_genes"], errors="coerce").fillna(0).to_numpy() >= int(min_genes)
        elif "n_genes_by_counts" in candidates.columns:
            qc_mask &= pd.to_numeric(candidates["n_genes_by_counts"], errors="coerce").fillna(0).to_numpy() >= int(min_genes)

        if min_counts is not None:
            if "total_counts" in candidates.columns:
                qc_mask &= pd.to_numeric(candidates["total_counts"], errors="coerce").fillna(0).to_numpy() >= int(min_counts)

        if mito_used is not None:
            qc_mask &= pd.to_numeric(candidates[mito_used], errors="coerce").fillna(0).to_numpy() <= float(max_mito)

    candidates_qc = candidates.loc[qc_mask].copy()

    # Optional strict mature exclusion (RNA-level; conservative)
    if strict_mature_exclusion:
        # Remove cells expressing mature OL markers (very strict)
        # Use adata.X for these genes; evaluate only candidate_qc cells
        cell_ids = candidates_qc["cell_id"].tolist()
        # Build mask in full space then subset
        mbp_pos = _gene_positive_mask(adata, "Mbp", threshold=0.0)
        plp1_pos = _gene_positive_mask(adata, "Plp1", threshold=0.0)
        mature_pos = mbp_pos | plp1_pos

        mature_lookup = pd.Series(mature_pos, index=adata.obs.index.astype(str))
        candidates_qc = candidates_qc.loc[~candidates_qc["cell_id"].map(mature_lookup).fillna(False).to_numpy()].copy()

    # ----------------------------
    # High-confidence selection by V5_0_1
    # ----------------------------
    if len(candidates_qc) == 0:
        # Still export empties with headers
        opc_candidates_path = os.path.join(output_dir, "opc_candidates.csv")
        opc_highconf_path = os.path.join(output_dir, "opc_highconf.csv")

        candidates.to_csv(opc_candidates_path, index=False)
        pd.DataFrame(columns=list(candidates.columns)).to_csv(opc_highconf_path, index=False)

        return {
            "opc_candidates_path": opc_candidates_path,
            "opc_highconf_path": opc_highconf_path,
            "n_candidates": int(len(candidates)),
            "n_candidates_after_qc": 0,
            "n_highconf": 0,
            "highconf_fraction_of_candidates": 0.0,
            "used_v5_threshold": v5_threshold,
            "used_top_pct": float(top_pct),
            "used_marker_gate": bool(use_marker_gate),
            "marker_gate_quantile": float(marker_gate_quantile),
            "strict_mature_exclusion": bool(strict_mature_exclusion),
            "mito_col_used": mito_used,
            "selected_clusters": selected_clusters,
            "anchor_refine_text": ANCHOR_REFINE_TEXT,
        }

    v5_01 = pd.to_numeric(candidates_qc["V5_0_1"], errors="coerce").fillna(0.0).to_numpy(dtype=float)

    if v5_threshold is not None:
        v5_keep = v5_01 >= float(v5_threshold)
        used_thr = float(v5_threshold)
        used_top_pct = float(top_pct)
    else:
        # Keep top X% by V5_0_1
        pct = float(top_pct)
        pct = min(max(pct, 0.1), 100.0)
        cutoff = np.quantile(v5_01, 1.0 - (pct / 100.0))
        v5_keep = v5_01 >= float(cutoff)
        used_thr = None
        used_top_pct = pct

    pre_marker = candidates_qc.loc[v5_keep].copy()

    # ----------------------------
    # Marker gate (optional): refine within the V5-filtered candidates
    # ----------------------------
    # Use a small OPC/oligo lineage RNA marker set as a "confidence" gate.
    # This is NOT a wet-lab gate; it’s computational cleanup.
    marker_genes = ["Sox10", "Olig1", "Olig2", "Nkx2-2", "Pdgfra", "Cspg4"]

    if use_marker_gate and len(pre_marker) > 0:
        # Score in full space, then map onto pre_marker rows
        score_vec = _mean_expr_score(adata, marker_genes)
        score_lookup = pd.Series(score_vec, index=adata.obs.index.astype(str))
        pre_marker["marker_score"] = pre_marker["cell_id"].map(score_lookup).astype(float)

        q = float(marker_gate_quantile)
        q = min(max(q, 0.0), 1.0)
        gate_cut = pre_marker["marker_score"].quantile(q)
        highconf = pre_marker.loc[pre_marker["marker_score"] >= gate_cut].copy()
        used_marker_gate = True
    else:
        if len(pre_marker) > 0 and "marker_score" not in pre_marker.columns:
            pre_marker["marker_score"] = np.nan
        highconf = pre_marker
        used_marker_gate = False

    # ----------------------------
    # Export
    # ----------------------------
    opc_candidates_path = os.path.join(output_dir, "opc_candidates.csv")
    opc_highconf_path = os.path.join(output_dir, "opc_highconf.csv")

    # Add cluster column for convenience (if available)
    if cluster_col is not None and cluster_col in obs.columns:
        # ensure it's present in exported frames
        if cluster_col not in candidates.columns:
            candidates[cluster_col] = obs.loc[candidates.index, cluster_col]
        if cluster_col not in candidates_qc.columns:
            candidates_qc[cluster_col] = obs.loc[candidates_qc.index, cluster_col]
        if cluster_col not in highconf.columns and len(highconf) > 0:
            highconf[cluster_col] = obs.loc[highconf.index, cluster_col]

    # Keep exports tidy (key columns first)
    key_cols = ["cell_id"]
    if cluster_col is not None and cluster_col in candidates.columns:
        key_cols.append(cluster_col)
    for c in ["V5_OPC_score", "V5_z", "V5_0_1", "marker_score"]:
        if c in candidates.columns and c not in key_cols:
            key_cols.append(c)

    def _reorder(df: pd.DataFrame) -> pd.DataFrame:
        cols = list(df.columns)
        front = [c for c in key_cols if c in cols]
        rest = [c for c in cols if c not in front]
        return df.loc[:, front + rest]

    _reorder(candidates).to_csv(opc_candidates_path, index=False)
    if len(highconf) > 0:
        _reorder(highconf).to_csv(opc_highconf_path, index=False)
    else:
        pd.DataFrame(columns=list(candidates.columns)).to_csv(opc_highconf_path, index=False)

    n_candidates = int(len(candidates))
    n_after_qc = int(len(candidates_qc))
    n_highconf = int(len(highconf))
    frac = float(n_highconf / n_after_qc) if n_after_qc > 0 else 0.0

    # Persist anchor text for consistent downstream messaging
    if not hasattr(adata, "uns") or adata.uns is None:
        adata.uns = {}
    adata.uns["opcode_anchor_refine_text"] = ANCHOR_REFINE_TEXT

    return {
        "opc_candidates_path": opc_candidates_path,
        "opc_highconf_path": opc_highconf_path,
        "n_candidates": n_candidates,
        "n_candidates_after_qc": n_after_qc,
        "n_highconf": n_highconf,
        "highconf_fraction_of_candidates": frac,
        "used_v5_threshold": used_thr,
        "used_top_pct": used_top_pct,
        "used_marker_gate": used_marker_gate,
        "marker_gate_quantile": float(marker_gate_quantile),
        "strict_mature_exclusion": bool(strict_mature_exclusion),
        "mito_col_used": mito_used,
        "selected_clusters": selected_clusters,
        "anchor_refine_text": ANCHOR_REFINE_TEXT,
    }