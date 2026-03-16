import numpy as np
from tqdm import trange

OPCODE_VERSION = "1.1.0"

def score_v5_opc(adata):

    print("\n==============================")
    print("V5 OPC SCORING ENGINE")
    print("==============================")

    # --------------------------------------------------
    # DETECT GENE FORMAT
    # --------------------------------------------------

    if "gene_symbol" in adata.var.columns:
        gene_symbols = adata.var["gene_symbol"].values
        print("Using gene_symbol column.")
    elif "gene_name" in adata.var.columns:
        gene_symbols = adata.var["gene_name"].values
        print("Using gene_name column.")
    else:
        gene_symbols = adata.var_names.values
        print("Using var_names.")

    # --------------------------------------------------
    # DEFINE V5 GENE SETS
    # --------------------------------------------------

    opc_genes = ["Pdgfra", "Cspg4", "Ptprz1", "Sox10", "Olig1", "Olig2"]
    oligo_genes = ["Mbp", "Plp1", "Mag", "Mog", "Opalin", "Enpp6"]
    neuron_genes = ["Snap25", "Rbfox3", "Syt1", "Tubb3", "Map2", "Slc17a7"]
    immune_genes = ["Ptprc", "Cd68", "Aif1", "Lyz2", "Cx3cr1"]

    opc_found = [g for g in opc_genes if g in gene_symbols]
    oligo_found = [g for g in oligo_genes if g in gene_symbols]
    neuron_found = [g for g in neuron_genes if g in gene_symbols]
    immune_found = [g for g in immune_genes if g in gene_symbols]

    print(f"OPC genes found: {len(opc_found)}")
    print(f"Oligo genes found: {len(oligo_found)}")
    print(f"Neuron genes found: {len(neuron_found)}")
    print(f"Immune genes found: {len(immune_found)}")

    if len(opc_found) == 0 and len(oligo_found) == 0:
        raise ValueError("No OPC/Oligo genes detected in dataset.")

    # --------------------------------------------------
    # GET INDICES
    # --------------------------------------------------

    def get_indices(gene_list):
        return [np.where(gene_symbols == g)[0][0] for g in gene_list]

    opc_idx = get_indices(opc_found)
    oligo_idx = get_indices(oligo_found)
    neuron_idx = get_indices(neuron_found)
    immune_idx = get_indices(immune_found)

    # --------------------------------------------------
    # CHUNKED SCORING
    # --------------------------------------------------

    scores = np.zeros(adata.n_obs)
    chunk_size = 50000

    print("\nScoring cells (chunked)...")

    for start in trange(0, adata.n_obs, chunk_size):
        end = min(start + chunk_size, adata.n_obs)

        X_chunk = adata.X[start:end, :]

        if not isinstance(X_chunk, np.ndarray):
            X_chunk = X_chunk.toarray()

        opc_score = X_chunk[:, opc_idx].mean(axis=1) if opc_idx else 0
        oligo_score = X_chunk[:, oligo_idx].mean(axis=1) if oligo_idx else 0
        neuron_score = X_chunk[:, neuron_idx].mean(axis=1) if neuron_idx else 0
        immune_score = X_chunk[:, immune_idx].mean(axis=1) if immune_idx else 0

        scores[start:end] = (
            opc_score
            + oligo_score
            - neuron_score
            - immune_score
        )

    adata.obs["V5_OPC_score"] = scores

    print("Scoring complete.")
    print("Non-zero cells:", (scores != 0).sum())
    print("Max score:", scores.max())

    return adata
