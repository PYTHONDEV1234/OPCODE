import os
import scanpy as sc
import pandas as pd
from tqdm import tqdm

from scoring_engine.v5_opc_scoring_engine import score_v5_opc
from metadata.metadata_utils import merge_metadata
from analysis.effect_size import cohens_d


def score_folder(folder_path, metadata_path=None, output_dir="outputs"):

    os.makedirs(output_dir, exist_ok=True)

    summary_rows = []

    files = [f for f in os.listdir(folder_path) if f.endswith(".h5ad")]

    for file in files:

        file_path = os.path.join(folder_path, file)
        print(f"\nProcessing: {file}")

        adata = sc.read_h5ad(file_path)

        # Score V5
        adata = score_v5_opc(adata)

        # Merge metadata if provided
        if metadata_path:
            adata = merge_metadata(adata, metadata_path)

        # Save per-file scores
        output_file = os.path.join(
            output_dir,
            f"v5_{file.replace('.h5ad','')}_scores.csv"
        )

        adata.obs.to_csv(output_file, index=False)

        # Compute summary metrics if class column exists
        if "class" in adata.obs.columns:

            oligo = adata.obs[adata.obs["class"].str.contains("Oligo", case=False, na=False)]["V5_OPC_score"]
            glut = adata.obs[adata.obs["class"].str.contains("Glut", case=False, na=False)]["V5_OPC_score"]
            gaba = adata.obs[adata.obs["class"].str.contains("GABA", case=False, na=False)]["V5_OPC_score"]

            if len(oligo) > 0 and len(glut) > 0:
                d_og = cohens_d(oligo, glut)
            else:
                d_og = None

            if len(oligo) > 0 and len(gaba) > 0:
                d_ob = cohens_d(oligo, gaba)
            else:
                d_ob = None

            summary_rows.append({
                "Region": file,
                "Oligo_Mean": oligo.mean() if len(oligo) > 0 else None,
                "Glut_Mean": glut.mean() if len(glut) > 0 else None,
                "GABA_Mean": gaba.mean() if len(gaba) > 0 else None,
                "Effect_Oligo_vs_Glut": d_og,
                "Effect_Oligo_vs_GABA": d_ob,
                "Dynamic_Range": adata.obs["V5_OPC_score"].max() - adata.obs["V5_OPC_score"].min()
            })

    # Save summary table
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_path = os.path.join(output_dir, "v5_batch_summary.csv")
        summary_df.to_csv(summary_path, index=False)
        print(f"\nBatch summary saved to {summary_path}")

    print("\nBatch scoring complete.")