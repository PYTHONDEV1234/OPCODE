import argparse
import os
import scanpy as sc
import numpy as np
from tqdm import trange

from scoring_engine.v5_opc_scoring_engine import score_v5_opc


# -------------------------------------------------------
# AUTO DATASET SIZE DETECTION
# -------------------------------------------------------

def should_use_backed(input_path, size_threshold_gb=2):
    file_size_gb = os.path.getsize(input_path) / (1024 ** 3)
    return file_size_gb >= size_threshold_gb


# -------------------------------------------------------
# MAIN
# -------------------------------------------------------

def main():

    parser = argparse.ArgumentParser(description="OPC Detection Tool — V5 CLI")

    parser.add_argument("--input", required=True, help="Path to .h5ad file")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--large", action="store_true",
                        help="Force backed chunk mode")

    args = parser.parse_args()

    input_path = args.input
    output_dir = args.output
    force_large = args.large

    print("\n====================================")
    print("OPC Detection Tool — V5 CLI")
    print("====================================\n")

    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input file not found: {input_path}")

    os.makedirs(output_dir, exist_ok=True)

    # ---------------------------------------------------
    # Decide memory mode
    # ---------------------------------------------------

    use_backed = force_large or should_use_backed(input_path)

    if use_backed:
        print("Detected large dataset → Using backed mode.")
        adata = sc.read_h5ad(input_path, backed="r")
    else:
        print("Dataset small → Using in-memory mode.")
        adata = sc.read_h5ad(input_path)

    print("\nDataset summary:")
    print(adata)

    # ---------------------------------------------------
    # Scoring
    # ---------------------------------------------------

    print("\nRunning V5 scoring...")
    adata = score_v5_opc(adata)

    # ---------------------------------------------------
    # Save
    # ---------------------------------------------------

    output_path = os.path.join(output_dir, "v5_scored_output.csv")

    print("\nSaving results...")
    adata.obs.to_csv(output_path, index=False)

    print("\n====================================")
    print("Scoring Complete")
    print("====================================")
    print("Total cells:", adata.n_obs)
    print("Mean score:", float(np.mean(adata.obs["V5_OPC_score"])))
    print("Max score:", float(np.max(adata.obs["V5_OPC_score"])))
    print("Output saved to:", output_path)


if __name__ == "__main__":
    main()