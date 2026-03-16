import pandas as pd


def merge_metadata(adata, metadata_file):
    """
    Merge external metadata CSV into AnnData object.
    Requires a column named 'cell_barcode' in metadata.
    """

    df = pd.read_csv(metadata_file)

    if "cell_barcode" not in df.columns:
        raise ValueError("Metadata file must contain 'cell_barcode' column.")

    df = df.set_index("cell_barcode")

    # Align with AnnData index
    adata.obs = adata.obs.join(df, how="left")

    return adata