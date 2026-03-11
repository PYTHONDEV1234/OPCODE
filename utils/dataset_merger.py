def merge_datasets(adata_list, dataset_names):
    """
    Merge multiple AnnData objects.
    Adds 'dataset' column to obs.
    """
    if len(adata_list) == 1:
        return adata_list[0]

    combined = adata_list[0].concatenate(
        adata_list[1:],
        batch_key="dataset",
        batch_categories=dataset_names
    )

    return combined