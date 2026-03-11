import os

def find_h5ad_files(folder_path):
    """
    Scan a folder for .h5ad files.
    Returns full paths.
    """
    h5ad_files = []

    if not os.path.isdir(folder_path):
        return []

    for file in os.listdir(folder_path):
        if file.endswith(".h5ad"):
            h5ad_files.append(os.path.join(folder_path, file))

    return h5ad_files