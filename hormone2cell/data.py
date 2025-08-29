import importlib.resources
import pandas as pd
import scanpy as sc


def load_hormone_producing_file():
    """
    Load a hormone receptor data file (pickle format) packaged within the current module.
    """
    # Access the resource file within the current package and open it in binary mode
    with importlib.resources.files(__package__).joinpath('IHD_HCA_freeze6.4_hormones_20250722.pkl').open("rb") as f:
        # Use pandas to load the pickled object
        return pd.read_pickle(f)

        
def load_hormone_receptor_file():
    """
    Load a hormone receptor data file (pickle format) packaged within the current module.
    """
    # Access the resource file within the current package and open it in binary mode
    with importlib.resources.files(__package__).joinpath('IHD_HCA_freeze6.4_receptors_20250212.pkl').open("rb") as f:
        # Use pandas to load the pickled object
        return pd.read_pickle(f)

def load_hormone_file(): 
    """
    Load a hormone data file (pickle format) packaged within the current module that contains all the hormone genes.
    """
    # Access the resource file within the current package and open it in binary mode
    with importlib.resources.files(__package__).joinpath('Hormone_info_list_freeze6.4_all.pkl').open("rb") as f:
        # Use pandas to load the pickled object
        return pd.read_pickle(f)

def load_precomputed_maxvalue(assay: str) -> pd.DataFrame:
    """
    Load precomputed max average expression values for hormones.

    Parameters
    ----------
    assay : str
        Either "cell" or "nucleus".

    Returns
    -------
    pd.DataFrame
        DataFrame loaded from the corresponding pickle file.
    """
    if assay == 'cell':
        file = 'HormoneCellAtlas_v3_max_value_cell.pkl'
    elif assay == 'nucleus':
        file = 'HormoneCellAtlas_v3_max_value_nucleus.pkl'  # 注意这里和 cell 对称
    else:
        raise ValueError("assay must be either 'cell' or 'nucleus'.")

    with importlib.resources.files(__package__).joinpath(file).open("rb") as f:
        dt = pd.read_pickle(f)

    return dt



def load_pancreas_data(): 
    """
    Load a sampled pancreas data as the query dataset, which includes both single-cell and single-nucleus data. .
    """
    # Get the file path within the current package
    file_path = importlib.resources.files(__package__).joinpath('pancreas_downsample200CT.h5ad')
    # Use scanpy to read directly from the path
    return sc.read(file_path)
