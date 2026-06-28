from __future__ import annotations
from .data import load_hormone_file
import os
import gc
from typing import List

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.sparse import csr_matrix


def get_exp_percentage(
    sub: AnnData,
    gene_ids: List[str],
    clusteruse: str
) -> pd.DataFrame:
    """
    Calculate, for each cluster (cell type), the average expression of each gene,
    the number of expressing cells, and the percentage of expressing cells.
    Efficiently operates on sparse matrices to avoid expensive DataFrame conversions.

    Parameters
    ----------
    sub : AnnData
        AnnData object containing the expression matrix `X` and cell metadata `obs`.
    gene_ids : List[str]
        List of gene IDs corresponding to `sub.var.index`.
    clusteruse : str
        The column name in `sub.obs` used to group cells into clusters.

    Returns
    -------
    result_df : pd.DataFrame
        A DataFrame summarizing average expression, expressed cell counts, 
        total cell counts, and expression percentage per cluster and gene.
    """
    # Ensure sub.X is a CSR sparse matrix
    obs_matrix = sub.X if isinstance(sub.X, csr_matrix) else sub.X.tocsr()
    
    # Extract cluster labels
    cluster_labels = sub.obs[clusteruse].values
    unique_clusters = np.unique(cluster_labels)
    
    # Initialize result containers
    average_obs_list = []      # Stores mean expression values
    expressed_cell_list = []   # Stores counts of expressing cells (> 0)
    total_cell_list = []       # Stores total cell counts per cluster

    # Iterate through each cluster and compute statistics
    for cluster in unique_clusters:
        mask = cluster_labels == cluster       # Boolean mask for current cluster
        cluster_matrix = obs_matrix[mask]      # Subset sparse matrix rows

        # 1. Compute mean expression
        cluster_mean = cluster_matrix.mean(axis=0)
        average_obs_list.append(np.array(cluster_mean).flatten())

        # 2. Count cells with expression > 0
        cluster_expressed = (cluster_matrix > 0).sum(axis=0)
        expressed_cell_list.append(np.array(cluster_expressed).flatten())

        # 3. Get total number of cells in this cluster
        total_cells = cluster_matrix.shape[0]
        total_cell_list.append(total_cells)
    
    # Convert results into DataFrames
    average_obs_df = pd.DataFrame(np.vstack(average_obs_list), columns=gene_ids, index=unique_clusters)
    expressed_cell_df = pd.DataFrame(np.vstack(expressed_cell_list), columns=gene_ids, index=unique_clusters)
    total_cell_df = pd.DataFrame(total_cell_list, columns=['TotalCellNumber'], index=unique_clusters)
    
    # Reshape to long format
    average_obs_melt = average_obs_df.reset_index().melt(id_vars='index', var_name='Gene', value_name='Expression')
    average_obs_melt.rename(columns={'index': clusteruse}, inplace=True)
    
    expressed_cell_melt = expressed_cell_df.reset_index().melt(id_vars='index', var_name='Gene', value_name='ExpressedCellNumber')
    expressed_cell_melt.rename(columns={'index': clusteruse}, inplace=True)
    
    # Merge all results
    result_df = pd.merge(average_obs_melt, expressed_cell_melt, on=[clusteruse, 'Gene'])
    result_df = pd.merge(result_df, total_cell_df, left_on=clusteruse, right_index=True)
    
    # Compute percentage of expressing cells
    result_df['Percentage'] = (result_df['ExpressedCellNumber'] / result_df['TotalCellNumber']) * 100
    
    return result_df



# Define a function to scale values in column A to the range 0-1, similar with scanpy function
# def scale_within_group(x):
#     """Scale values between 0 and 1 within a group."""
#     return (x - x.min()) / (x.max() - x.min()) if x.max() != x.min() else x

# calculate the mean expression and pct of each gene in each cell type.
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData


# calculate the mean expression and pct of each gene in each cell type.
def get_result_df(sub: AnnData, celltype_col: str, tissue_col: str) -> pd.DataFrame:
    """
    Calculate average gene expression and percentage per cell type, supporting multiple tissues.

    Parameters
    ----------
    sub : AnnData
        AnnData object containing expression matrix `X`, var (genes), and obs (cell metadata).
    celltype_col : str
        Column in `obs` that defines cell types (will be converted to categorical).
    tissue_col : str
        Column in `obs` that indicates tissue identity.

    Returns
    -------
    pd.DataFrame
        Aggregated expression statistics per (tissue, cell type, gene).
    """
    # 1) Filter out genes expressed in fewer than 3 cells
    sc.pp.filter_genes(sub, min_cells=3)

    # 2) Prepare gene IDs and cell type categories
    gene_ids = list(sub.var.index.values)
    sub.obs[celltype_col] = sub.obs[celltype_col].astype("category")

    # 3) Core stats from your existing helper
    result_df = get_exp_percentage(sub, gene_ids, clusteruse=celltype_col)

    # 4) Attach tissue per cell type via a de-duplicated mapping from obs
    tissue_map = sub.obs[[celltype_col, tissue_col]].drop_duplicates(celltype_col)
    result_df = result_df.merge(tissue_map, on=celltype_col, how="left")

    # 5) Standardize column names (only if needed)
    # If your get_exp_percentage already returns 'Percentage', this is unnecessary.
    result_df.columns = result_df.columns.str.replace("percentage", "Percentage", regex=False)

    # 6) Reorder/select columns (keep only those that exist)
    keep_cols = [
        tissue_col,                      # dynamic tissue column
        celltype_col,                    # dynamic celltype column
        "Gene",
        "ExpressedCellNumber",
        "TotalCellNumber",
        "Percentage",
        "Expression",
    ]
    result_df = result_df[[c for c in keep_cols if c in result_df.columns]]

    # 7) Log-transform average expression
    result_df["logExpression"] = np.log1p(result_df["Expression"])

    # Optional: a concise status message
    print(f"Tissues in input: {sub.obs[tissue_col].unique().tolist()} — average calculation done.")

    return result_df


# ## list the files and sort by file size
# def get_sorted_file_names(folder_path):
#     files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
#     files_sorted = sorted(files, key=lambda f: os.path.getsize(os.path.join(folder_path, f)), reverse=False)
#     return files_sorted


def _ensure_csr_x(adata: AnnData) -> AnnData:
    """Ensure that `adata.X` is a CSR sparse matrix for efficient operations."""
    if isinstance(adata.X, csr_matrix):
        return adata
    adata = adata.copy()
    adata.X = adata.X.tocsr()
    return adata


def _map_celltype_tissue_to_cluster(sub: AnnData, celltype_tissue_col: str) -> AnnData:
    """
    Create a new integer-encoded `Cluster` column based on `Celltype_tissue`.

    Uses pandas.Categorical codes (0..K-1) to represent unique categories.
    """
    sub = sub.copy()
    cats = pd.Categorical(sub.obs[celltype_tissue_col])
    sub.obs["Cluster"] = pd.Index(cats.codes, dtype="int64")
    return sub


def check_count_data(adata: AnnData) -> AnnData:
    """
    Check whether adata.X or adata.raw.X contains raw counts.
    Rules:
      - If max(adata.X.data) > 100 -> return a copy of adata
      - Else, if adata.raw exists and max(adata.raw.X.data) > 100 -> return adata.raw.to_adata()
      - Else -> print message in English and raise ValueError
    """
    # 1. Check adata.X (assumed sparse matrix)
    max_x = adata.X.data.max() if adata.X.data.size > 0 else 0.0
    if max_x > 100:
        return adata

    # 2. Check adata.raw if available
    if adata.raw is not None:
        max_raw = adata.raw.X.data.max() if adata.raw.X.data.size > 0 else 0.0
        if max_raw > 100:
            return adata.raw.to_adata()

    # 3. Neither X nor raw look like raw counts
    msg = (
        "Neither adata.X nor adata.raw.X seems to contain raw counts. "
        "Hormone2cell requires raw counts as input to compute celltype-level average expression."
    )
    print(msg)
    raise ValueError(msg)


def compute_aveExp_by_category(
    adata: AnnData,
    sc_sn_col: str = "suspension_type",
    celltype_col: str = "Celltype",
    tissue_col: str ='Tissue',
) -> pd.DataFrame:
    """
    Integrated workflow for average expression analysis by category (`cell` / `nucleus`).

    Steps
    -----
    - Check count data (raw counts required)
    - Perform QC: filter cells with fewer than `min_genes`
    - Normalize counts with `normalize_total(target_sum)`
    - Digitize `Celltype_tissue` into an integer-encoded `Cluster` column
    - For each category in `sc_sn_col` (e.g. 'cell', 'nucleus'), run get_result_df()
    - Return a list of result DataFrames, one per category
    """
    # Check input data type (must be raw counts)
    adata = check_count_data(adata)
    adata = _ensure_csr_x(adata)
    ## add a tissue level unique columns
    adata.obs['Celltype_tissue']=adata.obs[tissue_col].astype('str')+'____'+adata.obs[celltype_col].astype('str')
    #celltype_tissue_col='Celltype_tissue'
    adata = _map_celltype_tissue_to_cluster(adata, 'Celltype_tissue') 

    # Quality control and normalization
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.normalize_total(adata, target_sum=10000)

    # subset to hormone related genes
    gene_dt=load_hormone_file()
    genes_use=gene_dt['Gene'].unique().tolist()
    genes_use=list(set(genes_use).intersection(set(adata.var.index)))
    if len(genes_use) > 0:
        adata = adata[:, genes_use]
        print(f'{len(genes_use)} hormone-related genes are found.')
    else:
        raise ValueError("No hormone-related genes were found in the current AnnData object.")



    # Collect categories (e.g. ['cell'] or ['nucleus'])
    categories = list(pd.unique(adata.obs[sc_sn_col].astype(str)))
    results = []

    # Case 1: only one assay
    if len(categories) == 1:
        cat=categories[0]
        print(f'Computing average gene expression per cell type in the {cat} dataset.')
        result_df = get_result_df(adata,tissue_col=tissue_col, celltype_col="Cluster")
        

        # Merge Celltype_tissue annotation directly
        result_df = result_df.merge(
            adata.obs[['Celltype_tissue', "Cluster"]].drop_duplicates("Cluster"),
            on="Cluster",
            how="left"
        )
        result_df[celltype_col] = [i.split('____')[1] for i in result_df['Celltype_tissue'].values]
        # Keep consistent column order (only if they exist)
        #result_df.columns = result_df.columns.str.replace(celltype_tissue_col, "Celltype", regex=False)
        keep_cols = [
            "Tissue",celltype_col ,"Cluster", 'Celltype_tissue', "Gene",
            "ExpressedCellNumber", "TotalCellNumber", "Percentage", "logExpression"
        ]
        result_df = result_df[[c for c in keep_cols if c in result_df.columns]]
        

        # Add category column
        result_df[sc_sn_col] = cat
        results.append(result_df)
        
        res_df = pd.concat(results, ignore_index=True)
        mask=res_df['logExpression']>0
        res_df=res_df.loc[mask]
        return res_df
        
    # Case 2: multiple assays (e.g. 'cell' and 'nucleus')
    elif len(categories) > 1:        
        for cat in categories:
            print(f'Computing average gene expression per cell type in the {cat} dataset.')
            mask = adata.obs[sc_sn_col].astype(str) == str(cat)
            sub = adata[mask].copy()
            result_df = get_result_df(sub=sub,tissue_col=tissue_col, celltype_col="Cluster")
    
            # Merge Celltype_tissue annotation directly
            result_df = result_df.merge(
                adata.obs[['Celltype_tissue', "Cluster"]].drop_duplicates("Cluster"),
                on="Cluster",
                how="left"
            )
            result_df[celltype_col] = [i.split('____')[1] for i in result_df['Celltype_tissue'].values]
            keep_cols = [
                "Tissue",celltype_col ,"Cluster", 'Celltype_tissue', "Gene",
                "ExpressedCellNumber", "TotalCellNumber", "Percentage", "logExpression"
            ]
            result_df = result_df[[c for c in keep_cols if c in result_df.columns]]
            
    
            # Add category column
            result_df[sc_sn_col] = cat
            results.append(result_df)
            
        res_df = pd.concat(results, ignore_index=True)
        mask=res_df['logExpression']>0
        res_df=res_df.loc[mask]
        return res_df



