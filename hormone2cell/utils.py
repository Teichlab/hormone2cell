from __future__ import annotations
import pandas as pd
import numpy as np
import os
import time
from collections import Counter
from typing import Optional, Sequence, Iterable,Collection, Tuple

#from .data import load_hormone_file
from .data import load_precomputed_maxvalue


## hormone function
def AND_logic(df, and_genes):
    """Vectorized AND logic with gene exclusions"""
    sub_df = df.loc[and_genes]
    exclude_genes = ['PCSK1', 'PCSK2', 'CGA', 'TG', 'CYB5A']
    rows_to_drop = [g for g in exclude_genes if g in sub_df.index]
    
    if rows_to_drop:
        sub_df_adjusted = sub_df.drop(rows_to_drop, errors='ignore')
    else:
        sub_df_adjusted = sub_df
    
    valid_columns = (sub_df != 0).all(axis=0)
    result = pd.Series(0.0, index=sub_df.columns)
    means = sub_df_adjusted.mean(axis=0).fillna(0)
    result[valid_columns] = means[valid_columns]
    return result

def OR_logic(df, or_genes):
    """Vectorized OR logic"""
    if not or_genes:
        return pd.Series(0.0, index=df.columns)
        
    sub_df = df.loc[or_genes]
    valid_columns = (sub_df != 0).any(axis=0)
    result = pd.Series(0.0, index=sub_df.columns)
    result[valid_columns] = sub_df.mean(axis=0)[valid_columns]
    return result

## hormone rule for hormone with single gene
def assemble_hormone_expression_excludeDiff(ave_wide, 
                                            ave_wide_exclude,
                                            and_genes,
                                            exclude_genes):
    """Vectorized hormone assembly with exclusion rules"""
    if not and_genes:
        return pd.DataFrame(0, index=['hormone'], columns=ave_wide.columns)
        
    if len(and_genes) == 1:
        and_exp = ave_wide.loc[and_genes].copy()
    else:
        and_exp = pd.DataFrame(AND_logic(ave_wide, and_genes)).T

    if exclude_genes:
        valid_exclude = [g for g in exclude_genes if g in ave_wide_exclude.index]
        if valid_exclude:
            exclude_expr = pd.DataFrame(OR_logic(ave_wide_exclude, valid_exclude)).T
            and_exp = np.where(exclude_expr > 0, 0, and_exp)
            and_exp = pd.DataFrame(and_exp, columns=ave_wide.columns)
    return and_exp

## hormone rule for hormone with combined genes
def combine_expression_fair_excludeDiff(ave_wide, 
                                        ave_wide_exclude, 
                                        gene_sets,
                                        hormone_name="hormone"):
    """Optimized hormone combination with exclusion rules"""
    all_genes = [g for and_genes, _ in gene_sets for g in and_genes if g in ave_wide.index]
    gene_counts = Counter(all_genes)
    
    result = pd.Series(0.0, index=ave_wide.columns)
    
    for and_genes, exclude_genes in gene_sets:
        if not set(and_genes).issubset(ave_wide.index):
            continue
            
        and_exp = ave_wide.loc[and_genes]
        valid_cells = (and_exp > 1e-5).all(axis=0)
        
        if not valid_cells.any():
            continue
            
        if exclude_genes:
            valid_excl = [g for g in exclude_genes if g in ave_wide_exclude.index]
            if valid_excl:
                excl_exp = ave_wide_exclude.loc[valid_excl]
                exclusion_cells = (excl_exp > 1e-5).any(axis=0)
                valid_cells &= ~exclusion_cells
                
        if not valid_cells.any():
            continue
            
        weights = pd.Series({g: gene_counts.get(g, 1) for g in and_genes})
        weighted_expr = ave_wide.loc[and_genes].div(weights, axis=0).mean(axis=0)
        weighted_expr[~valid_cells] = 0
        result += weighted_expr
        
    return pd.DataFrame([result], index=[hormone_name])


def get_hormone_genes_single(row, 
                             columns_use):
    """Extract hormone genes from a single row"""
    genes = [g for g in row[columns_use].dropna().astype(str) if g != 'nan']
    return genes


def initial_filtering(ave_all: pd.DataFrame,
                      thresh_expr_low,
                      thresh_pct,
                      hormone_genes):
    """
    input: raw average expression file
    output: filtered expression file 
    """
    # filter cell type_gene pairs
    mask = (
        (ave_all['TotalCellNumber'] >= 10) &
        (ave_all['ExpressedCellNumber'] >= 3) &
        (ave_all['logExpression'] > thresh_expr_low) &
        (ave_all['Percentage'] > thresh_pct)
    )
    ave_positive = ave_all[mask].copy()

    # filter genes
    mask=ave_positive['Gene'].isin(hormone_genes)
    ave_positive = ave_positive.loc[mask]

    return ave_positive


# def calculate_gene_manner_threshold(ave_positive,
#                                     threshold_expr_gene=0.1
#                                    ):
#     #  Calculate expression thresholds
#     #print("Calculating expression thresholds...")
#     max_expression = (
#         ave_positive.groupby('Gene', observed=False)['logExpression']
#         .quantile(0.999)
#         .reset_index()
#         #.merge(hormone_genes, on='Gene')
#     )
#     max_expression.rename(columns={'logExpression': 'Max_logExp'}, inplace=True)
    
#     # Create initial thresholds
#     max_expression['Threshold'] = max_expression['Max_logExp'] * threshold_expr_gene

#     # filter target genes
#     # mask=max_expression['Gene'].isin(gene_use)
#     # max_expression=max_expression[mask]
    
#     # Return initial thresholds
#     max_expression=max_expression.set_index('Gene')#['Threshold'].to_dict()
#     return max_expression


# def calculate_gene_manner_threshold(ave_positive,
#                                     threshold_expr_gene=0.1,
#                                     use_precomputed:Optional[str] = None
#                                    ):
#     #  Calculate expression thresholds
#     #print("Calculating expression thresholds...")
#     if use_precomputed is None:
#         max_expression = (
#             ave_positive.groupby('Gene', observed=False)['logExpression']
#             .quantile(0.999)
#             .reset_index()
#         )
#         max_expression.rename(columns={'logExpression': 'Max_logExp'}, inplace=True)

#     elif use_precomputed == 'cell':
#         max_expression = load_precomputed_maxvalue(assay=use_precomputed)

#     elif use_precomputed == 'nucleus':
#         max_expression = load_precomputed_maxvalue(assay=use_precomputed)

#     else:
#         raise ValueError(
#             "Invalid value for use_precomputed. "
#             "It must be one of: 'cell', 'nucleus', or None."
#         )
    
#     # Create initial thresholds
#     max_expression['Threshold'] = max_expression['Max_logExp'] * threshold_expr_gene

#     # filter target genes
#     # mask=max_expression['Gene'].isin(gene_use)
#     # max_expression=max_expression[mask]
    
#     # Return initial thresholds
#     max_expression=max_expression.set_index('Gene')#['Threshold'].to_dict()
#     return max_expression

def calculate_gene_manner_threshold(
    ave_positive: pd.DataFrame,
    threshold_expr_gene: float = 0.1,
    use_precomputed: Optional[str] = None
) -> pd.DataFrame:
    """
    Calculate per-gene thresholds. If `use_precomputed` is provided ("cell"/"nucleus"),
    load precomputed Max_logExp; for genes absent from the precomputed table,
    fall back to computing from `ave_positive` (0.999 quantile of logExpression).
    """
    # Basic input checks (optional but helpful)
    required_cols = {"Gene", "logExpression"}
    if not required_cols.issubset(ave_positive.columns):
        missing = required_cols - set(ave_positive.columns)
        raise ValueError(f"`ave_positive` is missing required columns: {missing}")

    genes_all = ave_positive["Gene"].unique()

    if use_precomputed is None:
        # Compute Max_logExp from current data
        print('            Use Max_logExp from current data')
        max_expression = (
            ave_positive.groupby("Gene", observed=False)["logExpression"]
            .quantile(0.999)
            .reset_index()
            .rename(columns={"logExpression": "Max_logExp"})
        )

    elif use_precomputed in {"cell", "nucleus"}:
        # Load precomputed table
        print('            Use Max_logExp from precomputed data')
        pre_df = load_precomputed_maxvalue(assay=use_precomputed)

        # Normalize shape: ensure it has 'Gene' and 'Max_logExp'
        if "Gene" not in pre_df.columns:
            # assume index are gene ids
            pre_df = pre_df.reset_index().rename(columns={"index": "Gene"})
        if "Max_logExp" not in pre_df.columns:
            raise ValueError("Precomputed table must contain a 'Max_logExp' column.")

        # Keep only genes present in current data
        pre_df = pre_df[pre_df["Gene"].isin(genes_all)].copy()

        # Identify genes missing from the precomputed table
        missing_genes = np.setdiff1d(genes_all, pre_df["Gene"].unique())

        if missing_genes.size > 0:
            print(
                    f"[Warning] {len(missing_genes)} genes are missing in the "
                    f"precomputed '{use_precomputed}' dataset. "
                    f"Falling back to local calculation for these genes:\n"
                    f"{', '.join(missing_genes)}"
                )
            # Compute fallback Max_logExp for missing genes using current data
            fallback = (
                ave_positive[ave_positive["Gene"].isin(missing_genes)]
                .groupby("Gene", observed=False)["logExpression"]
                .quantile(0.999)
                .reset_index()
                .rename(columns={"logExpression": "Max_logExp"})
            )
            # Combine (precomputed first, then fallback), and drop any accidental dups
            max_expression = (
                pd.concat([pre_df[["Gene", "Max_logExp"]], fallback], ignore_index=True)
                .drop_duplicates(subset="Gene", keep="first")
            )
        else:
            max_expression = pre_df[["Gene", "Max_logExp"]].copy()
    
    else:
        raise ValueError(
            "Invalid value for use_precomputed. It must be one of: 'cell', 'nucleus', or None."
        )
    #print(max_expression)
    # Calculate thresholds
    max_expression["Threshold"] = max_expression["Max_logExp"] * threshold_expr_gene

    # Return indexed by Gene 
    return max_expression.set_index("Gene")


def create_expression_matrix(ave_positive,max_expression, celltype_column):
    """
    cast the long table into a wide dataframe for rule calculation
    """
    #print("Creating expression matrices...")
    ave_wide = ave_positive.pivot(index='Gene', columns=celltype_column, values='logExpression')
    ave_wide[pd.isna(ave_wide)]=0    
    genes=list(ave_wide.index)
    max_expression1=max_expression.loc[genes,:]
    ave_wide.loc[genes] = np.where(
        ave_wide.loc[genes] > max_expression1.loc[genes, 'Threshold'].values[:, None],
        ave_wide.loc[genes],
        0)
    return(ave_wide)



def adopt_hormone_rules(
    hormone_producing, 
    ave_wide, 
    ave_wide_exclude, 
    include_cols,
    exclude_cols
):
    """Hormone processing with AND and OR logic"""
    hormone_short = hormone_producing['hormone_short'].unique()
    holder_first = []
    gene_sets_map = {}
    
    # Precompute gene sets for all hormones
    for hormone in hormone_short:
        hormone_rows = hormone_producing[hormone_producing['hormone_short'] == hormone]
        hormone_sets = []
        
        for _, row in hormone_rows.iterrows():
            and_genes = get_hormone_genes_single(row, include_cols)
            exclude_genes = get_hormone_genes_single(row, exclude_cols)
            
            if set(and_genes).issubset(ave_wide.index):
                hormone_sets.append((and_genes, exclude_genes))
        
        gene_sets_map[hormone] = hormone_sets
    
    # First processing pass
    for hormone, gene_sets in gene_sets_map.items():
        if not gene_sets:
            continue
            
        if len(gene_sets) == 1:
            dt = assemble_hormone_expression_excludeDiff(
                ave_wide, ave_wide_exclude, gene_sets[0][0], gene_sets[0][1]
            )
            dt.index = [hormone]
            holder_first.append(dt)
        else:
            result = combine_expression_fair_excludeDiff(
                ave_wide, ave_wide_exclude, gene_sets, hormone_name=hormone
            )
            holder_first.append(result)
    
    first_pass_result = pd.concat(holder_first) if holder_first else pd.DataFrame()
    return first_pass_result


from typing import Optional, Set
import pandas as pd
import numpy as np

# def get_specificity_adjustment_hormones(
#     first_pass_result: pd.DataFrame,
#     specificity_threshold: Optional[float] = 0.8,
#     coverage_threshold: Optional[float] = 0.5,
#     cell_annos: Optional[pd.DataFrame] = None,
#     tissue_col: str = "Tissue",
#     celltype_col: str = "Celltype_unique",
# ) -> Optional[Set[str]]:
#     """
#     Identify hormones that need adjustment based on (1) specificity (e.g., τ index) and
#     (2) broad coverage across cell types within tissues.

#     Only if `cell_annos` is provided (not None) will the coverage criterion be evaluated.
#     If `cell_annos` is None, coverage-based adjustment is skipped entirely.

#     Parameters
#     ----------
#     first_pass_result : pd.DataFrame
#         Wide-format matrix with hormones as rows and cell types as columns. Values are expression.
#     specificity_threshold : Optional[float], default 0.8
#         If not None, hormones with specificity (τ) < threshold will be flagged.
#     coverage_threshold : Optional[float], default 0.5
#         If not None and `cell_annos` is provided, hormones expressed in > (threshold * #celltypes)
#         within any tissue will be flagged.
#     cell_annos : Optional[pd.DataFrame], default None
#         Cell annotations; must contain `tissue_col` and `celltype_col` if provided.
#     tissue_col : str, default "Tissue"
#         Column in `cell_annos` indicating tissue labels. Used only when `cell_annos` is not None.
#     celltype_col : str, default "Celltype_unique"
#         Column in `cell_annos` indicating cell type labels. Used only when `cell_annos` is not None.

#     Returns
#     -------
#     Optional[Set[str]]
#         Set of hormones to adjust, or None if empty.
#     """
#     hormones_to_adjust: Set[str] = set()

#     # ---- SPECIFICITY PHASE (always considered if threshold is given) ----
#     if specificity_threshold is not None:
#         hormone_specificity = {}
#         for hormone in first_pass_result.index:
#             # ensure 1D array for fTau
#             vals = np.asarray(first_pass_result.loc[hormone].values).ravel()
#             spec = fTau(vals)  
#             hormone_specificity[hormone] = spec

#         hormones_to_adjust.update(
#             h for h, spec in hormone_specificity.items()
#             if spec is not None and spec < specificity_threshold and spec != 0
#         )
#         print(f'The number of hormones adjusted because of the specifity: {len(hormones_to_adjust)}')
#         print(f"hormones_to_adjust: {hormones_to_adjust}")
#         print(f"specifity score: {hormone_specificity}")
#     # ---- COVERAGE PHASE (only if cell_annos is provided) ----
#     if (coverage_threshold is not None) and (cell_annos is not None):
#         # basic validation
#         missing_cols = [c for c in (tissue_col, celltype_col) if c not in cell_annos.columns]
#         if missing_cols:
#             raise ValueError(
#                 f"`cell_annos` is missing required columns: {missing_cols}. "
#                 "Provide the correct column names or set `cell_annos=None` to skip coverage filtering."
#             )

#         # number of unique cell types per tissue
#         # (use unique on the subset to avoid counting duplicates)
#         tissue_groups = cell_annos[[tissue_col, celltype_col]].dropna().drop_duplicates()
#         n_celltypes_per_tissue = tissue_groups.groupby(tissue_col)[celltype_col].nunique().to_dict()

#         for tissue, n_celltypes in n_celltypes_per_tissue.items():
#             # all celltypes for this tissue (as they should appear as columns in first_pass_result)
#             tissue_celltypes = tissue_groups.loc[
#                 tissue_groups[tissue_col] == tissue, celltype_col
#             ].tolist()

#             # keep only cell types that exist in first_pass_result columns
#             tissue_celltypes = [ct for ct in tissue_celltypes if ct in first_pass_result.columns]
#             if not tissue_celltypes:
#                 continue

#             #threshold_count = coverage_threshold * n_celltypes

#             # for each hormone, count how many celltypes in this tissue have expression > 0
#             expr_block = first_pass_result.loc[:, tissue_celltypes]
#             # boolean > 0 then sum across columns
#             for hormone, row in expr_block.iterrows():
#                 expressed_count = (row.values > 0).sum()
#                 #if expressed_count >= threshold_count:
#                 if expressed_count/n_celltypes > coverage_threshold:
#                     hormones_to_adjust.add(hormone)

#     print(f"hormones_to_adjust: {hormones_to_adjust}")
#     return hormones_to_adjust if len(hormones_to_adjust) > 0 else None



from typing import Optional, Dict, List
import pandas as pd
import numpy as np

def get_hormones_by_specificity(
    first_pass_result: pd.DataFrame,
    cell_annos: Optional[pd.DataFrame] = None,
    specificity_threshold: Optional[float] = 0.8,
    tissue_col: str = "Tissue",
    hormone_col: str = "Hormone"
) -> Dict[str, List[str]]:
    """
    Identify hormones to adjust based on specificity (τ index).
    Returns dict: {hormone -> list of tissues}.
    """
    hormone_to_tissues: Dict[str, List[str]] = {}
    if specificity_threshold is None or cell_annos is None:
        return hormone_to_tissues

    # calculate τ index for each hormone
    hormone_specificity = {}
    for hormone in first_pass_result.index:
        vals = np.asarray(first_pass_result.loc[hormone].values).ravel()
        spec = fTau(vals)  # assume fTau is defined elsewhere
        hormone_specificity[hormone] = spec

    # check which hormones fall below the specificity threshold
    for h, spec in hormone_specificity.items():
        if spec is not None and spec < specificity_threshold and spec != 0:
            # get tissues where this hormone appears in cell_annos
            #mask = cell_annos[hormone_col] == h
            tissues = cell_annos[tissue_col].unique().tolist()
            if tissues:
                hormone_to_tissues[h] = tissues

    return hormone_to_tissues


from typing import Optional, Dict, List
import pandas as pd

def get_hormones_by_coverage(
    first_pass_result: pd.DataFrame,
    cell_annos: Optional[pd.DataFrame],
    coverage_threshold: Optional[float] = 0.5,
    tissue_col: str = "Tissue",
    celltype_col: str = "Celltype_unique",
) -> Dict[str, List[str]]:
    """
    Identify hormones to adjust based on broad coverage across cell types in tissues.
    Returns a dict: {hormone -> list of tissues}.
    """
    hormone_to_tissues: Dict[str, List[str]] = {}
    if coverage_threshold is None or cell_annos is None:
        return hormone_to_tissues

    # group tissues and count number of unique cell types
    tissue_groups = cell_annos[[tissue_col, celltype_col]].dropna().drop_duplicates()
    n_celltypes_per_tissue = tissue_groups.groupby(tissue_col)[celltype_col].nunique().to_dict()

    for tissue, n_celltypes in n_celltypes_per_tissue.items():
        tissue_celltypes = tissue_groups.loc[
            tissue_groups[tissue_col] == tissue, celltype_col
        ].tolist()
        tissue_celltypes = [ct for ct in tissue_celltypes if ct in first_pass_result.columns]
        if not tissue_celltypes:
            continue

        # subset the expression matrix for this tissue
        expr_block = first_pass_result.loc[:, tissue_celltypes]

        # iterate through hormones and check coverage
        for hormone, row in expr_block.iterrows():
            expressed_count = (row.values > 0).sum()
            if expressed_count / n_celltypes >= coverage_threshold:
                if hormone not in hormone_to_tissues:
                    hormone_to_tissues[hormone] = []
                if tissue not in hormone_to_tissues[hormone]:
                    hormone_to_tissues[hormone].append(tissue)

    return hormone_to_tissues

# def fTau(x):
#     """Calculate tissue specificity using Tau metric"""
#     x = np.array(x)
#     if np.all(np.isnan(x)) or np.min(x) < 0:
#         return None
#     if np.max(x) == 0:
#         return 0
        
#     x_norm = 1 - (x / np.max(x))
#     return np.nansum(x_norm) / (len(x) - 1)

# def fTau(x):
#     """
#     Calculate the gene specificity based on the given vector of expression values.

#     Parameters:
#         x (list or numpy array): Vector with expression values of one gene in different tissues.

#     Returns:
#         float or None: Gene specificity value or None if calculation is not possible.
#     """
#     # Ensure input is a numpy array for easier handling
#     x = np.array(x)

#     # Check if all values are not NaN
#     if not np.all(np.isnan(x)):
#         # Check if all values are non-negative
#         if np.min(x) >= 0:
#             # Check if the maximum value is not 0
#             if np.max(x) != 0:
#                 x = 1 - (x / np.max(x))
#                 res = np.sum(x)
#                 res /= (len(x) - 1)
#             else:
#                 res = 0
#         else:
#             res = None  # Expression values have to be positive
#     else:
#         res = None  # No data available for this gene

#     return res

# import numpy as np

def fTau(x):
    """
    Calculate the gene specificity (τ index).

    Rules:
    - If all values are NaN → None
    - If any negative values → None
    - If all zeros → 0
    - If exactly one value > 0 → 1
    - Otherwise compute τ index normally
    """
    # Ensure float array
    x = np.array(x, dtype=float)

    # All NaN → None
    if np.all(np.isnan(x)):
        return None

    # Any negative values → None
    if np.min(x) < 0:
        return None

    # All zeros → 0
    if np.max(x) == 0:
        return 0

    # Special case: only one positive value
    if np.sum(x > 0) == 1:
        return 1.0

    # τ index calculation
    x = 1 - (x / np.max(x))
    res = np.sum(x) / (len(x) - 1)
    return res

def collect_hormone_genes(
    hormone_producing: pd.DataFrame,
    hormones_use: Optional[Collection[str]],
    include_cols: Sequence[str],
    exclude_cols: Sequence[str],
    *,
    deduplicate: bool = True
) -> Tuple[list[str], list[str], list[str]]:
    """
    From `hormone_producing`, collect included/excluded gene lists for the
    subset of rows whose 'hormone_short' is in `hormones_use`.

    Returns (included_genes, excluded_genes, all_genes).
    """

    if not hormones_use:
        return [], [], []

    if "hormone_short" not in hormone_producing.columns:
        raise ValueError("`hormone_producing` must contain 'hormone_short'.")

    hp = hormone_producing.loc[
        hormone_producing["hormone_short"].isin(hormones_use)
    ].copy()

    def _pull_genes(df: pd.DataFrame, cols: Sequence[str]) -> list[str]:
        cols = [c for c in cols if c in df.columns]
        if not cols:
            return []
        arr = df[cols].to_numpy().ravel()
        # drop NaN
        arr = arr[pd.notna(arr)]
        # to list[str]
        genes = [str(x) for x in arr]
        if deduplicate:
            genes = list(np.unique(genes))
        return genes

    included = _pull_genes(hp, include_cols)
    excluded = _pull_genes(hp, exclude_cols)
    all_genes = included + excluded

    return all_genes, excluded


import pandas as pd

def reshape_hormone_wide(
    hormone_wide: pd.DataFrame,
    celltype_column: str = "Celltype_unique",
    tissue_col: str = "Tissue",
    split='____'
) -> pd.DataFrame:
    """
    Convert wide-format hormone expression (hormones as rows, celltypes as columns)
    into long-format DataFrame with columns: Hormone, celltype, tissue, expression.
    Returns
    -------
    pd.DataFrame
        Long-format DataFrame with columns: Hormone, <celltype_column>, AveExpression, <tissue_col>.
    """
    hormone_long = hormone_wide.reset_index()
    hormone_long.columns.values[0] = "Hormone"

    hormone_long = pd.melt(
        hormone_long,
        id_vars=["Hormone"],
        var_name=celltype_column,
        value_name="Strength"
    )

    hormone_long[tissue_col] = [
        str(i).split(split)[0] for i in hormone_long[celltype_column]
    ]
    #print(hormone_long.head(2))

    return hormone_long



def annotate_hormone_long(
    hormone_long: pd.DataFrame,
    hormone_producing: pd.DataFrame,
    celltype_column: str='Celltype_unique',
    assay: str = None
) -> pd.DataFrame:
    """
    Annotate hormone-long dataframe with display info and remove duplicates.

    Parameters
    ----------
    hormone_long : pd.DataFrame
        Long-format dataframe with at least columns ['Hormone','Celltype_unique','Strength'].
    hormone_producing : pd.DataFrame
        Reference dataframe with columns ['hormone_display','hormone_short','hormone_full','Tier'].
    assay : str, optional
        If provided, add a column 'assay' with this value.

    Returns
    -------
    pd.DataFrame
        Annotated long dataframe without duplicates.
    """
    # merge annotations
    hormone_ct_long = pd.merge(
        hormone_long,
        hormone_producing.loc[:, ['hormone_display','hormone_short','hormone_full','Tier']],
        left_on='Hormone',
        right_on='hormone_short',
        how='left'
    )

    print(f"Before removing duplicates: {hormone_ct_long.shape}")

    # drop duplicates: Hormone + Celltype_unique + AveExpression
    hormone_ct_long['tmp'] = (
        hormone_ct_long['Hormone'].astype(str)
        + hormone_ct_long[celltype_column].astype(str)
        + hormone_ct_long['Strength'].astype(str)
    )
    hormone_ct_long = hormone_ct_long.drop_duplicates('tmp')
    hormone_ct_long = hormone_ct_long.drop(columns='tmp')

    print(f"After removing duplicates: {hormone_ct_long.shape}")

    # add assay info if provided
    if assay is not None:
        hormone_ct_long['assay'] = assay

    return hormone_ct_long


# def replace_adjusted_rows(
#     hormone_long1: pd.DataFrame,
#     hormone_long_adj_coverage: pd.DataFrame,
#     adjusted_hormones_coverage: dict
# ) -> pd.DataFrame:
#     """
#     Completely replace rows in hormone_long1 with rows from hormone_long_adj_coverage
#     for (Hormone, Tissue) pairs specified in adjusted_hormones_coverage.
#     If hormone_long_adj_coverage is empty, return hormone_long1 unchanged.
#     """
#     # If replacement table is empty → return original unchanged
#     if hormone_long_adj_coverage.empty:
#         return hormone_long1.copy()

#     # Build (Hormone, Tissue) pairs that need replacement
#     pairs = pd.DataFrame(
#         [(h, t) for h, ts in adjusted_hormones_coverage.items() for t in ts],
#         columns=["Hormone", "Tissue"]
#     )

#     # Drop those pairs from original
#     mask = hormone_long1.set_index(["Hormone", "Tissue"]).index.isin(
#         pairs.set_index(["Hormone","Tissue"]).index
#     )
#     hormone_long_new = hormone_long1.loc[~mask].copy()

#     # Pick replacement rows from hormone_long_adj_coverage
#     to_add = hormone_long_adj_coverage.merge(pairs, on=["Hormone", "Tissue"], how="inner")

#     # Concat back
#     hormone_long_new = pd.concat([hormone_long_new, to_add], ignore_index=True)

#     return hormone_long_new


import pandas as pd

def replace_adjusted_rows(
    hormone_long1: pd.DataFrame,
    hormone_long_adj_coverage: pd.DataFrame,
    adjusted_hormones_coverage: dict,
    coverage_fraction: float = 1/3
) -> pd.DataFrame:
    """
    Completely replace rows in hormone_long1 with rows from hormone_long_adj_coverage
    for (Hormone, Tissue) pairs specified in adjusted_hormones_coverage.
    
    If a hormone's adjusted tissues cover more than `coverage_fraction` of all tissues,
    then all tissues for that hormone will be replaced.

    If hormone_long_adj_coverage is empty, return hormone_long1 unchanged.
    """
    # If replacement table is empty → return original unchanged
    if hormone_long_adj_coverage is None:
        return hormone_long1.copy()

    # All unique tissues in the dataset
    all_tissues = hormone_long1["Tissue"].unique().tolist()
    n_total_tissues = len(all_tissues)

    # Expand adjusted_hormones_coverage according to coverage_fraction rule
    expanded_adjusted = {}
    print(f'    Hormones that have been adjusted includes: {adjusted_hormones_coverage.keys()}')
    for h, ts in adjusted_hormones_coverage.items():
        if len(ts) / n_total_tissues >= coverage_fraction:
            # replace all tissues for this hormone
            expanded_adjusted[h] = all_tissues
            #print(h)
        else:
            expanded_adjusted[h] = ts
    
    # Build (Hormone, Tissue) pairs that need replacement
    pairs = pd.DataFrame(
        [(h, t) for h, ts in expanded_adjusted.items() for t in ts],
        columns=["Hormone", "Tissue"]
    )

    # Drop those pairs from original
    mask = hormone_long1.set_index(["Hormone", "Tissue"]).index.isin(
        pairs.set_index(["Hormone","Tissue"]).index
    )
    hormone_long_new = hormone_long1.loc[~mask].copy()

    # Pick replacement rows from hormone_long_adj_coverage
    to_add = hormone_long_adj_coverage.merge(pairs, on=["Hormone", "Tissue"], how="inner")

    # Concat back
    hormone_long_new = pd.concat([hormone_long_new, to_add], ignore_index=True)

    return hormone_long_new


def calculate_hormone_strength(
    ave_all: pd.DataFrame,
    hormone_producing: pd.DataFrame,
    include_cols: Sequence[str],
    exclude_cols: Sequence[str],
    all_included_genes: Optional[Iterable[str]] = None,
    all_excluded_genes: Optional[Iterable[str]] = None,
    thresh_expr_low: float = 0.0,           # formerly THRESH_EXPR_LOW
    thresh_pct: float = 0.0,                # formerly THRESH_PCT
    celltype_column: str = "Celltype_unique",
    thresh_included: float = 0.1,           # threshold for included genes
    thresh_excluded: float = 0.5,            # threshold for excluded genes
    assay:Optional[str] = None,
    use_precomputed:Optional[str] = None
) -> pd.DataFrame:
    """
    Compute hormone strength by:
      (1) filtering the input expression table,
      (2) building gene–celltype matrices for included and excluded genes
          using different expression thresholds,
      (3) applying hormone rules to produce a wide output.

    Parameters
    ----------
    ave_all : pd.DataFrame
        Raw/aggregate expression table used as input to `initial_filtering`.
        Expected to contain columns needed by that function
        (e.g., 'TotalCellNumber', 'ExpressedCellNumber', 'logExpression', 'Percentage', 'Gene', etc.).
    hormone_producing : pd.DataFrame
        Table of hormone definitions/rules consumed by `adopt_hormone_rules`.
    include_cols : Sequence[str]
        Column names in `hormone_producing` that define the “include” rules.
    exclude_cols : Sequence[str]
        Column names in `hormone_producing` that define the “exclude” rules.
    all_included_genes : Optional[Iterable[str]], default None
        Optional gene list passed to `initial_filtering` (used to restrict/adjust included genes).
    all_excluded_genes : Optional[Iterable[str]], default None
        Optional gene list used to subset the excluded gene–celltype matrix after it’s built.
    thresh_expr_low : float, default 0.0
        Expression threshold for `initial_filtering` (formerly THRESH_EXPR_LOW).
    thresh_pct : float, default 0.0
        Percentage threshold for `initial_filtering` (formerly THRESH_PCT).
    celltype_column : str, default "Celltype_unique"
        Column name used by `create_expression_matrix` to form gene–celltype matrices.
    thresh_included : float, default 0.1
        Expression threshold for included genes when computing `max_expression`.
    thresh_excluded : float, default 0.5
        Expression threshold for excluded genes when computing `max_expression_exclude`.
    use_precomputed : cell or nucleus or None, default None.
        Load precomputed maximum log-expression values for hormone-related genes and indicate which assay to load the precomputed thresholds for.

    Returns
    -------
    pd.DataFrame
        A wide DataFrame (hormone × cell type) with computed hormone strength.

    Notes
    -----
    This function assumes the following helpers are available in scope:
      - initial_filtering(df, THRESH_EXPR_LOW, THRESH_PCT, hormone_genes)
      - calculate_gene_manner_threshold(df, threshold_expr_gene)
      - create_expression_matrix(df, max_expression_df, celltype_column)
      - adopt_hormone_rules(hormone_producing, expr_inc, expr_exc, include_cols, exclude_cols)
    """

    # 1) Filter the input expression table
    print('    1) Filter the input expression table')
    ave_positive = initial_filtering(
        ave_all=ave_all,
        thresh_expr_low=thresh_expr_low,
        thresh_pct=thresh_pct,
        hormone_genes=all_included_genes
    )
    
    # Early exit if nothing passes the filter
    if ave_positive is None or len(ave_positive) == 0:
        return pd.DataFrame()

    print('    2) Build gene–celltype matrix for included genes (lower threshold)')
    # 2) Build gene–celltype matrix for included genes (lower threshold)
    max_expression = calculate_gene_manner_threshold(
        ave_positive, threshold_expr_gene=thresh_included,use_precomputed=use_precomputed
    )
    #print(max_expression.head(2))
    expression_matrix = create_expression_matrix(
        ave_positive, max_expression, celltype_column=celltype_column
    )

    print('    3) Build gene–celltype matrix for excluded genes (higher threshold)')
    # 3) Build gene–celltype matrix for excluded genes (higher threshold)
    max_expression_exclude = calculate_gene_manner_threshold(
        ave_positive, threshold_expr_gene=thresh_excluded,use_precomputed=use_precomputed
    )
    expression_matrix_exclude = create_expression_matrix(
        ave_positive, max_expression_exclude, celltype_column=celltype_column
    )

    # If an explicit excluded-gene list is provided, subset safely to avoid KeyError
    if all_excluded_genes is not None:
        present = [g for g in all_excluded_genes if g in expression_matrix_exclude.index]
        expression_matrix_exclude = expression_matrix_exclude.loc[present]
    else:
        # Use an empty frame with matching columns if no excluded set is given
        expression_matrix_exclude = expression_matrix_exclude.iloc[0:0]

    # 4) Apply hormone rules to compute the final wide table
    print("    4) Apply hormone rules to compute the final table")
    # if combo:
    #     hormone_wide = adopt_hormone_rules_combo(
    #     hormone_producing,
    #     expression_matrix,
    #     expression_matrix_exclude,
    #     include_cols,
    #     exclude_cols
    #     )
    # else:
    #     hormone_wide = adopt_hormone_rules(
    #     hormone_producing,
    #     expression_matrix,
    #     expression_matrix_exclude,
    #     include_cols,
    #     exclude_cols) 
    hormone_wide = adopt_hormone_rules(hormone_producing,
                                       expression_matrix,
                                       expression_matrix_exclude,
                                       include_cols,
                                       exclude_cols)
    return hormone_wide
    # print("5) Remove duplicates and also annotate the hormones")
    # hormone_long1=reshape_hormone_wide(hormone_wide)
    # mask=hormone_long1['AveExpression']>0
    # hormone_long1=hormone_long1.loc[mask]
    # final_dt=annotate_hormone_long(hormone_long1,hormone_producing,assay=assay)
    # print(f'Finally hormone2cell finds: f{str(final_dt.shape[0])} hormone-celltype pairs.')
    # return final_dt





# adjusted the hormones in all tissues
# def calculate_hormone_strength_specificity(
#     ave_all: pd.DataFrame,
#     hormone_producing: pd.DataFrame,
#     include_cols: Sequence[str],
#     exclude_cols: Sequence[str],
#     all_included_genes: list[str] | None,
#     all_excluded_genes: list[str] | None,
#     celltype_column: str = "Celltype_unique",
#     # specificity / coverage thresholds
#     specificity_threshold: float | None = 0.8,
#     coverage_threshold: float | None = 0.5,
#     tissue_col: str = "Tissue",
#     # filtering thresholds (first pass)
#     thresh_expr_low: float = 0.01,
#     thresh_pct: float = 5,
#     thresh_included_initial: float = 0.1,
#     thresh_excluded_initial: float = 0.15,
#     # filtering thresholds (adjusted pass)
#     thresh_included_adjusted: float | None = 0.5,
#     thresh_excluded_adjusted: float | None = 0.15,
# ) -> pd.DataFrame:
#     """
#     Two-pass hormone strength pipeline:
#       1) Compute first-pass hormone × celltype matrix.
#       2) Identify hormones to adjust by specificity/coverage.
#       3) Recompute adjusted hormones with stricter thresholds.
#       4) Replace original rows with adjusted rows and return the merged table.
#     """

#     # --------------------------
#     # Pass 1: initial strengths
#     # --------------------------
#     hormone_wide = calculate_hormone_strength(
#         ave_all=ave_all,
#         hormone_producing=hormone_producing,
#         include_cols=include_cols,
#         exclude_cols=exclude_cols,
#         all_included_genes=all_included_genes,
#         all_excluded_genes=all_excluded_genes,
#         thresh_expr_low=thresh_expr_low,
#         thresh_pct=thresh_pct,
#         celltype_column=celltype_column,
#         thresh_included=thresh_included_initial,
#         thresh_excluded=thresh_excluded_initial,
#     )

#     if hormone_wide is None or hormone_wide.empty:
#         return pd.DataFrame()

#     # -----------------------------------
#     # Build cell annotations for coverage
#     # -----------------------------------
#     needed_cols = {celltype_column, tissue_col}
#     if not needed_cols.issubset(ave_all.columns):
#         raise ValueError(
#             f"`ave_all` must contain columns {needed_cols}, found {set(ave_all.columns)}"
#         )
#     cellanno = (
#         ave_all.loc[:, [celltype_column, "Tissue"]]
#         .drop_duplicates(subset=[celltype_column])
#         .reset_index(drop=True)
#     )

#     # -----------------------------------
#     # Identify hormones to adjust
#     # -----------------------------------
#     adjusted_hormones = get_specificity_adjustment_hormones(
#         first_pass_result=hormone_wide,
#         specificity_threshold=specificity_threshold,
#         coverage_threshold=coverage_threshold,
#         cell_annos=cellanno,
#         tissue_col=tissue_col,
#         celltype_col=celltype_column,
#     )

#     # If nothing to adjust, return first-pass result as-is
#     if not adjusted_hormones:
#         print('No hormone need to be adjusted')
#         return hormone_wide

#     # ---------------------------------------------------
#     # Collect genes corresponding to the adjusted hormones
#     # ---------------------------------------------------
#     if "hormone_short" not in hormone_producing.columns:
#         raise ValueError("`hormone_producing` must contain 'hormone_short'.")

#     # Filter rows for adjusted hormones
#     # mask = hormone_producing["hormone_short"].isin(adjusted_hormones)
#     # hormone_producing_adjust = hormone_producing.loc[mask].copy()
#     # all_included_genes_adjust = pd.unique(hormone_producing_adjust[include_cols].values.ravel())
#     # all_included_genes_adjust = [g for g in all_included_genes_adjust if pd.notna(g)]
#     # all_excluded_genes_adjust = pd.unique(hormone_producing_adjust[exclude_cols].values.ravel())
#     # all_excluded_genes_adjust = [g for g in all_excluded_genes_adjust if pd.notna(g)]
#     # all_hormone_genes_adjust=all_included_genes_adjust+all_excluded_genes_adjust
#     all_included_genes_adjust, all_excluded_genes_adjust = collect_hormone_genes(
#         hormone_producing=hormone_producing,
#         hormones_use=adjusted_hormones,
#         include_cols=include_cols,
#         exclude_cols=exclude_cols
#     )

    
#     # -----------------------------------------------
#     # Pass 2: recompute only the adjusted hormones
#     # -----------------------------------------------
#     hormone_wide_adj = calculate_hormone_strength(
#         ave_all=ave_all,
#         hormone_producing=hormone_producing,
#         include_cols=include_cols,
#         exclude_cols=exclude_cols,
#         all_included_genes=all_included_genes_adjust,
#         all_excluded_genes=all_excluded_genes_adjust,
#         thresh_expr_low=thresh_expr_low,
#         thresh_pct=thresh_pct,
#         celltype_column=celltype_column,
#         thresh_included=thresh_included_adjusted,
#         thresh_excluded=thresh_excluded_adjusted,
#     )
#     hormone_wide_adj = hormone_wide_adj.loc[list(adjusted_hormones),:]

#     if hormone_wide_adj is None or hormone_wide_adj.empty:
#         # Nothing got recomputed; keep the first pass
#         return hormone_wide

#     # ------------------------------------------------------
#     # Merge: replace original rows with adjusted hormone rows
#     # ------------------------------------------------------
#     all_cols = hormone_wide.columns.union(hormone_wide_adj.columns)
#     print(f'The number of adjusted hormone_wide: f{str(hormone_wide_adj.shape)}')
#     print(f'The number of all hormone_wide: f{str(hormone_wide.shape)}')
#     out = pd.concat(
#         [hormone_wide.reindex(columns=all_cols).loc[~hormone_wide.index.isin(hormone_wide_adj.index)],
#          hormone_wide_adj.reindex(columns=all_cols)],
#         axis=0
#     )
#     #out=hormone_wide_adj
#     print(f'The number of all hormone_wide after merge: f{str(out.shape)}')
#     # Preserve ordering (optional): sort by index or keep original order
#     # out = out.sort_index()

#     return out


## adjusted the hormones in specific tissue when the coverage fraction below 1/3.
def calculate_hormone_strength_specificity(
    ave_all: pd.DataFrame,
    hormone_producing: pd.DataFrame,
    include_cols: Sequence[str],
    exclude_cols: Sequence[str],
    all_included_genes: list[str] | None,
    all_excluded_genes: list[str] | None,
    celltype_column: str = "Celltype_unique",
    # specificity / coverage thresholds
    specificity_threshold: float | None = 0.8,
    coverage_threshold: float | None = 0.5,
    tissue_col: str = "Tissue",
    # filtering thresholds (first pass)
    thresh_expr_low: float = 0.01,
    thresh_pct: float = 5,
    thresh_included_initial: float = 0.1,
    thresh_excluded_initial: float = 0.15,
    # filtering thresholds (adjusted pass)
    thresh_included_adjusted: float | None = 0.5,
    thresh_excluded_adjusted: float | None = 0.15,
    assay:Optional[str] = None,
    use_precomputed:Optional[str] = None
) -> pd.DataFrame:
    """
    Two-pass hormone strength pipeline:
      1) Compute first-pass hormone × celltype matrix.
      2) Identify hormones to adjust by specificity/coverage.
      3) Recompute adjusted hormones with stricter thresholds.
      4) Replace original rows with adjusted rows and return the merged table.
    """

    # --------------------------
    # Pass 1: initial strengths
    # --------------------------
    print('I. Compute initial strengths')
    hormone_wide = calculate_hormone_strength(
        ave_all=ave_all,
        hormone_producing=hormone_producing,
        include_cols=include_cols,
        exclude_cols=exclude_cols,
        all_included_genes=all_included_genes,
        all_excluded_genes=all_excluded_genes,
        thresh_expr_low=thresh_expr_low,
        thresh_pct=thresh_pct,
        celltype_column=celltype_column,
        thresh_included=thresh_included_initial,
        thresh_excluded=thresh_excluded_initial,
        use_precomputed=use_precomputed
    )

    if hormone_wide is None or hormone_wide.empty:
        return pd.DataFrame()

    # -----------------------------------
    # Build cell annotations for coverage
    # -----------------------------------
    needed_cols = {celltype_column, tissue_col}
    if not needed_cols.issubset(ave_all.columns):
        raise ValueError(
            f"`ave_all` must contain columns {needed_cols}, found {set(ave_all.columns)}"
        )
    cellanno = (
        ave_all.loc[:, [celltype_column, "Tissue"]]
        .drop_duplicates(subset=[celltype_column])
        .reset_index(drop=True)
    )

    # -----------------------------------
    # Identify hormones to adjust and recompute only the adjusted hormones
    # -----------------------------------
    print('II. Identify hormones to adjust and recompute only the adjusted hormones')
    hormone_long_adj_coverage=None
    hormone_long_adj_specificity=None
    ## identify the hormones
    if specificity_threshold is not None:
        print('   1). Identify hormones to adjust based on the specificity')
        adjusted_hormones_specificity= get_hormones_by_specificity(
            first_pass_result = hormone_wide,
            cell_annos = cellanno,
            specificity_threshold = specificity_threshold,
            tissue_col = tissue_col,
            hormone_col = "Hormone"
        )
        
        if not adjusted_hormones_specificity:
            print('   No hormone need to be adjusted due to specificity calculation')
            all_included_genes_specificity_adjust=[]
            all_excluded_genes_specificity_adjust=[]
        else:
            ## list the genes
            all_included_genes_specificity_adjust, all_excluded_genes_specificity_adjust = collect_hormone_genes(
                hormone_producing=hormone_producing,
                hormones_use=adjusted_hormones_specificity.keys(),
                include_cols=include_cols,
                exclude_cols=exclude_cols
        )
            ## calculate the hormone strength of adjusted hormones 
            adjusted_hormones=adjusted_hormones_specificity.keys()
            mask=hormone_producing['hormone_short'].isin(adjusted_hormones)
            hormone_producing1=hormone_producing.loc[mask]
            hormone_wide_adj_specificity = calculate_hormone_strength(
                ave_all=ave_all,
                hormone_producing=hormone_producing1,
                include_cols=include_cols,
                exclude_cols=exclude_cols,
                all_included_genes=all_included_genes_specificity_adjust,
                all_excluded_genes=all_excluded_genes_specificity_adjust,
                thresh_expr_low=thresh_expr_low,
                thresh_pct=thresh_pct,
                celltype_column=celltype_column,
                thresh_included=thresh_included_adjusted,
                thresh_excluded=thresh_excluded_adjusted,
                use_precomputed=use_precomputed
            )
            ## wide to long
            hormone_long_adj_coverage=reshape_hormone_wide(hormone_wide_adj_specificity,
                                                           celltype_column = celltype_column,
                                                           tissue_col= tissue_col,
                                                           split='____')
    
    ## identify the hormones
    if coverage_threshold is not None:
        print('   2). Identify hormones to adjust based on the tissue coverage')
        adjusted_hormones_coverage= get_hormones_by_coverage(
        first_pass_result = hormone_wide,
        cell_annos = cellanno,
        coverage_threshold = coverage_threshold,
        tissue_col = tissue_col,
        celltype_col = celltype_column) 
    
        if not adjusted_hormones_coverage:
            print('   No hormone need to be adjusted due to tissue coverage')
        else:
            all_included_genes_coverage_adjust, all_excluded_genes_coverage_adjust = collect_hormone_genes(
                hormone_producing=hormone_producing,
                hormones_use=adjusted_hormones_coverage.keys(),
                include_cols=include_cols,
                exclude_cols=exclude_cols
        )
            ## calculate the hormone strength of adjusted hormones 
            adjusted_hormones=adjusted_hormones_coverage.keys()
            mask=hormone_producing['hormone_short'].isin(adjusted_hormones)
            hormone_producing1=hormone_producing.loc[mask]
            hormone_wide_adj_coverage = calculate_hormone_strength(
                ave_all=ave_all,
                hormone_producing=hormone_producing1,
                include_cols=include_cols,
                exclude_cols=exclude_cols,
                all_included_genes=all_included_genes_coverage_adjust,
                all_excluded_genes=all_excluded_genes_coverage_adjust,
                thresh_expr_low=thresh_expr_low,
                thresh_pct=thresh_pct,
                celltype_column=celltype_column,
                thresh_included=thresh_included_adjusted,
                thresh_excluded=thresh_excluded_adjusted,
                use_precomputed=use_precomputed
            )
            hormone_long_adj_coverage=reshape_hormone_wide(hormone_wide_adj_coverage,
                                                          celltype_column =celltype_column,
                                                           tissue_col= tissue_col,
                                                           split='____')
    # ------------------------------------------------------
    # Merge: replace original rows with adjusted hormone rows
    # ------------------------------------------------------
    print('III. Merge: replace original rows with adjusted hormone rows')
    hormone_long1=reshape_hormone_wide(hormone_wide,
                                       celltype_column = celltype_column,
                                       tissue_col= tissue_col,
                                       split='____')
    print(hormone_long1.shape)
    if (hormone_long_adj_coverage is None ) and (hormone_long_adj_specificity is None):
        # Nothing got recomputed; keep the first pass
        final_dt=annotate_hormone_long(hormone_long1,hormone_producing,celltype_column=celltype_column,assay=assay)
    else:
        ## replace the hormone-tissue pair data
        hormone_long_new_replace=replace_adjusted_rows(
            hormone_long1=hormone_long1,
            hormone_long_adj_coverage=hormone_long_adj_coverage,
            adjusted_hormones_coverage=adjusted_hormones_coverage,
            coverage_fraction=1/3)
        
        hormone_long_new_replace =replace_adjusted_rows(
            hormone_long1=hormone_long_new_replace,
            hormone_long_adj_coverage=hormone_long_adj_specificity,
            adjusted_hormones_coverage=adjusted_hormones_specificity,
            coverage_fraction=1
        )
        print(hormone_long_new_replace.shape)
    
        mask=hormone_long_new_replace['Strength']>0
        hormone_long_new_replace=hormone_long_new_replace.loc[mask]
        final_dt=annotate_hormone_long(hormone_long_new_replace,hormone_producing,celltype_column=celltype_column,assay=assay)
    
    print(f'Finally hormone2cell finds: {str(final_dt.shape[0])} hormone-celltype pairs.')
    return final_dt


# def process_final_hormone_dt(hormone_wide: pd.DataFrame,
#                              hormone_producing: pd.DataFrame,
#                              assay:Optional[str] = None,
#                             celltype_column: str ='Celltype_unique'
#                             ) -> pd.DataFrame:
#     """
#     Convert a wide hormone × celltype matrix into a long table and append
#     hormone annotations.

#     Parameters
#     ----------
#     hormone_wide : pd.DataFrame
#         Wide matrix where index are hormones and columns are cell types;
#         values are average expression.
#     hormone_producing : pd.DataFrame
#         Annotation table containing at least the key 'hormone_short' and
#         optionally 'hormone_display', 'hormone_full', 'Tier'.
#     assay : Optional[str], default None
#         If provided, a constant 'assay' column will be added to the output.

#     Returns
#     -------
#     pd.DataFrame
#         Long-form DataFrame with columns:
#         ['Hormone', 'Celltype_unique', 'AveExpression', annotation columns..., 'assay'(optional)].
#     """
#     ## wide to long
#     hormone_long = hormone_wide.reset_index()
#     hormone_long.columns.values[0]='Hormone'
#     hormone_long=pd.melt(hormone_long, id_vars=['Hormone'], var_name=celltype_column, value_name='AveExpression')
#     hormone_long = hormone_long.loc[hormone_long['AveExpression']>0]
    
#     ## annotate Hormone
#     hormone_ct_long = pd.merge(hormone_long,hormone_producing.loc[:,['hormone_display','hormone_short','hormone_full','Tier']],
#                                left_on='Hormone',right_on='hormone_short',how='left')

#     ## remove duplicates as some hormones in hormone2cell have several defination of hormone production.
#     # hormone_ct_long['tmp'] = hormone_ct_long['Hormone'].astype('str') +hormone_ct_long['Celltype_unique'].astype('str')
#     # hormone_ct_long = hormone_ct_long.sort_values('AveExpression',ascending=False)
#     # hormone_ct_long = hormone_ct_long.drop_duplicates('tmp')
#     # hormone_ct_long.drop('tmp',axis=1,inplace=True)
#     print(f'Before removing the duplicates: f{str(hormone_ct_long.shape)}')
    
#     hormone_ct_long['tmp'] = hormone_ct_long['Hormone'].astype('str') +hormone_ct_long['Celltype_unique'].astype('str')+hormone_ct_long['AveExpression'].astype('str')
#     hormone_ct_long = hormone_ct_long.drop_duplicates('tmp')
#     hormone_ct_long.drop('tmp',axis=1,inplace=True)
#     print(f'After removing the duplicates: f{str(hormone_ct_long.shape)}')


#     ## add assay info if needed
#     if assay is not None:
#         hormone_ct_long['assay'] =assay

#     return hormone_ct_long
    



