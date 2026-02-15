
from .utils import *
import pandas as pd
import numpy as np
from typing import Optional



def hormone_strength(
    ave_all: pd.DataFrame,
    hormone_producing: pd.DataFrame,
    celltype_column: str = "Celltype_unique",
    tissue_col: Optional[str] ="Tissue",
    adjustment: bool = True,
    assay: Optional[str] = "cell",
    include_cols: str = ['hormonegene_include1', 'hormonegene_include2', 
                'hormonegene_include3', 'hormonegene_include4'],
    exclude_cols:str =['hormonegene_exclude1', 'hormonegene_exclude2'],
    # initial thresholds
    thresh_expr_low=0.01,
    thresh_pct=5,
    # specificity / coverage thresholds
    specificity_threshold=0.8,
    coverage_threshold=0.5,
    # filtering thresholds (first pass)
    thresh_included_initial=0.1,
    thresh_excluded_initial=0.15,
    # filtering thresholds (adjusted pass)
    thresh_included_adjusted=0.5,
    thresh_excluded_adjusted=0.15,
    use_precomputed:Optional[str] = None,
    max_expression_file='gene_max_expression.csv'

)-> pd.DataFrame:
    """
    End-to-end pipeline to compute hormone strength per cell type and return a
    long-form, annotated table.

    Workflow
    --------
    1) Collect all included/excluded genes from `hormone_producing`.
    2) Compute a wide hormone × cell-type matrix of strengths.
       - If `adjustment=True`, run a specificity/coverage-adjusted pass and
         replace the affected hormones.
       - If `adjustment=False`, run a single unadjusted pass.
    3) Convert the wide matrix to long format and append hormone annotations
       (and the `assay` label if provided).

    Parameters
    ----------
    ave_all : pd.DataFrame
        Input table used by downstream steps. Must contain at least the
        column named by `celltype_column`; if `adjustment=True`, it must
        also contain `tissue_col`.
    hormone_producing : pd.DataFrame
        Hormone definition/annotation table. Must contain 'hormone_short';
        may optionally include 'hormone_display', 'hormone_full', 'Tier'.
        Columns listed in `include_cols` / `exclude_cols` should hold gene IDs.
    celltype_column : str, default "Celltype_unique"
        Column name identifying cell types.
    tissue_col : Optional[str], default "Tissue"
        Column name identifying tissue; only required when `adjustment=True`.
    adjustment : bool, default True
        Whether to run the specificity/coverage-adjusted second pass.
    assay : Optional[str], default "cell"
        If provided, added as a constant 'assay' column in the output.
    include_cols, exclude_cols
        Column lists in `hormone_producing` that define include/exclude gene sets.
    thresh_expr_low, thresh_pct : float, defaults 0.1 and 0.05
        Initial filtering thresholds used inside the calculation.
    specificity_threshold, coverage_threshold : float, defaults 0.8 and 0.5
        Thresholds for identifying hormones that need adjustment (τ specificity
        and coverage across cell types within tissues).
    thresh_included_initial, thresh_excluded_initial : float
        Expression thresholds for the first pass (included/excluded genes).
    thresh_included_adjusted, thresh_excluded_adjusted : float
        Expression thresholds used for the adjusted pass.
    use_precomputed : cell or nucleus or None, default None.
        Load precomputed maximum log-expression values for hormone-related genes and indicate which assay to load the precomputed thresholds for.


    Returns
    -------
    pd.DataFrame
        Long-form DataFrame with columns like:
        ['Hormone', <celltype_column>, 'Strength', 'hormone_short',
         'hormone_display', 'hormone_full', 'Tier', 'assay'(optional)].
        'Strength' is derived from the wide matrix’s values.
    """
    # get the gene names
    hormones_use = hormone_producing['hormone_short'].dropna().unique().tolist()
    all_included_genes, all_excluded_genes = collect_hormone_genes(
        hormone_producing=hormone_producing,
        hormones_use=hormones_use,
        include_cols=include_cols,
        exclude_cols=exclude_cols
    )
    # calculate the hormone strength
    if adjustment:
        hormone_long=calculate_hormone_strength_specificity(
            ave_all = ave_all ,
            hormone_producing = hormone_producing,
            all_included_genes = all_included_genes,
            all_excluded_genes = all_excluded_genes,
            celltype_column = celltype_column,
            tissue_col=tissue_col,
            include_cols=include_cols,
            exclude_cols=exclude_cols,
            # initial thresholds
            thresh_expr_low=thresh_expr_low,
            thresh_pct=thresh_pct,
            # specificity / coverage thresholds
            specificity_threshold=specificity_threshold,
            coverage_threshold=coverage_threshold,
            # filtering thresholds (first pass)
            thresh_included_initial=thresh_included_initial,
            thresh_excluded_initial=thresh_excluded_initial,
            # filtering thresholds (adjusted pass)
            thresh_included_adjusted=thresh_included_adjusted,
            thresh_excluded_adjusted=thresh_excluded_adjusted,
            assay =assay,
            use_precomputed=use_precomputed,
            max_expression_file=max_expression_file
        )
    else:
        hormone_wide=calculate_hormone_strength(
            ave_all = ave_all ,
            hormone_producing = hormone_producing,
            include_cols=include_cols,
            exclude_cols=exclude_cols,
            all_included_genes = all_included_genes,
            all_excluded_genes = all_excluded_genes,
            celltype_column = celltype_column,
            # initial thresholds
            thresh_expr_low=thresh_expr_low,
            thresh_pct=thresh_pct,
            # filtering thresholds (first pass)
            thresh_included=thresh_included_initial,
            thresh_excluded=thresh_excluded_initial,
            assay =assay,
            use_precomputed=use_precomputed,
            max_expression_file=max_expression_file
            ) 
        #print("5) Remove duplicates and also annotate the hormones")
        hormone_long1=reshape_hormone_wide(hormone_wide,
                                           celltype_column =celltype_column,
                                           tissue_col= tissue_col,
                                           split='____')
        mask=hormone_long1['Strength']>0
        hormone_long1=hormone_long1.loc[mask]
        hormone_long=annotate_hormone_long(hormone_long1,hormone_producing,assay=assay,celltype_column=celltype_column)
        print(f'Finally hormone2cell finds: {str(hormone_long.shape[0])} hormone-celltype pairs.')
        
    # wide to long and add annotation for hormones
    # hormone_long_annotated = process_final_hormone_dt(hormone_wide=hormone_wide,
    #                               hormone_producing=hormone_producing,
    #                               assay=assay,
    #                               celltype_column = celltype_column) 
    return hormone_long


