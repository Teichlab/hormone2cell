
from typing import Literal
import pandas as pd


def combine_assay(cell: pd.DataFrame, nucleus: pd.DataFrame,celltype_column: str) -> pd.DataFrame:
    """
    Combine cell and nucleus assay results:
    - Build Hormone_CT key = Hormone + '___' + Celltype_unique
    - Mark rows present only in cell / only in nucleus / in both
    - For pairs present in both, keep the row with the highest Strength
    - Drop optional columns 
    """
    required = {"Hormone", celltype_column, "Strength"}
    miss_cell = required - set(cell.columns)
    miss_nuc  = required - set(nucleus.columns)
    if miss_cell or miss_nuc:
        raise ValueError(f"Missing columns: cell {miss_cell}, nucleus {miss_nuc}")

    # Work on copies; build key without turning NaN into the literal string "nan"
    c = cell.copy()
    n = nucleus.copy()
    c["Hormone_CT"] = c["Hormone"].astype("string").str.cat(
        c[celltype_column].astype("string"), sep="___", na_rep=None
    )
    n["Hormone_CT"] = n["Hormone"].astype("string").str.cat(
        n[celltype_column].astype("string"), sep="___", na_rep=None
    )

    # Compute membership using Index set ops (fast, concise)
    idx_c = pd.Index(c["Hormone_CT"])
    idx_n = pd.Index(n["Hormone_CT"])
    only_c = idx_c.difference(idx_n)
    only_n = idx_n.difference(idx_c)
    both   = idx_c.intersection(idx_n)

    # Slice and tag
    dt_cell_only = c.loc[c["Hormone_CT"].isin(only_c)].assign(assay="cell_only")
    dt_nuc_only  = n.loc[n["Hormone_CT"].isin(only_n)].assign(assay="nucleus_only")

    # For keys in both, concat then keep the max AveExpression per Hormone_CT
    dt_both = pd.concat(
        [
            c.loc[c["Hormone_CT"].isin(both)].assign(assay="both"),
            n.loc[n["Hormone_CT"].isin(both)].assign(assay="both"),
        ],
        ignore_index=True,
    )
    if not dt_both.empty:
        dt_both = (dt_both
           .dropna(subset=["Strength"])
           .sort_values(["Hormone_CT", "Strength"], ascending=[True, False])
           .drop_duplicates("Hormone_CT", keep="first")
           .reset_index(drop=True))

    # Combine all parts
    out = pd.concat([dt_cell_only, dt_nuc_only, dt_both], ignore_index=True)

    # Clean up
    out = out.drop(columns=[c for c in ["Type", "tmp"] if c in out.columns], errors="ignore")
    out = out.loc[out["Strength"].notna()]#.rename(columns={"Strength": "Strength"})
    #out = out.drop('Hormone_CT',axis=1)

    return out
