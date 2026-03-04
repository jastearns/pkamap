"""
Quality filters for pKa fits.

These thresholds distinguish genuine protonation events from noise.  Most
nucleotides in an RNA sequence are *not* expected to exhibit pH-dependent
reactivity changes, so aggressive filtering is expected — the filters
enrich for residues that show a clear sigmoidal DMS-reactivity response.
"""

import pandas as pd

# Default thresholds (empirically determined)
DEFAULT_FILTERS = {
    "min_delta": 0.00024,       # high − low  (minimum signal change)
    "max_pka_error": 0.56,      # pKa fit uncertainty
    "min_high": 0.00067,        # minimum upper-plateau reactivity
    "max_high_error": 0.02471,  # upper-plateau fit uncertainty
    "pka_range": (4.0, 9.0),    # physiologically relevant pKa window
}


def apply_quality_filters(
    df: pd.DataFrame,
    min_delta: float = DEFAULT_FILTERS["min_delta"],
    max_pka_error: float = DEFAULT_FILTERS["max_pka_error"],
    min_high: float = DEFAULT_FILTERS["min_high"],
    max_high_error: float = DEFAULT_FILTERS["max_high_error"],
    pka_range: tuple[float, float] = DEFAULT_FILTERS["pka_range"],
    inplace: bool = False,
) -> pd.DataFrame:
    """Apply quality-control filters to a pKa DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns ``high``, ``low``, ``pKa``, ``pKa_error``,
        ``high_error``.
    min_delta : float
        Minimum difference between high and low plateaus.
    max_pka_error : float
        Maximum acceptable pKa fitting error.
    min_high : float
        Minimum upper-plateau (deprotonated) DMS reactivity.
    max_high_error : float
        Maximum acceptable upper-plateau fitting error.
    pka_range : tuple of float
        Allowed pKa range ``(low, high)``.
    inplace : bool
        If ``True``, filter the DataFrame in place and return it.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame.
    """
    mask = (
        ((df["high"] - df["low"]) > min_delta)
        & (df["pKa_error"] < max_pka_error)
        & (df["high"] > min_high)
        & (df["high_error"] < max_high_error)
        & df["pKa"].between(*pka_range)
    )

    before = len(df)
    if inplace:
        df.drop(df[~mask].index, inplace=True)
        print(f"Filtered: {before} → {len(df)} rows ({len(df)/before:.1%} retained)")
        return df
    else:
        out = df[mask].copy()
        print(f"Filtered: {before} → {len(out)} rows ({len(out)/before:.1%} retained)")
        return out
