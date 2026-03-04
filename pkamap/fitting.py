"""
pKa fitting: Henderson-Hasselbalch curve fitting and secondary-structure
annotation using ``rna_secstruct``.
"""

import numpy as np
import pandas as pd
from lmfit import Parameters, minimize
from rna_secstruct import SecStruct


# ---------------------------------------------------------------------------
# Henderson-Hasselbalch model
# ---------------------------------------------------------------------------

def henderson_hasselbalch(
    pH: np.ndarray, pKa: float, high: float, low: float
) -> np.ndarray:
    """Henderson-Hasselbalch model for DMS-MaPseq reactivity.

    At high pH the nucleobase is deprotonated and DMS-reactive (*high*);
    at low pH it is protonated and protected (*low*).

    See: https://pubs.acs.org/doi/10.1021/jacs.3c07736
    """
    return (high * 10 ** (-pKa + pH)) / (1 + 10 ** (-pKa + pH)) + low


def _objective(params, pH, data):
    """Residuals for least-squares minimisation."""
    predicted = henderson_hasselbalch(
        pH, params["pKa"], params["high"], params["low"]
    )
    return (data - predicted).flatten()


_NAN_RESULT = {
    "pKa": np.nan, "pKa_error": np.nan,
    "low": np.nan, "low_error": np.nan,
    "high": np.nan, "high_error": np.nan,
    "delta": np.nan,
}


def compute_pka(ph: np.ndarray, data: np.ndarray) -> dict:
    """Fit a Henderson-Hasselbalch curve and return pKa ± error.

    Parameters
    ----------
    ph : array-like
        pH values.
    data : array-like
        DMS reactivity values.

    Returns
    -------
    dict
        Keys: ``pKa``, ``pKa_error``, ``low``, ``low_error``, ``high``,
        ``high_error``, ``delta``.
    """
    mask = ~(np.isnan(ph) | np.isnan(data))
    ph, data = ph[mask], data[mask]

    if len(ph) < 4:
        return dict(_NAN_RESULT)

    try:
        fit_params = Parameters()
        fit_params.add("low", value=np.min(data))
        fit_params.add("high", value=np.max(data))
        fit_params.add("pKa", value=6.5, min=3, max=11)

        result = minimize(
            _objective, fit_params, args=(ph, data),
            nan_policy="propagate", max_nfev=5000,
        )
        return {
            "pKa": result.params["pKa"].value,
            "pKa_error": result.params["pKa"].stderr or np.nan,
            "low": result.params["low"].value,
            "low_error": result.params["low"].stderr or np.nan,
            "high": result.params["high"].value,
            "high_error": result.params["high"].stderr or np.nan,
            "delta": result.params["high"].value - result.params["low"].value,
        }
    except Exception as exc:
        print(f"Fitting failed: {exc}")
        return dict(_NAN_RESULT)


# ---------------------------------------------------------------------------
# Batch pKa calculation
# ---------------------------------------------------------------------------

def calculate_pka_values(df_res: pd.DataFrame) -> pd.DataFrame:
    """Fit pKa for every (position, construct, condition) combination.

    Parameters
    ----------
    df_res : pd.DataFrame
        Per-residue DataFrame from :func:`pkamap.processing.generate_residue_dataframe`.

    Returns
    -------
    pd.DataFrame
        One row per group with fitted pKa, error, and plateau values.
    """
    def _fit_group(group):
        return pd.Series(compute_pka(group["ph"].values, group["data"].values))

    group_cols = ["pos", "name", "buffer_conc", "temperature", "mg_conc"]
    pka_fitted = (
        df_res
        .groupby(group_cols)
        .apply(_fit_group, include_groups=False)
        .reset_index()
    )
    return pka_fitted


# ---------------------------------------------------------------------------
# Structure annotation
# ---------------------------------------------------------------------------

def _annotate_position(row, secstruct_lookup: dict, nt_window: int = 2) -> dict:
    """Annotate a position with motif and sequence-context information."""
    seq = row["sequence"]
    i = int(row["pos"]) - 1  # 0-based
    seq_obj = secstruct_lookup.get(row["name"])

    motif_keys = [
        "motif_type", "motif_id", "motif_token", "motif_sequence",
        "motif_structure", "motif_strands", "motif_positions",
        "motif_start_pos", "motif_end_pos",
    ]
    motif_info = {k: None for k in motif_keys}

    if seq_obj:
        for motif in seq_obj:
            if i in motif.positions:
                motif_info.update({
                    "motif_type": motif.m_type,
                    "motif_id": motif.m_id,
                    "motif_token": getattr(motif, "token", None),
                    "motif_sequence": motif.sequence,
                    "motif_structure": motif.structure,
                    "motif_strands": motif.strands,
                    "motif_positions": motif.positions,
                    "motif_start_pos": motif.start_pos,
                    "motif_end_pos": motif.end_pos,
                })
                break

    # Nucleotide context window
    L = len(seq)
    context = {
        f"nt_{offset:+d}": seq[i + offset] if 0 <= i + offset < L else ""
        for offset in range(-nt_window, nt_window + 1)
    }
    return {**motif_info, **context}


def add_structure_annotations(
    pka_fitted: pd.DataFrame,
    df: pd.DataFrame,
    nt_window: int = 2,
) -> pd.DataFrame:
    """Add secondary-structure motif and sequence-context annotations.

    Parameters
    ----------
    pka_fitted : pd.DataFrame
        Output of :func:`calculate_pka_values`.
    df : pd.DataFrame
        Construct-level DataFrame with ``sequence`` and ``structure`` columns.
    nt_window : int
        Number of flanking nucleotides to include in context columns.
    """
    df_unique = df.drop_duplicates(subset=["name"], keep="first")
    secstruct_lookup = {
        row["name"]: SecStruct(row["sequence"], row["structure"])
        for _, row in df_unique.iterrows()
    }

    merged = pka_fitted.merge(
        df_unique[["name", "sequence", "structure"]], on="name", how="left"
    )
    annotations = pd.DataFrame(
        [_annotate_position(row, secstruct_lookup, nt_window) for _, row in merged.iterrows()],
        index=merged.index,
    )
    annotated = pd.concat([merged, annotations], axis=1)

    # Convert positions to 1-based indexing
    annotated["motif_strands"] = annotated["motif_strands"].apply(
        lambda x: [[p + 1 for p in strand] for strand in x] if x else None
    )
    annotated["motif_positions"] = annotated["motif_positions"].apply(
        lambda x: [p + 1 for p in x] if x else None
    )
    for col in ("motif_start_pos", "motif_end_pos"):
        annotated[col] = annotated[col].apply(lambda x: x + 1 if x is not None else None)

    return annotated
