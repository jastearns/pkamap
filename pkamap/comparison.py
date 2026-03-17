"""
Average (weighted) pKa summary across constructs and experimental conditions.
"""

import numpy as np
import pandas as pd


# Condition helpers

Condition = tuple  # (buffer_conc, temperature, mg_conc, color, label)


def _condition_mask(df: pd.DataFrame, buffer, temp, mg) -> pd.Series:
    return (
        (df["buffer_conc"] == buffer)
        & (df["temperature"] == temp)
        & (df["mg_conc"] == mg)
    )


# Weighted pKa summary

def summarize(
    pka_filtered_df: pd.DataFrame,
    pka_unfiltered_df: pd.DataFrame,
    conditions: list[Condition],
    protonation_sites: list[dict],
) -> pd.DataFrame:
    """Compute inverse-variance-weighted pKa per motif, site, and condition.

    Parameters
    ----------
    pka_filtered_df : pd.DataFrame
        Filtered pKa values for protonation sites of interest.
    pka_unfiltered_df : pd.DataFrame
        Unfiltered pKa values (used for denominator counts).
    conditions : list of Condition
        Each entry: ``(buffer_conc, temperature, mg_conc, color, label)``.
    protonation_sites : list of dict
        Protonation site definitions (see :mod:`pkamap.nmr_reference_C01HK`).

    Returns
    -------
    pd.DataFrame
        Weighted pKa, error, and construct counts for each site × condition.
    """
    rows = []

    for site in protonation_sites:
        motif_seq = site["motif_sequence"]
        nt_m1, nt_0, nt_p1 = site["nt_-1"], site["nt_+0"], site["nt_+1"]

        good_mask = (
            (pka_filtered_df["motif_sequence"] == motif_seq)
            & (pka_filtered_df["nt_-1"] == nt_m1)
            & (pka_filtered_df["nt_+0"] == nt_0)
            & (pka_filtered_df["nt_+1"] == nt_p1)
        )
        good_df = pka_filtered_df.loc[good_mask]
        all_df = pka_unfiltered_df[pka_unfiltered_df["motif_sequence"] == motif_seq]

        for buffer, temp, mg, _color, label in conditions:
            # Denominator: total unique constructs for this motif + condition
            all_cond = all_df[_condition_mask(all_df, buffer, temp, mg)]
            n_total = all_cond["name"].nunique()

            # Numerator: constructs with valid fits
            good_cond = good_df[_condition_mask(good_df, buffer, temp, mg)]
            good_cond = good_cond[
                good_cond["pKa"].notna()
                & good_cond["pKa_error"].notna()
                & (good_cond["pKa_error"] > 0)
            ]
            per_construct = good_cond.groupby("name").first().reset_index()
            n_good = len(per_construct)

            base = {
                "motif_sequence": motif_seq,
                "nt_-1": nt_m1, "nt_+0": nt_0, "nt_+1": nt_p1,
                "condition": label,
                "buffer_conc": buffer, "temperature": temp, "mg_conc": mg,
                "n": f"{n_good}/{n_total}",
            }

            if n_good == 0:
                base.update({"weighted_pKa": np.nan, "weighted_pKa_error": np.nan})
            else:
                pkas = per_construct["pKa"].values
                errors = per_construct["pKa_error"].values
                weights = 1.0 / (errors ** 2)
                base["weighted_pKa"] = round(np.sum(weights * pkas) / np.sum(weights), 4)
                base["weighted_pKa_error"] = round(1.0 / np.sqrt(np.sum(weights)), 4)

            rows.append(base)

    return pd.DataFrame(rows)