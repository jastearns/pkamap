"""
Comparison to NMR reference values and weighted pKa summary.
"""

import numpy as np
import pandas as pd
from scipy import stats


# ---------------------------------------------------------------------------
# Condition helpers
# ---------------------------------------------------------------------------

Condition = tuple  # (buffer_conc, temperature, mg_conc, color, label)


def _condition_mask(df: pd.DataFrame, buffer, temp, mg) -> pd.Series:
    return (
        (df["buffer_conc"] == buffer)
        & (df["temperature"] == temp)
        & (df["mg_conc"] == mg)
    )


# ---------------------------------------------------------------------------
# Weighted pKa summary
# ---------------------------------------------------------------------------

def summarize_pka(
    pkamap_df: pd.DataFrame,
    pka_all_df: pd.DataFrame,
    conditions: list[Condition],
    protonation_sites: list[dict],
) -> pd.DataFrame:
    """Compute inverse-variance-weighted pKa per motif, site, and condition.

    Parameters
    ----------
    pkamap_df : pd.DataFrame
        Filtered pKa values for protonation sites of interest.
    pka_all_df : pd.DataFrame
        Unfiltered pKa values (used for denominator counts).
    conditions : list of Condition
        Each entry: ``(buffer_conc, temperature, mg_conc, color, label)``.
    protonation_sites : list of dict
        Protonation site definitions (see :mod:`pkamap.reference`).

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
            (pkamap_df["motif_sequence"] == motif_seq)
            & (pkamap_df["nt_-1"] == nt_m1)
            & (pkamap_df["nt_+0"] == nt_0)
            & (pkamap_df["nt_+1"] == nt_p1)
        )
        good_df = pkamap_df.loc[good_mask]
        all_df = pka_all_df[pka_all_df["motif_sequence"] == motif_seq]

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


# ---------------------------------------------------------------------------
# Comparison to NMR (correlation plot + stats)
# ---------------------------------------------------------------------------

def _parse_n_fraction(n_str):
    """Parse '4/5' into (4, 5). Returns (None, None) on failure."""
    try:
        num, denom = n_str.split("/")
        return int(num), int(denom)
    except (ValueError, AttributeError):
        return None, None


def _plot_single_comparison(pkamap_data, nmr_data, label, color, ax, panel_letter):
    """Plot one pKaMaP-vs-NMR panel. Returns stats dict or None."""
    # Drop invalid rows
    pkamap_data = pkamap_data[
        pkamap_data["weighted_pKa"].notna() & pkamap_data["weighted_pKa_error"].notna()
    ]

    merged = pd.merge(
        pkamap_data, nmr_data,
        on=["nt_-1", "nt_+0", "nt_+1", "motif_sequence"],
        how="inner",
    )
    if merged.empty:
        ax.set_visible(False)
        return None

    x = merged["nmr_pKa"].values
    y = merged["weighted_pKa"].values
    x_err = merged["nmr_pKa_error"].values
    y_err = merged["weighted_pKa_error"].values
    motif_labels = merged["motif_sequence"].values
    n_labels = merged["n"].values if "n" in merged.columns else None

    # Count motifs: how many made it onto the plot / how many in NMR reference
    n_motifs_plotted = len(merged)
    n_motifs_total = len(nmr_data)

    valid = ~(np.isnan(x) | np.isnan(y) | np.isinf(x) | np.isinf(y))
    x, y, x_err, y_err = x[valid], y[valid], x_err[valid], y_err[valid]
    motif_labels = motif_labels[valid]
    if n_labels is not None:
        n_labels = n_labels[valid]

    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        ax.set_visible(False)
        return None

    slope, intercept, r, p, se = stats.linregress(x, y)
    r_sq = r ** 2

    # Data points
    ax.errorbar(
        x, y, xerr=x_err, yerr=y_err, fmt="o", color=color,
        alpha=0.6, capsize=3, capthick=1.5, markersize=6, elinewidth=1.5,
    )

    # Label each point with motif name and construct fraction
    if n_labels is not None:
        for xi, yi, motif, n_frac in zip(x, y, motif_labels, n_labels):
            ax.annotate(
                f"{motif} ({n_frac})",
                (xi, yi), textcoords="offset points", xytext=(6, 6),
                fontsize=6, alpha=0.8, color="0.3",
            )
    else:
        for xi, yi, motif in zip(x, y, motif_labels):
            ax.annotate(
                motif, (xi, yi), textcoords="offset points", xytext=(6, 6),
                fontsize=6, alpha=0.8, color="0.3",
            )

    # Best-fit line
    x_line = np.linspace(x.min(), x.max(), 100)
    ax.plot(x_line, slope * x_line + intercept, "darkorange", ls="--", lw=2, label="Best fit")

    # y = x reference
    lims = [min(x.min(), y.min()), max(x.max(), y.max())]
    buf = (lims[1] - lims[0]) * 0.05
    ax.plot(lims, lims, "k--", alpha=0.3, lw=1, label="y = x")
    ax.set_xlim(lims[0] - buf, lims[1] + buf)
    ax.set_ylim(lims[0] - buf, lims[1] + buf)

    # Annotations
    n_display = f"{n_motifs_plotted}/{n_motifs_total}"
    ax.text(
        0.05, 0.95,
        f"$R^2 = {r_sq:.4f}$\n$y = {slope:.2f}x + {intercept:.2f}$\n$n = {n_display}$",
        transform=ax.transAxes, fontsize=12, va="top",
        bbox=dict(boxstyle="round,pad=0.5", fc="wheat", alpha=0.9),
    )
    ax.annotate(
        panel_letter, xy=(0, 1), xycoords="axes fraction",
        xytext=(-40, 20), textcoords="offset points",
        fontsize=22, fontweight="bold", va="bottom", ha="left",
        annotation_clip=False,
    )
    ax.set_xlabel("NMR Reference pKa", fontsize=14)
    ax.set_ylabel(f"pKaMaP pKa ({label})", fontsize=14)
    ax.set_title(f"pKaMaP ({label}) vs. NMR", fontsize=14, fontweight="bold", pad=30)
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3)

    return {
        "condition": label, "n": n_display, "R_squared": r_sq,
        "slope": slope, "intercept": intercept, "p_value": p, "std_err": se,
    }


def compare_to_nmr(
    summary_df: pd.DataFrame,
    nmr_ref_df: pd.DataFrame,
    conditions: list[Condition] | None = None,
    motif_filter: str | list[str] | None = None,
    plot_combined: bool = True,
) -> pd.DataFrame:
    """Plot pKaMaP vs. NMR correlation for each experimental condition.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Output of :func:`summarize_pka`.
    nmr_ref_df : pd.DataFrame
        NMR reference values (see :data:`pkamap.reference.NMR_REFERENCE`).
    conditions : list of Condition, optional
        Experimental conditions to plot separately.
    motif_filter : str or list of str, optional
        Restrict to specific motif(s).
    plot_combined : bool
        If ``True``, add a panel combining all conditions.

    Returns
    -------
    pd.DataFrame
        Per-condition R², slope, intercept, and sample size.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    if motif_filter is not None:
        motifs = [motif_filter] if isinstance(motif_filter, str) else motif_filter
        summary_df = summary_df[summary_df["motif_sequence"].isin(motifs)]
        nmr_ref_df = nmr_ref_df[nmr_ref_df["motif_sequence"].isin(motifs)]

    panels = []  # (data_subset, label, color)

    if conditions is not None:
        for buffer, temp, mg, color, label in conditions:
            mask = _condition_mask(summary_df, buffer, temp, mg)
            panels.append((summary_df[mask], label, color))
        if plot_combined and len(conditions) > 1:
            mask = pd.Series(False, index=summary_df.index)
            for buffer, temp, mg, _, _ in conditions:
                mask |= _condition_mask(summary_df, buffer, temp, mg)
            panels.append((summary_df[mask], "All conditions", "royalblue"))
    else:
        panels.append((summary_df, "All data", "royalblue"))

    n_panels = len(panels)
    ncols = min(n_panels, 3)
    nrows = int(np.ceil(n_panels / ncols))

    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(nrows, ncols, figsize=(8 * ncols, 8 * nrows), squeeze=False)
    axes_flat = axes.flatten()
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    results = []
    for i, (data, label, color) in enumerate(panels):
        result = _plot_single_comparison(data, nmr_ref_df, label, color, axes_flat[i], letters[i])
        if result:
            results.append(result)

    for j in range(n_panels, len(axes_flat)):
        axes_flat[j].set_visible(False)

    plt.tight_layout()
    plt.show()

    return pd.DataFrame(results)