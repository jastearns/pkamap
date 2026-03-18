"""
Additional processing for 'C01HK' (pKa comparison library) 
"""

import numpy as np
import pandas as pd
from scipy import stats

from pkamap.comparison import _condition_mask, Condition

# Published NMR pKa values

_NMR_DATA = {
    "motif_sequence": [
        "GAC&GCC", "GAU&ACC", "CGAAG&CCGAG", "GAG&CCC", "GAA&UCC",
        "UAU&ACA", "UAA&UCA", "CCC&GAUG", "CCGAGC&GGAG", "GAC&GC",
        "AGA&UAU", "GAG&CAU", "CAGAAACAG&CGUAUAUUACG",
    ],
    "nt_-1": ["G", "G", "A", "G", "G", "U", "U", "G", "G", "G", "U", "C", "U"],
    "nt_+0": ["A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A"],
    "nt_+1": ["C", "U", "G", "G", "A", "U", "A", "U", "G", "C", "U", "U", "U"],
    "nmr_pKa": [8.1, 7.84, 6.4, 8.00, 7.28, 7.09, 6.51, 6.0, 6.4, 5.84, 6.21, 5.65, 5.8],
    "nmr_pKa_error": [0.06, 0.05, 0.1, 0.06, 0.08, 0.03, 0.04, 0.1, 0.1, 0.08, 0.05, 0.05, 0.2],
    "source": [
        "https://pubs.acs.org/doi/full/10.1021/bi400768q",
        "https://pubs.acs.org/doi/full/10.1021/bi400768q",
        "https://rnajournal.cshlp.org/content/6/12/1821.long",
        "https://pubs.acs.org/doi/full/10.1021/bi400768q",
        "https://pubs.acs.org/doi/full/10.1021/bi400768q",
        "https://pubs.acs.org/doi/full/10.1021/bi400768q",
        "https://pubs.acs.org/doi/full/10.1021/bi400768q",
        "https://www.nature.com/articles/nsb800",
        "https://pubs.acs.org/doi/10.1021/ja9640051",
        "https://www.nature.com/articles/s41589-020-00667-5",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC9116285/",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC9116285/",
        "https://pubs.acs.org/doi/10.1021/bi001976r",
    ],
}

NMR_REFERENCE: pd.DataFrame = pd.DataFrame(_NMR_DATA)
"""Published NMR pKa values used for comparison"""

# Motif lists

MOTIFS_PKA: list[str] = [
    "GAC&GCC",
    "GAU&ACC",
    "CGAAG&CCGAG",
    "GAG&CCC",
    "GAA&UCC",
    "UAU&ACA",
    "UAA&UCA",
    "CCC&GAUG",
    "CCGAGC&GGAG",
    "GAC&GC",
    "AGA&UAU",
    "GAG&CAU",
    "CAGAAACAG&CGUAUAUUACG",
]
"""Motifs in the library that have published NMR pKa values."""

MOTIFS_PROTONATION: list[str] = MOTIFS_PKA + [
    "ACAAG&UCCAU",   # NMR pka value not yet public... NMR confirms protonation of 4th and 9th residues!
    "GCU&AAC",       # Protonation confirmed by NMR but no published pKa value... the 5th residue is protonated!
]
"""All motifs with NMR-confirmed protonation events."""

# NMR-verified protonation sites

PROTONATION_SITES: list[dict] = [
    {"motif_sequence": "GAC&GCC",                   "nt_-1": "G", "nt_+0": "A", "nt_+1": "C"},
    {"motif_sequence": "GAU&ACC",                   "nt_-1": "G", "nt_+0": "A", "nt_+1": "U"},
    {"motif_sequence": "CGAAG&CCGAG",               "nt_-1": "A", "nt_+0": "A", "nt_+1": "G"},
    {"motif_sequence": "GAG&CCC",                   "nt_-1": "G", "nt_+0": "A", "nt_+1": "G"},
    {"motif_sequence": "GAA&UCC",                   "nt_-1": "G", "nt_+0": "A", "nt_+1": "A"},
    {"motif_sequence": "UAU&ACA",                   "nt_-1": "U", "nt_+0": "A", "nt_+1": "U"},
    {"motif_sequence": "UAA&UCA",                   "nt_-1": "U", "nt_+0": "A", "nt_+1": "A"},
    {"motif_sequence": "CCC&GAUG",                  "nt_-1": "G", "nt_+0": "A", "nt_+1": "U"},
    {"motif_sequence": "CCGAGC&GGAG",               "nt_-1": "G", "nt_+0": "A", "nt_+1": "G"},
    {"motif_sequence": "GAC&GC",                    "nt_-1": "G", "nt_+0": "A", "nt_+1": "C"},
    {"motif_sequence": "AGA&UAU",                   "nt_-1": "U", "nt_+0": "A", "nt_+1": "U"},
    {"motif_sequence": "GAG&CAU",                   "nt_-1": "C", "nt_+0": "A", "nt_+1": "U"},
    {"motif_sequence": "CAGAAACAG&CGUAUAUUACG",     "nt_-1": "U", "nt_+0": "A", "nt_+1": "U"},
]
"""Specific residue positions with NMR-confirmed protonation events."""

PROTONATION_SITES_ALL: list[dict] = PROTONATION_SITES + [
    {"motif_sequence": "ACAAG&UCCAU", "nt_-1": "A", "nt_+0": "A", "nt_+1": "G"},
    {"motif_sequence": "ACAAG&UCCAU", "nt_-1": "C", "nt_+0": "A", "nt_+1": "U"},
    {"motif_sequence": "GCU&AAC",     "nt_-1": "A", "nt_+0": "A", "nt_+1": "C"},
]
"""All protonation sites including those without published pKa values."""


def filter_to_protonation_sites(
    df: pd.DataFrame,
    sites: list[dict] | None = None,
) -> pd.DataFrame:
    """Filter a pKa DataFrame to rows matching known protonation sites.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns ``motif_sequence``, ``nt_-1``, ``nt_+0``, ``nt_+1``.
    sites : list of dict, optional
        Protonation site definitions.  Defaults to :data:`PROTONATION_SITES`.
    """
    if sites is None:
        sites = PROTONATION_SITES

    frames = []
    for site in sites:
        mask = (
            (df["motif_sequence"] == site["motif_sequence"])
            & (df["nt_-1"] == site["nt_-1"])
            & (df["nt_+0"] == site["nt_+0"])
            & (df["nt_+1"] == site["nt_+1"])
        )
        frames.append(df[mask])
    return pd.concat(frames, ignore_index=True)


# Comparison to NMR (correlation plot + stats)

def _plot_single_comparison(summary_data, nmr_data, label, color, ax, panel_letter, n_conditions=1):
    """Plot one pKaMaP-vs-NMR panel. Returns stats dict or None."""
    summary_data = summary_data[
        summary_data["weighted_pKa"].notna() & summary_data["weighted_pKa_error"].notna()
    ]

    merged = pd.merge(
        summary_data, nmr_data,
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

    n_motifs_plotted = len(merged)
    n_motifs_total = len(nmr_data) * n_conditions

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

    ax.errorbar(
        x, y, xerr=x_err, yerr=y_err, fmt="o", color=color,
        alpha=0.6, capsize=3, capthick=1.5, markersize=6, elinewidth=1.5,
    )

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

    x_line = np.linspace(x.min(), x.max(), 100)
    ax.plot(x_line, slope * x_line + intercept, "darkorange", ls="--", lw=2, label="Best fit")

    lims = [min(x.min(), y.min()), max(x.max(), y.max())]
    buf = (lims[1] - lims[0]) * 0.05
    ax.plot(lims, lims, "k--", alpha=0.3, lw=1, label="y = x")
    ax.set_xlim(lims[0] - buf, lims[1] + buf)
    ax.set_ylim(lims[0] - buf, lims[1] + buf)

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
        Output of :func:`pkamap.comparison.summarize`.
    nmr_ref_df : pd.DataFrame
        NMR reference values (e.g. :data:`NMR_REFERENCE`).
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

    panels = []

    if conditions is not None:
        for buffer, temp, mg, color, label in conditions:
            mask = _condition_mask(summary_df, buffer, temp, mg)
            panels.append((summary_df[mask], label, color, 1))
        if plot_combined and len(conditions) > 1:
            mask = pd.Series(False, index=summary_df.index)
            for buffer, temp, mg, _, _ in conditions:
                mask |= _condition_mask(summary_df, buffer, temp, mg)
            panels.append((summary_df[mask], "All conditions", "royalblue", len(conditions)))
    else:
        panels.append((summary_df, "All data", "royalblue", 1))

    n_panels = len(panels)
    ncols = min(n_panels, 3)
    nrows = int(np.ceil(n_panels / ncols))

    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(nrows, ncols, figsize=(8 * ncols, 8 * nrows), squeeze=False)
    axes_flat = axes.flatten()
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    results = []
    for i, (data, label, color, n_cond) in enumerate(panels):
        result = _plot_single_comparison(data, nmr_ref_df, label, color, axes_flat[i], letters[i], n_cond)
        if result:
            results.append(result)

    for j in range(n_panels, len(axes_flat)):
        axes_flat[j].set_visible(False)

    plt.tight_layout()
    plt.show()

    return pd.DataFrame(results)
