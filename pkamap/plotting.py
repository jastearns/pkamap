"""
Plotting: titration curves for individual residue fits, showing DMS reactivity
vs. pH with fitted to Henderson-Hasselbalch equation... error shown as transluscent border to fit
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pkamap.fitting import henderson_hasselbalch


# Internal helpers

def _get_condition_data(df_res, pka_df, name, pos, buffer, temp, mg):
    """Extract raw data and fitted parameters for one position + condition."""
    res_mask = (
        (df_res["name"] == name) & (df_res["pos"] == pos)
        & (df_res["buffer_conc"] == buffer) & (df_res["temperature"] == temp)
        & (df_res["mg_conc"] == mg)
    )
    res_subset = df_res.loc[res_mask]
    if len(res_subset) == 0:
        return None

    ph = res_subset["ph"].values
    data = res_subset["data"].values
    valid = ~(np.isnan(ph) | np.isnan(data))
    ph, data = ph[valid], data[valid]
    if len(ph) < 4:
        return None

    pka_mask = (
        (pka_df["name"] == name) & (pka_df["pos"] == pos)
        & (pka_df["buffer_conc"] == buffer) & (pka_df["temperature"] == temp)
        & (pka_df["mg_conc"] == mg) & (pka_df["pKa"].notna())
    )
    pka_subset = pka_df.loc[pka_mask]
    if len(pka_subset) == 0:
        return None

    row = pka_subset.iloc[0]
    return {
        "ph": ph, "data": data,
        "pKa": row["pKa"], "pKa_error": row["pKa_error"],
        "low": row["low"], "high": row["high"],
    }


def _collect_position_data(pka_df, df_res, motif_pos_idx, conditions):
    """Collect data for all constructs at one motif position."""
    pos_df = pka_df[(pka_df["motif_pos_idx"] == motif_pos_idx) & (pka_df["pKa"].notna())]
    result = {}
    for name in pos_df["name"].unique():
        pos = pos_df[pos_df["name"] == name].iloc[0]["pos"]
        cond_list = []
        for buffer, temp, mg, color, label in conditions:
            cond_data = _get_condition_data(df_res, pka_df, name, pos, buffer, temp, mg)
            if cond_data:
                cond_data["color"] = color
                cond_data["label"] = label
                cond_list.append(cond_data)
        if cond_list:
            result[name] = cond_list
    return result


def _count_total(pka_df, motif_pos_idx, conditions):
    """Count total constructs per condition (for denominator display)."""
    pos_df = pka_df[pka_df["motif_pos_idx"] == motif_pos_idx]
    counts = {}
    for buffer, temp, mg, _color, label in conditions:
        n = pos_df[
            (pos_df["buffer_conc"] == buffer)
            & (pos_df["temperature"] == temp)
            & (pos_df["mg_conc"] == mg)
        ]["name"].nunique()
        counts[label] = n
    return counts


def _count_filtered(pos_data, conditions):
    """Count constructs with good fits per condition."""
    counts = {}
    for *_, label in conditions:
        n = sum(
            1 for cond_list in pos_data.values()
            if any(c["label"] == label for c in cond_list)
        )
        counts[label] = n
    return counts


def _assign_motif_pos_idx(df: pd.DataFrame) -> pd.DataFrame:
    """Add a ``motif_pos_idx`` column: position's index within its motif."""
    df = df.copy()
    df["motif_pos_idx"] = df.apply(
        lambda r: (
            r["motif_positions"].index(r["pos"])
            if r["motif_positions"] and r["pos"] in r["motif_positions"]
            else None
        ),
        axis=1,
    )
    return df


def _plot_construct(ax, construct_name, cond_list):
    """Plot one construct panel with overlaid conditions."""
    for c in cond_list:
        color, pKa, pKa_err = c["color"], c["pKa"], c["pKa_error"]

        ax.scatter(
            c["ph"], c["data"], c=color, s=50, zorder=3,
            edgecolors="white", linewidths=0.5, alpha=0.65,
        )
        ph_smooth = np.linspace(c["ph"].min() - 0.5, c["ph"].max() + 0.5, 200)
        label_text = (
            f"{c['label']}: pKa={pKa:.2f}±{pKa_err:.2f}"
            if pd.notna(pKa_err) else f"{c['label']}: pKa={pKa:.2f}"
        )
        ax.plot(
            ph_smooth,
            henderson_hasselbalch(ph_smooth, pKa, c["high"], c["low"]),
            color=color, linewidth=2.5, label=label_text,
        )
        ax.axvline(x=pKa, color=color, linestyle="--", linewidth=1.5, alpha=0.8)
        if pd.notna(pKa_err):
            ax.axvspan(pKa - pKa_err, pKa + pKa_err, color=color, alpha=0.10)

    ax.set_title(construct_name[:30], fontsize=11, fontweight="bold")
    ax.set_xlabel("pH", fontsize=11)
    ax.set_ylabel("DMS Reactivity", fontsize=11)
    ax.grid(True, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=9, loc="best")


# Public API

def plot_titration_curves(
    df_res: pd.DataFrame,
    pka_df: pd.DataFrame,
    conditions: list[tuple],
    pka_unfiltered: pd.DataFrame,
    motif_filter: str | list[str] | None = None,
):
    """Plot titration curves for each motif position across constructs.

    Parameters
    ----------
    df_res : pd.DataFrame
        Per-residue data.
    pka_df : pd.DataFrame
        Annotated pKa DataFrame (filtered).
    conditions : list of tuple
        Each entry: ``(buffer_conc, temperature, mg_conc, color, label)``.
    pka_unfiltered : pd.DataFrame
        Unfiltered pKa DataFrame (used for denominator in construct fractions).
    motif_filter : str or list, optional
        Restrict to specific motif(s).
    """
    all_motifs = pka_df["motif_sequence"].dropna().unique()
    if motif_filter is None:
        motifs = all_motifs
    elif isinstance(motif_filter, str):
        motifs = [motif_filter]
    else:
        motifs = motif_filter

    for motif_seq in motifs:
        if motif_seq not in all_motifs:
            continue

        pka_sub = _assign_motif_pos_idx(pka_df[pka_df["motif_sequence"] == motif_seq])

        pka_unfilt_sub = _assign_motif_pos_idx(pka_unfiltered[pka_unfiltered["motif_sequence"] == motif_seq])

        residues = list(motif_seq.replace("&", ""))
        for idx, res in enumerate(residues):
            data = _collect_position_data(pka_sub, df_res, idx, conditions)
            if not data:
                continue

            filtered_counts = _count_filtered(data, conditions)
            total_counts = _count_total(pka_unfilt_sub, idx, conditions)
            count_parts = [
                f"{label}: {filtered_counts[label]}/{total_counts.get(label, '?')}"
                for *_, label in conditions
            ]

            n = len(data)
            ncols = min(3, n)
            nrows = int(np.ceil(n / ncols))

            fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 5 * nrows))
            axes = np.atleast_1d(axes).flatten()

            fig.suptitle(
                f"{motif_seq} — Position {idx + 1}: {res}",
                fontsize=16, fontweight="bold", y=1.03,
            )
            fig.text(0.5, 0.98, "  |  ".join(count_parts), ha="center", fontsize=10, style="italic")

            for ax, (name, cond_list) in zip(axes, data.items()):
                _plot_construct(ax, name, cond_list)
            for ax in axes[n:]:
                ax.set_visible(False)

            plt.tight_layout(rect=[0, 0, 1, 0.96])
            plt.show()