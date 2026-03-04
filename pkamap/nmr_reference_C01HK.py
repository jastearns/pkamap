"""
NMR reference data from publications and their associated motifs.

All published pKa values and the motifs inserted into the pKaMaP library (C01HK)
are defined here. Condition information and original publications for
the NMR measurements can be found in the accompanying Google Doc:
https://docs.google.com/document/d/1WZEvpQKaSlWuZmMnYhkm5khPLuqFN3tMpkSqTuPGkl4
"""

import pandas as pd

# ---------------------------------------------------------------------------
# Published NMR pKa values
# ---------------------------------------------------------------------------

_NMR_DATA = {
    "motif_sequence": [
        "GAC&GCC", "GAU&ACC", "CGAAG&CCGAG", "GAG&CCC", "GAA&UCC",
        "UAU&ACA", "UAA&UCA", "CCC&GAUG", "CCGAGC&GGAG", "GAC&GC",
        "AGA&UAU", "GAG&CAU",
    ],
    "nt_-1": ["G", "G", "A", "G", "G", "U", "U", "G", "G", "G", "U", "C"],
    "nt_+0": ["A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A"],
    "nt_+1": ["C", "U", "G", "G", "A", "U", "A", "U", "G", "C", "U", "U"],
    "nmr_pKa": [8.1, 7.84, 6.4, 8.00, 7.28, 7.09, 6.51, 6.0, 6.4, 5.84, 6.21, 5.65],
    "nmr_pKa_error": [0.06, 0.05, 0.1, 0.06, 0.08, 0.03, 0.04, 0.1, 0.1, 0.08, 0.05, 0.05],
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
    ],
}

NMR_REFERENCE: pd.DataFrame = pd.DataFrame(_NMR_DATA)
"""Published NMR pKa values used as ground truth for benchmarking."""

# ---------------------------------------------------------------------------
# Motif lists
# ---------------------------------------------------------------------------

MOTIFS_WITH_NMR: list[str] = [
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

ALL_LIBRARY_MOTIFS: list[str] = MOTIFS_WITH_NMR + [
    "ACAAG&UCCAU",   # NMR pka data not yet public... NMR confirms protonation of 4th and 9th residues!
    "GCU&AAC",       # Protonation confirmed by NMR but no published pKa value... the 5th residue is protonated!
]
"""All motifs inserted into the pKaMaP comparison library."""

# ---------------------------------------------------------------------------
# NMR-verified protonation sites (motif + specific residue position)
# ---------------------------------------------------------------------------

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
