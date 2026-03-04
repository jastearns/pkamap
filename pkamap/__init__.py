"""
pKaMaP: pKa Calculation from DMS-MaPseq

Determine RNA nucleotide pKa values from changes in DMS reactivity in response to pH.
Dependencies:
    - seq_tools (https://github.com/jyesselm/seq_tools)
    - rna_secstruct (https://github.com/jyesselm/rna_secstruct)
"""

from pkamap.processing import process_json_files, generate_residue_dataframe
from pkamap.fitting import calculate_pka_values, add_structure_annotations
from pkamap.filtering import apply_quality_filters
from pkamap.comparison import summarize_pka, compare_to_nmr
from pkamap.nmr_reference_C01HK import NMR_REFERENCE, MOTIFS_WITH_NMR, PROTONATION_SITES

__version__ = "0.1.0"
