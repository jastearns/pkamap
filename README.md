# pKaMaP

**RNA pKa Calculation Using DMS-MaPseq**

Calculate RNA nucleotide pKa values from pH-dependent DMS-MaPseq reactivity data (@N1 and N3 of A and C residues).

## Overview

pKaMaP fits the Henderson-Hasselbalch equation to DMS reactivity measured across a pH titration series to extract per-nucleotide pKa values. The pipeline:

1. Loads sequencing data from run JSON files (produced by [dreem](https://github.com/jyesselm/dreem) or equivalent)
2. Extracts experimental conditions (pH, buffer, temperature, Mg²⁺) from construct naming conventions
3. Trims primer sequences and converts to per-residue format
4. Fits pKa values using nonlinear least-squares minimization
5. Annotates each position with secondary-structure motif information via [rna_secstruct](https://github.com/jyesselm/rna_secstruct)
6. Filters for good fits using quality thresholds
7. Further analyze and contextualize data
8. Create plots comparing pKaMaP data to NMR results (for the C01HK library)

## Installation

```bash
pip install -e .
```

### Dependencies

- Python ≥ 3.10
- [seq_tools](https://github.com/jyesselm/seq_tools)
- [rna_secstruct](https://github.com/jyesselm/rna_secstruct)
- numpy, pandas, matplotlib, seaborn, scipy, lmfit

## Usage

```python
import pkamap
from pkamap.reference import NMR_REFERENCE, CONDITIONS, PROTONATION_SITES_ALL

# Load and process data
df = pkamap.process_json_files(json_files, codes=["C01HK"])
df_res = pkamap.generate_residue_dataframe(df)

# Fit pKa values
pka_fits = pkamap.calculate_pka_values(df_res)
pka_all = pkamap.add_structure_annotations(pka_fits, df)

# Filter for genuine protonation events
pka_filtered = pkamap.apply_quality_filters(pka_all)

# Summarize and compare to NMR
summary = pkamap.summarize_pka(pka_filtered, pka_all, conditions, PROTONATION_SITES_ALL)
stats = pkamap.compare_to_nmr(summary, NMR_REFERENCE, conditions)
```

See [`notebooks/analysis.ipynb`](notebooks/analysis.ipynb) for a complete workflow.

## Data

Sequencing run JSON files are not included in this repository. Place them in `data/` and update the paths in the analysis notebook. The NMR reference pKa values and motif definitions are bundled in `pkamap/reference.py`.

## Package Structure

```
pkamap/
├── __init__.py        # Public API
├── processing.py      # JSON loading, condition extraction, primer trimming
├── fitting.py         # Henderson-Hasselbalch fitting, structure annotation
├── filtering.py       # Quality-control filters
├── comparison.py      # Weighted pKa summary, NMR correlation plots
├── plotting.py        # Titration curve visualization
└── reference.py       # NMR reference data, motif definitions
```

## References

- Original PDMaP paper: https://pubs.acs.org/doi/10.1021/jacs.3c07736
- List of all papers used in CO1HK library: https://docs.google.com/document/d/1j1FX7y3rpGniN0VKyja1IOi7_4Fc90A0l-gywzMYwKg/edit?usp=sharing

