# pKaMaP

**RNA pKa Calculation Using DMS-MaPseq**

## Overview

pKaMaP fits the Henderson-Hasselbalch equation to reactivity from DMS-MaPseq experiments done across multiple pH levels to calculate the pKa values of individual RNA nucleotides.

The analysis pipeline:

1. Loads sequencing data from run JSON files (produced by [dreem](https://github.com/jyesselm/dreem) or equivalent) and extracts experimental conditions (e.g., temperature, pH, Mg²⁺) from naming conventions
2. Trims primer sequences and converts to per-residue dataframe
3. Fits pKa values using nonlinear least-squares minimization
4. Annotates each position with secondary-structure motif information via [rna_secstruct](https://github.com/jyesselm/rna_secstruct)
5. Filters for significant fits using quality thresholds based upon comparison with known protonation events
6. Data plotting and summarization
7. In-depth analysis comparing pKaMaP values to NMR results (for the C01HK 'pka comparison' library)

## Installation

```bash
pip install -e .
```

### Dependencies

- Python ≥ 3.10
- [seq_tools](https://github.com/jyesselm/seq_tools)
- [rna_secstruct](https://github.com/jyesselm/rna_secstruct)
- numpy, pandas, matplotlib, seaborn, scipy, lmfit

## Quick Start
from pkamap import process_json_files, generate_residue_dataframe
from pkamap import calculate_pka_values, add_structure_annotations

```python
df = process_json_files(["yourdata/sequencingrun.json", codes=["C01HK"])
df_res = generate_residue_dataframe(df, min_aligned=2000)

pka = add_structure_annotations(calculate_pka_values(df_res), df)
```

See [`analysis.ipynb`](analysis.ipynb) for a complete workflow.

## Input Data
 
### Sequencing JSON files
 
pKaMaP expects run JSON files where each row contains a construct with fields including `name`, `code`, `construct`, `sequence`, `structure`, `data` (per-nucleotide DMS reactivity values), `num_aligned`, and `sn`.
 
Experimental conditions are encoded in the `construct` field name:
 
| Convention | Example | Parsed value |
| :--- | :--- | :--- |
| pH | `pH750` | 7.50 |
| Buffer concentration | `50mM` in `pool_50mM_25C` | 50 mM (defaults to 300 mM if absent) |
| Temperature | `25C` in `pool_50mM_25C` | 25 °C (defaults to `"room"` if absent) |
| Mg²⁺ concentration | `0mMMgCl2` | 0 mM (defaults to 10 mM if absent) |

### 5′ primer CSV files
 
Each run JSON requires a corresponding primer file in the same directory. pKaMaP looks for `<run_name>_p5_sequences.csv` first, then falls back to a shared `p5_sequences.csv`. The CSV must contain a `sequence` column with 5′ primer sequences. 20 nucleotides are also trimmed from the 3′ end.
 
### Experimental conditions
 
Several pipeline functions require a `conditions` list that defines which experimental conditions to analyze. Each entry is a tuple of `(buffer_mM, temperature, Mg_mM, color, label)`:
 
```python
CONDITIONS = [
    (50,  25, 10, "green",  "50mM 25C"),
    (50,  25,  0, "purple", "50mM 25C 0mM Mg"),
]
```

## Quality Filtering
 
`apply_quality_filters()` removes fits that have DMS-reactivity responses to pH that have an absent, nonlinear, or noisy response to pH. Most nucleotides are not expected to show pH-dependent reactivity, so many residues are expected to be excluded from the filtered dataset. Default thresholds:
 
| Filter | Default | Description |
| :--- | :--- | :--- |
| `min_delta` | 0.00024 | Minimum high − low difference |
| `max_pka_error` | 0.56 | Maximum pKa fit uncertainty |
| `min_high` | 0.00067 | Minimum upper DMS reactivity |
| `max_high_error` | 0.02471 | Maximum upper fit uncertainty |
| `pka_range` | (4.0, 9.0) | Allowed pKa range |
 
All thresholds are adjustable — pass keyword arguments to `apply_quality_filters()` to override.

## pKa Summary
 
`summarize_pka()` computes an inverse-variance-weighted pKa for each unqiue residue (i.e., unique motif and secondary structure) that takes into account each repeat of the motif across the evaluated sequences. When multiple constructs contain the same motif and residues, their individual pKa fits are combined using `1/σ²` weighting, so higher-precision fits contribute more to the final value.
 
```python
from pkamap import summarize_pka
 
summary = summarize_pka(
    pka_filtered,       # filtered pKa values
    pka_all,            # unfiltered (used for denominator counts)
    CONDITIONS,         # experimental conditions
    protonation_sites,  # list of dicts with motif_sequence, nt_-1, nt_+0, nt_+1
)
```
 
The output DataFrame contains the weighted pKa, its error, and a construct fraction (`n` column, e.g. `"4/7"`) showing how many constructs passed filtering out of the total available for that motif and condition.


## NMR Comparison (C01HK Library)
 
The C01HK (pka comparison) library contains RNA motifs with published NMR pKa values. [`pkamap/nmr_reference_C01HK.py`](pkamap/nmr_reference_C01HK.py) bundles:
 
- **`NMR_REFERENCE`** — Published NMR pKa values and sources for included motifs
- **`MOTIFS_WITH_NMR`** / **`MOTIFS`** — Motif sequence lists
- **`PROTONATION_SITES`** — Specific residue positions with NMR-confirmed protonation events
 
Use `compare_to_nmr()` to generate correlation plots and R² statistics between pKaMaP and NMR values.
 
## Package Structure
 
```text
pkamap/
├── __init__.py              # Public API
├── processing.py            # JSON loading, condition extraction, primer trimming
├── fitting.py               # Henderson-Hasselbalch fitting, structure annotation
├── filtering.py             # Quality-control filters
├── comparison.py            # Weighted pKa summary, NMR correlation plots
├── plotting.py              # Titration curve visualization... visualize multiple conditions on the same plot!
└── nmr_reference_C01HK.py   # NMR reference data, motif definitions (C01HK library)
```
 
## References
 
- Characterizing Protonation-Coupled Conformational Ensembles in RNA via pH-Differential Mutational Profiling with DMS Probing
Edgar M. Faison, Amrita Nallathambi, and Qi Zhang
Journal of the American Chemical Society 2023 145 (34), 18773-18777
DOI: 10.1021/jacs.3c07736
- NMR reference values and sources included in pkamap/nmr_reference_C01HK.py
