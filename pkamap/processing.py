"""
Data processing: Take sequencing results form JSON files and extract experimental conditions from construct names,
trim primers, retain sequences with desired conditions, and generate the per-residue dataframe
"""

import pathlib
import re

import numpy as np
import pandas as pd
import seq_tools


# Condition extraction from construct names

def extract_ph(construct: str) -> float | None:
    """Extract pH from construct name (e.g., 'pH750' -> 7.50)."""
    idx = construct.find("pH")
    if idx != -1:
        try:
            return float(construct[idx + 2 : idx + 5]) / 100
        except ValueError:
            return None
    return None


def extract_buffer_concentration(name: str) -> int:
    """Extract buffer concentration in mM (e.g., '50mM' from 'pool_50mM_25C').

    Returns 300 mM if no explicit value is found (legacy naming convention).
    """
    match = re.search(r"_(\d+)mM_", name)
    return int(match.group(1)) if match else 300


def extract_temperature(name: str) -> int | str:
    """Extract temperature in °C (e.g., 25 from '50mM_25C_0mM').

    Returns ``'room'`` when the construct name does not encode a temperature.
    """
    match = re.search(r"_(\d+)C_", name)
    return int(match.group(1)) if match else "room"


def extract_magnesium(name: str) -> int:
    """Extract MgCl₂ concentration in mM (e.g., 0 from '0mMMgCl2').

    Returns 10 mM if no explicit value is found (legacy naming convention).
    """
    match = re.search(r"_(\d+)mMMgCl2", name)
    return int(match.group(1)) if match else 10


# Primer trimming

def _trim(df: pd.DataFrame, start: int, end: int) -> pd.DataFrame:
    """Trim sequence, structure, and data columns by *start* and *end* positions."""

    def _trim_str(col, s, e):
        if s == 0 and e == 0:
            return col
        return col.str[s:] if e == 0 else col.str[s:-e]

    def _trim_list(x, s, e):
        if x is None:
            return None
        try:
            return x[s:-e] if e != 0 else x[s:]
        except (TypeError, IndexError):
            return x

    df = df.copy()
    for col in ["sequence", "structure"]:
        if col in df.columns:
            df[col] = _trim_str(df[col], start, end)
    if "data" in df.columns:
        df["data"] = df["data"].apply(lambda x: _trim_list(x, start, end))
    return df


def _infer_p5_csv(json_path: pathlib.Path) -> pathlib.Path:
    """Locate the 5ʹ primer CSV for a given run JSON.

    Searches for ``<run_name>_p5_sequences.csv`` first, then falls back to
    a shared ``p5_sequences.csv`` in the same directory.
    """
    json_dir = json_path.parent
    candidate = json_dir / f"{json_path.stem}_p5_sequences.csv"
    if candidate.exists():
        return candidate
    default = json_dir / "p5_sequences.csv"
    if default.exists():
        return default
    raise FileNotFoundError(f"No p5 CSV found for {json_path}")


def _trim_p5_and_p3(df: pd.DataFrame, p5_csv_path: pathlib.Path) -> pd.DataFrame:
    """Remove 5ʹ primer and 20 nt from the 3ʹ end."""
    df_p5 = pd.read_csv(p5_csv_path)
    for p5_seq in df_p5["sequence"]:
        if seq_tools.has_5p_sequence(df, p5_seq):
            return _trim(df, len(p5_seq), 20)
    raise ValueError("No matching 5ʹ primer sequence found")


# Filtering

def _process_raw_df(
    df: pd.DataFrame, code: str, p5_csv_path: pathlib.Path
) -> pd.DataFrame:
    """Filter for *code*, extract conditions, trim primers, convert T→U."""
    df = df.query("code == @code").copy()
    if df.empty:
        return df

    df["ph"] = df["construct"].apply(extract_ph)
    df["buffer_conc"] = df["construct"].apply(extract_buffer_concentration)
    df["temperature"] = df["construct"].apply(extract_temperature)
    df["mg_conc"] = df["construct"].apply(extract_magnesium)
    df = df.dropna(subset=["ph"])

    try:
        df = _trim_p5_and_p3(df, p5_csv_path)
    except ValueError:
        pass

    df["sequence"] = df["sequence"].str.replace("T", "U")
    return df


# Public API

def process_json_files(
    json_files: list[str | pathlib.Path],
    codes: list[str],
) -> pd.DataFrame:
    """Load and process one or more sequencing-run JSON files.

    Parameters
    ----------
    json_files : list of str or Path
        Paths to run JSON files produced by the sequencing pipeline.
    codes : list of str
        Library construct codes to extract (e.g. ``["C01HK"]``).

    Returns
    -------
    pd.DataFrame
        Concatenated, primer-trimmed DataFrame with columns for pH, buffer
        concentration, temperature, and Mg²⁺ concentration.
    """
    frames = []
    for json_file in json_files:
        json_path = pathlib.Path(json_file)
        if not json_path.exists():
            print(f"Warning: {json_path} not found, skipping")
            continue
        df_raw = pd.read_json(json_path)
        p5_csv = _infer_p5_csv(json_path)
        for code in codes:
            result = _process_raw_df(df_raw, code, p5_csv)
            if not result.empty:
                print(f"{json_path.stem} | {code} → {len(result)} rows")
                frames.append(result)
    if not frames:
        raise RuntimeError("No data loaded — check paths and codes")
    return pd.concat(frames, ignore_index=True)


def generate_residue_dataframe(
    df: pd.DataFrame,
    min_aligned: int = 2000,
    include_all_bases: bool = True,
) -> pd.DataFrame:
    """Expand construct-level data to per-residue rows.

    Parameters
    ----------
    df : pd.DataFrame
        Construct-level DataFrame (output of :func:`process_json_files`).
    min_aligned : int
        Minimum number of aligned reads to keep a construct.
    include_all_bases : bool
        If ``False``, exclude U and G residues (DMS primarily modifies A and C).
    """
    records = []
    for _, row in df.iterrows():
        seq, ss = row["sequence"], row["structure"]
        for j, d in enumerate(row["data"]):
            if not include_all_bases and seq[j] in ("U", "G"):
                continue
            records.append(
                {
                    "seq": seq[j],
                    "ss": ss[j],
                    "pos": j + 1,
                    "data": d,
                    "ph": row["ph"],
                    "buffer_conc": row["buffer_conc"],
                    "temperature": row["temperature"],
                    "mg_conc": row["mg_conc"],
                    "name": row["name"],
                    "construct": row["construct"],
                    "sn": row["sn"],
                    "num_aligned": row["num_aligned"],
                }
            )
    df_res = pd.DataFrame(records)
    df_res = df_res[df_res["num_aligned"] >= min_aligned]
    return df_res