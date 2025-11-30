#!/usr/bin/env python3
"""
compute_lumi_from_csv_pxlhfoc.py

Purpose
-------
Compute the effective recorded integrated luminosity corresponding to the
(run, lumi) pairs present in the analysis ntuple, by matching them to
per-lumisection luminosity information from two official 2011 CSV sources:
the "pxl" (pixel-based) and "hfoc" (HF online calibration) luminosity records.

The script:
  1. Loads the (run, ls) set from a JSON file produced by `runls_to_json.py`,
     which contains ranges of lumisections for each run.
  2. Reads the per-lumisection luminosity from the 2011 `lumibyls` CSV files
     (pxl and hfoc), in the "record 1053" format.
  3. For each (run, ls) in the ntuple, uses pxl if available, otherwise hfoc.
  4. Sums the recorded luminosity in /ub and converts the total to /pb.

Inputs (defaults)
-----------------
- RUNLS_JSON : "runls_muonia.json"
    JSON with the (run, ls) coverage of the merged ntuple, in the standard
    `run: [ [ls_start, ls_end], ... ]` format.
- PXL_CSV    : "2011lumibyls_pxl.csv"
- HFOC_CSV   : "2011lumibyls_hfoc.csv"

  Both CSV files are expected to follow the CERN BRIL "record 1053" format
  used in the open data documentation, with at least the following columns:

    - '#run:fill'      : run and fill number, e.g. '160431:1615'
    - 'ls'             : lumisection range, e.g. '10:12'
    - 'recorded(/ub)'  : recorded luminosity in microbarn^{-1} (/ub)

Outputs
-------
- Prints to stdout:
    - How many (run, ls) pairs are matched using pxl vs hfoc.
    - How many (run, ls) are missing in both sources (orphans).
    - The total recorded luminosity in /ub and /pb.

Usage
-----
  python compute_lumi_from_csv_pxlhfoc.py

If needed, you can change `RUNLS_JSON`, `PXL_CSV`, and `HFOC_CSV` at the top
of the script to point to different files.

Notes
-----
- Units:
    The CSV files provide luminosity in /ub, so the script converts the total
    to /pb by dividing by 1e6.
- Priority:
    For each (run, ls) pair, pxl is used if available; otherwise hfoc is used.
"""

import json
import csv

# Input files
RUNLS_JSON = "runls_muonia.json"
PXL_CSV = "2011lumibyls_pxl.csv"
HFOC_CSV = "2011lumibyls_hfoc.csv"

# Column names according to the 2011 "record 1053" CSV format
RUNFILL_COL = "#run:fill"       # column with 'run:fill'
LS_COL = "ls"                   # column with 'ls_start:ls_end'
REC_COL = "recorded(/ub)"       # recorded luminosity in /ub


def load_runls_from_json(json_file):
    """
    Load (run, ls) pairs from the JSON (with LS ranges) and return a set of (run, ls).
    """
    with open(json_file) as f:
        data = json.load(f)

    runls_set = set()
    for run_str, ranges in data.items():
        run = int(run_str)
        for start, end in ranges:
            for ls in range(int(start), int(end) + 1):
                runls_set.add((run, ls))

    return runls_set


def load_lumi_by_ls(csv_file):
    """
    Load a per-lumisection luminosity CSV with header '#run:fill,ls,...'
    and return a dict[(run, ls)] = recorded_lumi_in_ub (float).

    If a row covers an LS range (e.g. '10:12'), the recorded luminosity is
    split uniformly among those lumisections.
    """
    lumi_dict = {}

    # Read all lines and look for the header line starting with '#run:fill'
    with open(csv_file) as f:
        lines = f.readlines()

    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith("#run:fill"):
            header_idx = i
            break

    if header_idx is None:
        raise RuntimeError(
            "No encontré una línea de cabecera que empiece por '#run:fill' en %s"
            % csv_file
        )

    # Use csv.DictReader starting from the header line
    reader = csv.DictReader(lines[header_idx:], skipinitialspace=True)

    # Check that the expected columns exist
    missing = [c for c in (RUNFILL_COL, LS_COL, REC_COL) if c not in reader.fieldnames]
    if missing:
        raise RuntimeError(
            "En el archivo %s no encontré las columnas: %s\nCabecera encontrada: %s"
            % (csv_file, missing, reader.fieldnames)
        )

    for row in reader:
        # Some rows may be incomplete; skip them
        runfill_raw = row.get(RUNFILL_COL)
        ls_raw = row.get(LS_COL)
        rec_raw = row.get(REC_COL)

        if runfill_raw is None or ls_raw is None or rec_raw is None:
            continue

        runfill = runfill_raw.strip()
        lsrange = ls_raw.strip()
        rec_str = rec_raw.strip()

        if not runfill or not lsrange or not rec_str:
            continue

        # runfill = '160431:1615' -> keep only the run
        try:
            run_str = runfill.split(":")[0]
            run = int(run_str)
        except ValueError:
            continue

        # lsrange = '19:19' or '10:12' -> LS range
        try:
            ls_start_str, ls_end_str = lsrange.split(":")
            ls_start = int(ls_start_str)
            ls_end = int(ls_end_str)
        except ValueError:
            # If not in start:end format, ignore
            continue

        try:
            rec = float(rec_str)
        except ValueError:
            continue

        # If the range spans more than one LS, split luminosity uniformly
        nls = ls_end - ls_start + 1
        if nls <= 0:
            continue
        rec_per_ls = rec / float(nls)

        for ls in range(ls_start, ls_end + 1):
            key = (run, ls)
            # If for some reason LS already exists, accumulate
            lumi_dict[key] = lumi_dict.get(key, 0.0) + rec_per_ls

    return lumi_dict


def main():
    # 1) Load (run, ls) pairs present in the merged ntuple
    runls_set = load_runls_from_json(RUNLS_JSON)
    print("Total de pares (run,ls) en", RUNLS_JSON, ":", len(runls_set))

    # 2) Load per-LS luminosity from pxl and hfoc CSV files
    print("Leyendo", PXL_CSV)
    pxl_lumi = load_lumi_by_ls(PXL_CSV)
    print("  Entradas (run,ls) pxl:", len(pxl_lumi))

    print("Leyendo", HFOC_CSV)
    hfoc_lumi = load_lumi_by_ls(HFOC_CSV)
    print("  Entradas (run,ls) hfoc:", len(hfoc_lumi))

    # 3) Sum luminosity with priority pxl > hfoc
    total_ub = 0.0
    used_pxl = 0
    used_hfoc = 0
    missing = []

    for run, ls in sorted(runls_set):
        key = (run, ls)
        if key in pxl_lumi:
            total_ub += pxl_lumi[key]
            used_pxl += 1
        elif key in hfoc_lumi:
            total_ub += hfoc_lumi[key]
            used_hfoc += 1
        else:
            # LS with no info in pxl nor hfoc; record as orphan
            missing.append(key)

    total_pb = total_ub / 1.0e6  # from /ub to /pb

    print("\n=== Resultado ===")
    print("LS con pxl usados   :", used_pxl)
    print("LS con hfoc usados  :", used_hfoc)
    print("LS sin info (órfanos):", len(missing))
    if missing:
        print("Ejemplo(s) de órfanos (primeros 10):", missing[:10])
    print("\nLuminosidad total registrada:")
    print("  %g /ub" % total_ub)
    print("  %g /pb" % total_pb)


if __name__ == "__main__":
    main()
