#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
runls_to_json.py

Purpose
-------
Convert a plain-text list of (run, lumiSection) pairs into a JSON file
with run → [ [ls_start, ls_end], ... ] ranges, suitable for use with
luminosity tools (e.g. brilcalc or CSV-based workflows).

Inputs (defaults)
-----------------
- infile  : "runls_muonia.txt"
    Text file produced by `dumpRunLumiFromTree.C`, containing one
    "run lumi" pair per line (whitespace-separated).

- outfile : "runls_muonia.json"
    JSON file with the following structure:
        {
          "160431": [[1, 10], [15, 18]],
          "160578": [[3, 7]],
          ...
        }

Usage
-----
  python runls_to_json.py

You can also modify `infile` and `outfile` at the top of this script
if you want to process a different run–lumi list.
"""

import json
import collections

infile = "runls_muonia.txt"
outfile = "runls_muonia.json"

runs = collections.OrderedDict()

with open(infile) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        run_str, ls_str = line.split()
        run = int(run_str)
        ls = int(ls_str)
        runs.setdefault(run, set()).add(ls)

# convertir a rangos [ls_start, ls_end]
json_dict = {}
for run in sorted(runs.keys()):
    ls_list = sorted(runs[run])
    ranges = []
    start = prev = ls_list[0]
    for ls in ls_list[1:]:
        if ls == prev + 1:
            prev = ls
        else:
            ranges.append([start, prev])
            start = prev = ls
    ranges.append([start, prev])
    json_dict[str(run)] = ranges

with open(outfile, "w") as out:
    json.dump(json_dict, out, sort_keys=True)

print("Escribí", outfile)
