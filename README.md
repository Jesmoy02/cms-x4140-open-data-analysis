# cms-x4140-open-data-analysis

Search for the χc₁(4140) in B⁺ → J/ψ φ K⁺ decays using CMS Run 2011A MuOnia Open Data.

This repository contains the analysis workflow developed to study the CMS Run 2011A MuOnia open data in a search for the χc₁(4140) candidate in B⁺ → J/ψ φ K⁺ decays. It comprises the main steps of the analysis, including event selection, reconstruction of intermediate states, and exploration of the B⁺ mass spectrum and the J/ψ φ invariant mass distribution.

Note: several scripts include the original local paths and directory structures used in this project to call other scripts, access input files, or store newly produced outputs. These paths are meant as a reference and may need to be adapted to the user’s local environment.

---

### Software environment

The analysis was performed in the following software environment:

- CMSSW_5_3_32 (CMS software release), running inside a Docker container.
- ROOT 5.32/00 (as provided by CMSSW_5_3_32).
- C++ analysis code and EDAnalyzers compiled within the CMSSW_5_3_32 area.
- Python 2.6 (inside CMSSW) for CMSSW-related tools and auxiliary scripts.
- Python 3.x (system / Anaconda) for the luminosity scripts under `luminosity/`.
- `brilcalc` (brilws) for luminosity calculations, using normtags `tag19v3` and `tag24v2`.
- Additional Python packages for luminosity handling: `uproot`, `pandas` (and standard scientific Python stack).

The repository assumes an existing CMSSW_5_3_32 working area with the custom `JPsiTrkTrkTrk` EDAnalyzer and the configuration files located under `cmssw_cfg/`. The detailed installation and setup of CMSSW and Docker are described in the thesis and in the official CMS / CERN Open Data documentation.

---

### Workflow overview

In short, the analysis workflow proceeds as follows:

1. **Select input AOD files**  
   Use the XRootD file index lists in `muonia_links/` to define the subset of the CMS Run2011A MuOnia AOD dataset used in this analysis. Each text file contains one AOD ROOT file URL per line.

2. **Download or process AOD files**  
   - **Option A (local AOD):**  
     Run `batch_scripts/01_download_muonia2011A_AOD_rootfiles.sh` to download AOD ROOT files into `root_files/`, then  
     `batch_scripts/02_run_jpsitrktrktrk_local_files.sh` to process the local AOD files with `cmssw_cfg/cmssw_cfg_jpsitrktrktrk.py`, producing skimmed JPsiTrkTrkTrk ntuples in `output_files_jpsitrktrktrk/`.  
   - **Option B (remote XRootD):**  
     Run `batch_scripts/03_run_jpsitrktrktrk_xrootd_remote.sh` to process the same AOD inputs directly over XRootD, without storing them locally. The script reads from `muonia_links/` and writes the skimmed ntuples to `output_files_jpsitrktrktrk/`.

3. **Merge skimmed outputs**  
   Run `batch_scripts/04_merge_jpsitrktrktrk_outputs.sh` to merge all per-file `*_jpsitrktrktrk.root` outputs into a single final ntuple:  
   `muonia_all_jpsitrktrktrk_merged.root`. This merged file is the starting point for the final physics plots and fits.

4. **Compute effective integrated luminosity**  
   Use the scripts under `luminosity/` to derive the effective integrated luminosity corresponding to the events present in `muonia_all_jpsitrktrktrk_merged.root`:
   - `dumpRunLumiFromTree.C` extracts the list of `(run, lumi)` pairs → `runls_muonia.txt`.  
   - `runls_to_json.py` converts this list into `runls_muonia.json`, in a format suitable for `brilcalc`.  
   - `brilcalc` and `compute_effective_lumi_from_csv.py` combine this JSON with the official 2011 `lumibyls` tables to obtain the effective recorded luminosity (in `/ub` and `/pb`).  
   The intermediate CSV and text products are stored in `luminosity/muonia_lumi/`.

5. **Perform final fits**  
   Use the ROOT/RooFit macros in `root_macros/` to extract the final physics observables:
   - `fit_Bplus_mass_dm_roofit.C` performs an unbinned fit to the B⁺ candidate mass spectrum.  
   - `fit_dm_twoRelBW_ps.C` performs an extended fit to the Δm distribution with a two-resonance relativistic Breit–Wigner plus phase-space model.  
   These fits use `muonia_all_jpsitrktrktrk_merged.root` as input and produce the B⁺ mass and Δm plots shown in the thesis.

---

### XRootD file index lists (`muonia_links/`)

This directory stores text files with XRootD URLs pointing to the CMS Run 2011A MuOnia AOD sample from the CERN Open Data portal. Each file

- follows the naming pattern  
  `CMS_Run2011A_MuOnia_AOD_12Oct2013-v1_<block>_file_index.txt`, and  
- contains one AOD ROOT file per line, e.g.  
  `root://eospublic.cern.ch//eos/opendata/cms/Run2011A/MuOnia/AOD/12Oct2013-v1/00000/00210445-2A44-E311-B69F-00304867918A.root`

These index files are used by the batch scripts in `batch_scripts/` to either:

- download subsets of the MuOnia AOD dataset locally into `root_files/`, or  
- run the `JPsiTrkTrkTrk` CMSSW configuration directly over remote inputs via XRootD.

The underlying dataset is the CMS Run 2011A MuOnia primary dataset in AOD format, as released through the CERN Open Data portal.

---

### CMSSW configuration

This directory contains the CMSSW configuration files used to run the custom `JPsiTrkTrkTrk` EDAnalyzer on the CMS Run 2011A MuOnia AOD sample.

- `cmssw_cfg/cmssw_cfg_jpsitrktrktrk.py`  
  CMSSW configuration used to run the `JPsiTrkTrkTrk` EDAnalyzer on the CMS Run 2011A MuOnia AOD sample. It handles command-line input and output files, selects the appropriate GlobalTag for data or MC, applies a good-lumi JSON when available, defines the J/ψ HLT pre-filter, configures the ROOT output via `TFileService`, clones and configures the `JPsiTrkTrkTrk` analyser, prints the primary selection parameters, and defines the processing `Path`.

- `cmssw_cfg/JPsiTrkTrkTrk_cfi.py`  
  Default configuration of the `JPsiTrkTrkTrk` EDAnalyzer. It defines the input collections (muons, tracks, primary vertices and HLT results), the list of J/ψ trigger paths, and the complete set of preselection and final selection cuts used in the thesis analysis. The module implements the six-step selection strategy for the B⁺ → J/ψ φ K⁺ candidates (J/ψ dimuon mass window, kaon pT/η/hits, φ mass window, vertex probabilities and Lxy significance). This `cfi` is cloned and slightly customised in `cmssw_cfg_jpsitrktrktrk.py`.

---

### Batch scripts

This directory contains helper shell scripts that automate the main stages of the workflow: downloading AOD files, running the `JPsiTrkTrkTrk` CMSSW configuration (locally or via XRootD), and merging the skimmed outputs into a single ROOT file.

- `batch_scripts/01_download_muonia2011A_AOD_rootfiles.sh`  
  Interactive script to download subsets of the CMS Run 2011A MuOnia AOD sample via XRootD. It reads file-index lists from `muonia_links/`, allows the user to select which list and which line range to process, downloads the corresponding AOD ROOT files into `root_files/` with a consistent naming scheme (`muonia_<slice>_<index>.root`), and records a timestamped log in `root_files/logs_downloads/` summarising successful, failed, and skipped downloads.

- `batch_scripts/02_run_jpsitrktrktrk_local_files.sh`  
  Batch runner for the `JPsiTrkTrkTrk` CMSSW configuration on locally stored ROOT files. It scans `root_files/`, allows the user to select a range of input files (or processes all of them when `SKIP_PROMPT=1`), and for each file calls `cmsRun` with `cmssw_cfg/cmssw_cfg_jpsitrktrktrk.py`, producing skimmed outputs in `output_files_jpsitrktrktrk/` with a `_jpsitrktrktrk.root` suffix and per-file log files. The script tracks the processing time per file, controls overwriting via the `OVERWRITE` flag, and prints a final summary of successful, failed, and skipped jobs.

- `batch_scripts/03_run_jpsitrktrktrk_xrootd_remote.sh`  
  Online batch runner for the `JPsiTrkTrkTrk` CMSSW configuration using remote AOD files accessed via XRootD. It reads XRootD URLs from `muonia_links/*_file_index.txt`, lets the user select an index file and a line range (or uses the full range when `SKIP_PROMPT=1`), and for each URL calls `cmsRun` with `cmssw_cfg/cmssw_cfg_jpsitrktrktrk.py`, writing the skimmed outputs to `output_files_jpsitrktrktrk/` with a `_jpsitrktrktrk.root` suffix and per-job logs in `output_files_jpsitrktrktrk/logs/`. The script checks output size, controls overwriting via `OVERWRITE`, and prints a final OK/FAIL/SKIP summary.

- `batch_scripts/04_merge_jpsitrktrktrk_outputs.sh`  
  Single-pass merger for all `*_jpsitrktrktrk.root` output files. It scans `output_files_jpsitrktrktrk/` for per-file JPsiTrkTrkTrk ROOT outputs (excluding any file that already contains `merged` in its name), writes a manifest with the full input list in `output_files_jpsitrktrktrk/merged/`, and runs a single `hadd -f` (ROOT 5.32, CMSSW_5_3_32 environment) to produce the final merged file `muonia_all_jpsitrktrktrk_merged.root`. A short preview of the inputs and a size check of the merged output are printed.

---

### Luminosity

This directory contains the scripts and auxiliary tables used to derive the effective integrated luminosity associated with the merged ntuple `muonia_all_jpsitrktrktrk_merged.root`. The resulting luminosity value is the one quoted in the final B⁺ mass and Δm plots in the thesis.

- `luminosity/dumpRunLumiFromTree.C`  
  ROOT macro that scans the final JPsiTrkTrkTrk TTree in `muonia_all_jpsitrktrktrk_merged.root` and extracts all unique `(run, lumi)` pairs that contribute to the selected candidates. It writes them to a plain-text file (by default `runls_muonia.txt`), one `run lumi` pair per line, sorted. This file is then used as input to the subsequent scripts (JSON / `brilcalc` / CSV workflow) to compute the effective integrated luminosity used in the thesis.

- `luminosity/runls_to_json.py`  
  Python helper that converts the flat list of `(run, ls)` pairs from `runls_muonia.txt` into a JSON file in the format required by `brilcalc` (per-run lists of contiguous LS ranges). The default output is `runls_muonia.json`, which is used by both the CSV-based and the `brilcalc`-based luminosity determinations.

- `luminosity/compute_effective_lumi_from_csv.py`  
  Python script that computes the effective recorded integrated luminosity for the analysed MuOnia sample. It reads the `(run, ls)` coverage from `runls_muonia.json`, matches each pair to the 2011 `lumibyls` CSV files (`2011lumibyls_pxl.csv` and `2011lumibyls_hfoc.csv`), and sums the recorded luminosity in `/ub` using a pxl > hfoc priority. The result is printed both in `/ub` and `/pb`.

All raw luminosity tables and `brilcalc` outputs used in this workflow are stored under `luminosity/muonia_lumi/`:

- `luminosity/muonia_lumi/muonia_lumi_brilcalc_tag19v3_lumibyls.csv`  
  Per-lumisection luminosity output from `brilcalc` for the MuOnia Run2011A selection, using the historical 2011 Open Data luminosity tag (e.g. `--normtag tag19v3 ... --byls`). Contains `run:fill`, LS ranges and recorded luminosity in `/ub`, restricted to the `(run, ls)` pairs listed in `runls_muonia.json`.

- `luminosity/muonia_lumi/muonia_lumi_brilcalc_tag24v2_lumibyls.csv`  
  Same as above, but produced with a more recent luminosity tag (e.g. `--normtag tag24v2 ... --byls`) and used to cross-check the stability of the effective integrated luminosity against the 19v3 configuration, LS by LS.

- `luminosity/muonia_lumi/muonia_lumi_brilcalc_tag19v3_summary.csv`  
  CSV summary produced with `brilcalc lumi` for the MuOnia Run2011A selection, using the historical 2011 Open Data normtag (`tag19v3`) and the JSON file `runls_muonia.json` as input. Each row corresponds to a single `run:fill` and includes the total number of lumisections (`nls`), the number of CMS lumisections (`ncms`), and the delivered/recorded integrated luminosity in `/pb`. This file is primarily a machine-readable counterpart to the plain-text `muonia_lumi_brilcalc_tag19v3_summary.txt` summary.

- `luminosity/muonia_lumi/muonia_lumi_brilcalc_tag24v2_summary.csv`  
  Same `brilcalc lumi` CSV summary as above, but obtained with the more recent luminosity normtag (`tag24v2`). It provides a per-run comparison of delivered and recorded luminosities in `/pb` between the legacy 19v3 and the updated 24v2 calibration, for precisely the same `(run, ls)` coverage used in the analysis.

- `luminosity/muonia_lumi/muonia_lumi_brilcalc_tag19v3_summary.txt`  
  Plain-text summary produced with `brilcalc lumi` for the MuOnia Run2011A selection, using the historical 2011 Open Data normtag (`tag19v3`) and the JSON file `runls_muonia.json` as input. The file contains an ASCII table with one row per `run:fill`, listing the number of lumisections (`nls`), the number of CMS lumisections (`ncms`), and the delivered/recorded integrated luminosity in `/pb`. It is primarily kept as a human-readable reference, complementing the machine-readable CSV file `muonia_lumi_brilcalc_tag19v3_summary.csv`.

- `luminosity/muonia_lumi/muonia_lumi_brilcalc_tag24v2_summary.txt`  
  Same `brilcalc lumi` plain-text summary as above, but computed with the more recent luminosity normtag (`tag24v2`). It provides a per-run comparison of delivered and recorded luminosities in `/pb` for the same `(run, ls)` coverage. It is used in the thesis to quote the effective luminosity corresponding to the updated calibration.

- `luminosity/muonia_lumi/runls_muonia.txt`  
  Plain-text list of all `(run, lumi)` pairs that contribute to the selected JPsiTrkTrkTrk candidates in `muonia_all_jpsitrktrktrk_merged.root`. It is produced by the ROOT macro `dumpRunLumiFromTree.C` and contains one `run lumi` pair per line (no weights, no duplication), already sorted. This file serves as the starting point for the effective-luminosity workflow.

- `luminosity/muonia_lumi/runls_muonia.json`  
  JSON version of the same `(run, ls)` information as `runls_muonia.txt`, but compressed into contiguous lumisection ranges for each run. It is generated by the helper script `runls_to_json.py` and is used as input both for the `brilcalc lumi` commands and for the CSV-based luminosity computation in `compute_effective_lumi_from_csv.py`.

---

### ROOT macros (fits)

This directory contains the ROOT/RooFit macros used to perform the final fits to the B⁺ mass spectrum and to the Δm distribution.

- `root_macros/fit_Bplus_mass_dm_roofit.C`  
  Unbinned (non-extended) RooFit fit of the B⁺ candidate mass (`cand_mass_fit`) using the final merged JPsiTrkTrkTrk ntuple (`muonia_all_jpsitrktrktrk_merged.root`). The macro reconstructs the J/ψ φ K⁺ system from the stored kinematics, applies the Δm window (`m(μμKK) − m(μμ)` in [1.008, 1.568] GeV) and an explicit B-mass window ([5.15, 5.45] GeV), and then fits the resulting `cand_mass_fit` distribution with a Gaussian signal plus a 2nd-order Chebyshev background. Two scenarios are produced: (i) floating mean and sigma, and (ii) mean fixed to the PDG B⁺ mass with sigma floating. The macro saves PNG/PDF plots and prints a summary of the fitted signal fraction, estimated signal yield, mean, and width.

- `root_macros/fit_dm_twoRelBW_ps.C`  
  RooFit model for the Δm spectrum, defined as \(m(\mu^+\mu^-K^+K^-) - m(\mu^+\mu^-)\), using the merged JPsiTrkTrkTrk ntuple. It reconstructs \(J/\psi \phi\) candidates, builds Δm in the window [1.008, 1.568] GeV, optionally applies an efficiency weight from `efficiency_vs_dm.root`, and performs an extended fit with two S-wave relativistic Breit–Wigner signals (mass-dependent widths, convolved with a common Gaussian resolution) plus a three-body phase-space background. It produces the plots `dm_twoRelBW_ps.{png,pdf}` and stores the rebinned histogram in `filtered_dm.root`.

---

### Final merged ROOT output

Due to GitHub size constraints, the final merged JPsiTrkTrkTrk ntuple used for the plots and fits in the thesis is stored externally and not tracked in this repository.

- File name: `muonia_all_jpsitrktrktrk_merged.root`  
- Size: ~161 MB  
- Location: [Download link (Google Drive)](https://drive.google.com/file/d/17Tk8BTXR8BlpdjVnIjPG_XVqiS_Au43Y/view?usp=drive_link)

This file is obtained by running `batch_scripts/04_merge_jpsitrktrktrk_outputs.sh` over all `*_jpsitrktrktrk.root` files produced in `output_files_jpsitrktrktrk/`. The detailed list of branches and the selection criteria applied in the analysis are described in the thesis.
