# cms-x4140-open-data-analysis
Search for the χc1(4140) in B⁺→J/ψφK⁺ decays using CMS Run 2011A MuOnia Open Data.

This repository contains the analysis workflow developed to study CMS Run 2011A MuOnia open data in a search for the χc₁(4140) candidate in B⁺→J/ψφK⁺ decays. It comprises the main steps of the analysis, including event selection, reconstruction of intermediate states, and exploration of the B⁺ mass spectrum and the J/ψφ invariant mass distribution.

- `batch_scripts/01_download_muonia2011A_AOD_rootfiles.sh`  
Interactive script to download subsets of the CMS Run 2011A MuOnia AOD sample via XRootD.  It reads file-index lists from `muonia_links/`, allows the user to select which list, and which line range to process, downloads the corresponding ROOT files into `root_files/` with a consistent naming scheme (`muonia_<slice>_<index>.root`), and records a timestamped log in `root_files/logs_downloads/` summarising successful, failed, and skipped downloads.

- `batch_scripts/02_run_jpsitrktrktrk_local_files.sh`  
Batch runner for the `JPsiTrkTrkTrk` CMSSW configuration on locally stored ROOT files. It scans `root_files/`, allows the user to select a range of input files (or processes all of them when `SKIP_PROMPT=1`), and for each file calls `cmsRun` with `poet_cfg_jpsitrktrktrk.py`, producing skimmed outputs in `output_files_jpsitrktrktrk/`  with a `_jpsitrktrktrk.root` suffix and per-file log files. The script tracks the processing time per file, controls overwriting via the `OVERWRITE` flag, and prints a final summary of successful, failed, and skipped jobs.

- `batch_scripts/03_run_jpsitrktrktrk_xrootd_remote.sh`  
Online batch runner for the `JPsiTrkTrkTrk` CMSSW configuration using remote AOD files accessed via XRootD. It reads XRootD URLs from `muonia_links/*_file_index.txt`, lets the user select an index file and a line range (or uses the full range when `SKIP_PROMPT=1`), and for each URL calls `cmsRun` with `poet_cfg_jpsitrktrktrk.py`, writing the skimmed outputs to `output_files_jpsitrktrktrk/` with a `_jpsitrktrktrk.root` suffix and per-job logs in `output_files_jpsitrktrktrk/logs/`. The script checks output size, controls overwriting via `OVERWRITE`, and prints a final OK/FAIL/SKIP summary.

- `cmssw_cfg/cmssw_cfg_jpsitrktrktrk.py`  
CMSSW configuration used to run the `JPsiTrkTrkTrk` EDAnalyzer on the CMS Run 2011A MuOnia AOD sample. It handles command-line input/output files, selects the appropriate GlobalTag for data or MC, applies a good-lumi JSON when available, defines the J/ψ HLT pre-filter, configures the ROOT output via `TFileService`, clones and configures the `JPsiTrkTrkTrk` analyser, prints the primary selection parameters, and defines the processing `Path`.

- `cmssw_cfg/JPsiTrkTrkTrk_cfi.py`
Default configuration of the `JPsiTrkTrkTrk` EDAnalyzer. It defines the input collections (muons, tracks, primary vertices and HLT results), the list of J/ψ trigger paths, and the complete set of preselection and final selection cuts used in the thesis analysis. The module implements the six-step selection strategy for the B⁺ → J/ψ φ K⁺ candidates (J/ψ dimuon mass window, kaon pT/η/hits, φ mass window, vertex probabilities and Lxy significance). This cfi is cloned and slightly customised in `cmssw_cfg_jpsitrktrktrk.py`.

- `batch_scripts/04_merge_jpsitrktrktrk_outputs.sh`  
Single-pass merger for all `*_jpsitrktrktrk.root` output files. It scans `output_files_jpsitrktrktrk/` for per-file JPsiTrkTrkTrk ROOT outputs (excluding any file that already contains `merged` in its name), writes a manifest with the full input list in `output_files_jpsitrktrktrk/merged/`, and runs a single `hadd -f` (ROOT 5.32, CMSSW_5_3_32 environment) to produce the final merged file `muonia_all_jpsitrktrktrk_merged.root`. A short preview of the inputs and a size check of the merged output are printed.

### ROOT outputs
- `root_outputs/muonia_all_jpsitrktrktrk_merged.root`
Final merged JPsiTrkTrkTrk ntuple used for the analysis plots and fits in the thesis (see the thesis for the detailed list of branches and selection criteria).

