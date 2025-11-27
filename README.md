# cms-x4140-open-data-analysis
Search for the χc1(4140) in B⁺→J/ψφK⁺ decays using CMS Run 2011A MuOnia Open Data.

This repository contains the analysis workflow developed to study CMS Run 2011A MuOnia open data in a search for the χc1(4140) candidate in B⁺→J/ψφK⁺ decays. It comprises the main steps of the analysis, including event selection, reconstruction of intermediate states, and exploration of the B⁺ mass spectrum and the J/ψφ invariant mass distribution.

`batch_scripts/01_download_muonia2011A_AOD_rootfiles.sh`  
Interactive script to download subsets of the CMS Run 2011A MuOnia AOD sample via XRootD.  It reads file-index lists from `muonia_links/`, allows the user to select which list, and which line range to process, downloads the corresponding ROOT files into `root_files/` with a consistent naming scheme (`muonia_<slice>_<index>.root`), and records a timestamped log in `root_files/logs_downloads/` summarising successful, failed, and skipped downloads.
