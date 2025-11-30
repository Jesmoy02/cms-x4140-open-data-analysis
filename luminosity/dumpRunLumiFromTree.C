#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <set>
#include <utility>

// dumpRunLumiFromTree.C
//
// Purpose:
//   Scan a CMS analysis TTree and extract the unique (run, lumi)
//   pairs corresponding to the events that survive the final selection.
//   These (run, lumi) pairs are then used as input to compute the
//   effective integrated luminosity (e.g. with brilcalc or CSV-based tools).
//
// Inputs (defaults):
//   - infile    : "muonia_all_jpsitrktrktrk_merged.root"
//       Merged JPsiTrkTrkTrk ntuple produced after the full selection.
//   - treename  : "myJPsiTrkTrkTrk/JPsiTrkTrkTrkTree"
//       Path to the analysis TTree inside the ROOT file.
//   - runBranch : name of the run branch (UInt_t).
//   - lumiBranch: name of the lumi section branch (UInt_t).
//   - outfile   : "runls_muonia.txt"
//       Output ASCII file with one "run lumi" pair per line.
//
// Output:
//   - A plain text file (by default `runls_muonia.txt`) containing
//     all unique (run, lumi) pairs, sorted. This file is the basis
//     for building JSONs or run–lumi lists used in the effective
//     luminosity calculation for the thesis analysis.
//
// Usage (ROOT batch):
//   root -l -b -q 'dumpRunLumiFromTree.C()'
//   root -l -b -q 'dumpRunLumiFromTree.C("myMerged.root",
//                                        "myTreePath",
//                                        "run", "lumi",
//                                        "runls_output.txt")'

void dumpRunLumiFromTree(const char* infile  = "muonia_all_jpsitrktrktrk_merged.root",
                         const char* treename = "myJPsiTrkTrkTrk/JPsiTrkTrkTrkTree",
                         const char* runBranch  = "run",
                         const char* lumiBranch = "lumi",
                         const char* outfile = "runls_muonia.txt")
{
  TFile *f = TFile::Open(infile);
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: no pude abrir el archivo " << infile << std::endl;
    return;
  }

  TTree *t = dynamic_cast<TTree*>(f->Get(treename));
  if (!t) {
    std::cerr << "ERROR: no encontré el árbol " << treename << std::endl;
    f->Close();
    return;
  }

  UInt_t run = 0;
  UInt_t lumi = 0;

  if (t->GetBranch(runBranch) == 0 || t->GetBranch(lumiBranch) == 0) {
    std::cerr << "ERROR: no encontré las ramas " << runBranch
              << " o " << lumiBranch << std::endl;
    f->Close();
    return;
  }

  t->SetBranchAddress(runBranch, &run);
  t->SetBranchAddress(lumiBranch, &lumi);

  std::set< std::pair<UInt_t,UInt_t> > runlsSet;

  Long64_t nentries = t->GetEntries();
  std::cout << "Leyendo " << nentries << " entradas..." << std::endl;

  for (Long64_t i = 0; i < nentries; ++i) {
    t->GetEntry(i);
    runlsSet.insert(std::make_pair(run, lumi));
  }

  std::ofstream out(outfile);
  if (!out.is_open()) {
    std::cerr << "ERROR: no pude abrir el archivo de salida " << outfile << std::endl;
    f->Close();
    return;
  }

  for (std::set< std::pair<UInt_t,UInt_t> >::const_iterator it = runlsSet.begin();
       it != runlsSet.end(); ++it) {
    out << it->first << " " << it->second << "\n";  // run lumi
  }

  out.close();
  f->Close();

  std::cout << "Guardados " << runlsSet.size()
            << " pares run-lumi en " << outfile << std::endl;
}
