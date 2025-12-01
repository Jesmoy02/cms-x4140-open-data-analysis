// fit_dm_twoRelBW_ps.C
//
// Purpose:
//   Build the Δm spectrum,
//     Δm = m(μ⁺μ⁻K⁺K⁻) − m(μ⁺μ⁻),
//   from the merged JPsiTrkTrkTrk ntuple and perform an extended RooFit
//   with two signal hypotheses plus a three-body phase-space background.
//
// Model:
//   - Observable:
//       dm  ≡  m(μ⁺μ⁻K⁺K⁻) − m(μ⁺μ⁻), in a Δm window [dm_min, dm_max].
//   - Kinematics:
//       J/ψ and φ are reconstructed from the JPsiTrkTrkTrk tree;
//       dm is computed event-by-event using 4-vectors built from
//       (pT, η, φ, m).
//   - Signals (two-component):
//       • sig1: relativistic S-wave Breit–Wigner with mass-dependent width
//               Γ₁(m), convolved with a Gaussian resolution G(dm; 0, σ);
//       • sig2: same structure with independent m₀₂ and Γ₀₂, also
//               convolved with the same Gaussian resolution.
//   - Background:
//       • ps: three-body phase-space inspired shape
//             ~ (dm − dm_min)^{α_L} (dm_max − dm)^{α_U}.
//   - Extended fit:
//       model(dm) = N₁·sig1(dm) + N₂·sig2(dm) + N_B·ps(dm),
//       with (N₁, N₂, N_B) floating yields.
//
// Inputs (defaults):
//   - infile   : "muonia_all_jpsitrktrktrk_merged.root"
//   - treePath : "myJPsiTrkTrkTrk/JPsiTrkTrkTrkTree"
//   - dm_min   : 1.008 GeV (lower edge of Δm window, as in the CMS note)
//   - dm_max   : 1.568 GeV (upper edge of Δm window)
//   - binw     : 0.020 GeV (20 MeV bin width for plotting)
//
// Optional inputs:
//   - efficiency_vs_dm.root (if present in the working directory):
//       • histogram "heff" is interpreted as ε(dm);
//       • each candidate is given a weight w = 1/ε(dm), truncated at w ≤ 50;
//       • the RooDataSet is built with WeightVar(w) and a 1D histogram
//         "hdata" is filled with the same weights.
//
// Output:
//   - Plots:
//       • dm_twoRelBW_ps.png
//       • dm_twoRelBW_ps.pdf
//   - ROOT file:
//       • filtered_dm.root containing the binned Δm histogram (TH1F hdata).
//   - Console summary:
//       • Best-fit values and uncertainties for (m₀₁, Γ₀₁, N₁),
//         (m₀₂, Γ₀₂, N₂), resolution σ, and N_B.
//
// Requirements:
//   - ROOT ≥ 6 with RooFit enabled.
//   - The merged JPsiTrkTrkTrk file and tree must be present under the
//     names specified by `infile` and `treePath`.
//   - If efficiency weighting is desired, the file efficiency_vs_dm.root
//     with TH1F "heff" should be available in the current directory.

#include <vector>
#include <cmath>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLorentzVector.h>
#include <TH1F.h>
#include <TLegend.h>

#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"

using namespace RooFit;

static const double M_K     = 0.493677;   // GeV
static const double M_JPSI  = 3.096916;   // GeV (PDG)
static const double M_PHI   = 1.019461;   // GeV (PDG)

// 4-vector from (pt, eta, phi, m)
static TLorentzVector p4(float pt,float eta,float phi,double m){
  TLorentzVector v; const double px=pt*std::cos(phi), py=pt*std::sin(phi), pz=pt*std::sinh(eta);
  v.SetPxPyPzE(px,py,pz, std::sqrt(px*px+py*py+pz*pz+m*m)); return v;
}

// Choose the KK pair that forms the φ using mKK_min (with tolerance)
static void chooseKK(double m12,double m13,double m23,double mmin,int& i,int& j){
  const double tol=0.005; i=j=-1;
  if (std::fabs(m12-mmin)<tol){i=1;j=2;return;}
  if (std::fabs(m13-mmin)<tol){i=1;j=3;return;}
  if (std::fabs(m23-mmin)<tol){i=2;j=3;return;}
  const double d12=std::fabs(m12-mmin), d13=std::fabs(m13-mmin), d23=std::fabs(m23-mmin);
  if (d12<=d13 && d12<=d23){i=1;j=2;} else if (d13<=d12 && d13<=d23){i=1;j=3;} else {i=2;j=3;}
}

// ----------------------------- MAIN ---------------------------------
void fit_dm_twoRelBW_ps(const char* infile="muonia_all_jpsitrktrktrk_merged.root",
                        const char* treePath="myJPsiTrkTrkTrk/JPsiTrkTrkTrkTree",
                        double dm_min=1.008, double dm_max=1.568,  // CMS note window
                        double binw=0.020)                          // 20 MeV
{
  // --- Open ROOT file & tree
  TFile* f=TFile::Open(infile);
  if(!f||f->IsZombie()){ printf("ERROR: no pude abrir %s\n", infile); return; }
  TTree* t=(TTree*)f->Get(treePath);
  if(!t){ printf("ERROR: no encontré %s\n", treePath); f->Close(); return; }

  // Small helper: set vector<float>* with alternative branch names (a or b)
  auto setVec=[&](const char* a, const char* b, std::vector<float>** p)->bool{
    if (t->GetBranch(a)) { t->SetBranchAddress(a,p); return true; }
    if (b && t->GetBranch(b)) { t->SetBranchAddress(b,p); return true; }
    printf("ERROR: no existe rama '%s'%s%s\n", a, b?" ni '":"", b?b:"");
    return false;
  };

  // Branch pointers
  std::vector<float> *dimu_pt_raw=0, *dimu_eta_raw=0, *dimu_phi_raw=0, *dimu_mass_raw=0;
  std::vector<float> *trk1_pt=0, *trk1_eta=0, *trk1_phi=0;
  std::vector<float> *trk2_pt=0, *trk2_eta=0, *trk2_phi=0;
  std::vector<float> *trk3_pt=0, *trk3_eta=0, *trk3_phi=0;
  std::vector<float> *m12=0, *m13=0, *m23=0, *mmin=0;

  t->SetBranchAddress("dimu_pt_raw",   &dimu_pt_raw);
  t->SetBranchAddress("dimu_eta_raw",  &dimu_eta_raw);
  if (!setVec("dimu_phl_raw","dimu_phi_raw",&dimu_phi_raw)) { f->Close(); return; }
  t->SetBranchAddress("dimu_mass_raw", &dimu_mass_raw);
  t->SetBranchAddress("trk1_pt",&trk1_pt); t->SetBranchAddress("trk1_eta",&trk1_eta); if(!setVec("trk1_phl","trk1_phi",&trk1_phi)) { f->Close(); return; }
  t->SetBranchAddress("trk2_pt",&trk2_pt); t->SetBranchAddress("trk2_eta",&trk2_eta); if(!setVec("trk2_phl","trk2_phi",&trk2_phi)) { f->Close(); return; }
  t->SetBranchAddress("trk3_pt",&trk3_pt); t->SetBranchAddress("trk3_eta",&trk3_eta); if(!setVec("trk3_phl","trk3_phi",&trk3_phi)) { f->Close(); return; }
  t->SetBranchAddress("mKK_12",&m12);
  t->SetBranchAddress("mKK_13",&m13);
  t->SetBranchAddress("mKK_23",&m23);
  t->SetBranchAddress("mKK_min",&mmin);

  // RooFit dataset (with optional efficiency weights)
  RooRealVar dm("dm","#Delta m", dm_min, dm_max, "GeV");
  RooRealVar w("w","weight", 1.0, 0.0, 1e6);
  RooDataSet data("data","data", RooArgSet(dm,w), WeightVar(w));

  // Optional efficiency: heff(dm) from efficiency_vs_dm.root
  TH1F* heff=nullptr;
  { TFile fe("efficiency_vs_dm.root");
    if (!fe.IsZombie()) { fe.GetObject("heff", heff); if (heff) heff->SetDirectory(nullptr); }
  }

  // 1D histogram for reference / export
  const int nbins = std::max(1,int(std::floor((dm_max-dm_min)/binw + 0.5)));
  TH1F hdata("hdata",";m_{#mu#mu K^{+}K^{-}}-m_{#mu#mu}  [GeV];N / 20 MeV", nbins, dm_min, dm_max);

  // Fill Δm with optional 1/ε(dm) weighting
  const Long64_t nEv = t->GetEntries();
  Long64_t nCand=0;
  for (Long64_t ie=0; ie<nEv; ++ie){
    t->GetEntry(ie);
    size_t n = mmin->size();
    auto shrink=[&](const std::vector<float>* v){ n = std::min(n, v->size()); };
    shrink(dimu_pt_raw); shrink(dimu_eta_raw); shrink(dimu_phi_raw); shrink(dimu_mass_raw);
    shrink(trk1_pt); shrink(trk1_eta); shrink(trk1_phi);
    shrink(trk2_pt); shrink(trk2_eta); shrink(trk2_phi);
    shrink(trk3_pt); shrink(trk3_eta); shrink(trk3_phi);
    shrink(m12); shrink(m13); shrink(m23);

    for (size_t ic=0; ic<n; ++ic){
      TLorentzVector pdimu = p4(dimu_pt_raw->at(ic), dimu_eta_raw->at(ic), dimu_phi_raw->at(ic), dimu_mass_raw->at(ic));
      TLorentzVector pk1   = p4(trk1_pt->at(ic), trk1_eta->at(ic), trk1_phi->at(ic), M_K);
      TLorentzVector pk2   = p4(trk2_pt->at(ic), trk2_eta->at(ic), trk2_phi->at(ic), M_K);
      TLorentzVector pk3   = p4(trk3_pt->at(ic), trk3_eta->at(ic), trk3_phi->at(ic), M_K);
      int a=-1,b=-1; chooseKK(m12->at(ic), m13->at(ic), m23->at(ic), mmin->at(ic), a, b);
      TLorentzVector pphi = (a==1&&b==2)? pk1+pk2 : (a==1&&b==3? pk1+pk3 : pk2+pk3);

      const double dmv = (pdimu+pphi).M() - pdimu.M();
      if (dmv<=dm_min || dmv>=dm_max) continue;

      double ww = 1.0;
      if (heff){
        int bin = heff->FindBin(dmv);
        double eps = std::max(1e-3, std::min(0.999, (double)heff->GetBinContent(bin)));
        ww = 1.0/eps; if (ww>50) ww=50;
      }
      dm.setVal(dmv); w.setVal(ww);
      data.add(RooArgSet(dm,w), ww);
      hdata.Fill(dmv, ww);
      ++nCand;
    }
  }
  printf("Candidatos Δm usados: %lld (bins=%d)\n", nCand, nbins);

  // ===== Model definition =====
  // mX = dm + mJpsi  (X(→J/psi phi) mass as function of Δm)
  RooConstVar mJ("mJ","mJ", M_JPSI);
  RooFormulaVar mX("mX","@0+@1", RooArgList(dm,mJ));
  RooConstVar mPhi("mPhi","mPhi", M_PHI);

  // ----- Peak 1 (S-wave relativistic BW with mass-dependent width Γ(m)) -----
  RooRealVar m01("m01","m0_1", 1.051, 1.000, 1.100);
  RooRealVar g01("g01","Gamma0_1", 0.030, 0.005, 0.090);
  RooFormulaVar MX1("MX1","@0+@1", RooArgList(m01,mJ));
  RooFormulaVar q1("q1",
    "sqrt(max(0.0,(@0*@0-pow(@1+@2,2))*(@0*@0-pow(@1-@2,2))))/(2*@0)",
    RooArgList(mX,mJ,mPhi));
  RooFormulaVar q01("q01",
    "sqrt(max(0.0,(@0*@0-pow(@1+@2,2))*(@0*@0-pow(@1-@2,2))))/(2*@0)",
    RooArgList(MX1,mJ,mPhi));
  RooFormulaVar Gam1("Gam1","@0*(@1/@2)*(@3/@4)", RooArgList(g01,q1,q01, MX1, mX));
  RooGenericPdf rbw1("rbw1",
    "1.0/((@0*@0-@1*@1)*(@0*@0-@1*@1) + @1*@1*@2*@2)",
    RooArgList(mX, MX1, Gam1));

  // Common resolution in Δm (Gaussian centered at 0)
  RooRealVar sigma("sigma","resolution", 0.010, 0.006, 0.030);
  RooGaussian G("G","G", dm, RooConst(0.0), sigma);

  RooFFTConvPdf sig1("sig1","rbw1 (*) G", dm, rbw1, G);
  sig1.setBufferFraction(1.0);

  // ----- Peak 2 (same structure, different m0 and Γ0) -----
  RooRealVar m02("m02","m0_2", 1.220, 1.160, 1.300);
  RooRealVar g02("g02","Gamma0_2",0.035, 0.005, 0.100);
  RooFormulaVar MX2("MX2","@0+@1", RooArgList(m02,mJ));
  RooFormulaVar q2("q2",
    "sqrt(max(0.0,(@0*@0-pow(@1+@2,2))*(@0*@0-pow(@1-@2,2))))/(2*@0)",
    RooArgList(mX,mJ,mPhi));
  RooFormulaVar q02("q02",
    "sqrt(max(0.0,(@0*@0-pow(@1+@2,2))*(@0*@0-pow(@1-@2,2))))/(2*@0)",
    RooArgList(MX2,mJ,mPhi));
  RooFormulaVar Gam2("Gam2","@0*(@1/@2)*(@3/@4)", RooArgList(g02,q2,q02, MX2, mX));
  RooGenericPdf rbw2("rbw2",
    "1.0/((@0*@0-@1*@1)*(@0*@0-@1*@1) + @1*@1*@2*@2)",
    RooArgList(mX, MX2, Gam2));

  RooFFTConvPdf sig2("sig2","rbw2 (*) G", dm, rbw2, G);
  sig2.setBufferFraction(1.0);

  // ----- Background: three-body phase space (shape only) -----
  RooConstVar dmL("dmL","dmL", dm_min);
  RooConstVar dmU("dmU","dmU", dm_max);
  RooRealVar aL("aL","alpha_L", 1.0, 0.0, 6.0);
  RooRealVar aU("aU","alpha_U", 1.0, 0.0, 6.0);
  RooGenericPdf ps("ps","pow(@0-@1,@2)*pow(@3-@0,@4)", RooArgList(dm,dmL,aL,dmU,aU));

  // Extended mixture: two signals + PS background
  RooRealVar N1("N1","yield1",  300, 0, 1e6);
  RooRealVar N2("N2","yield2",  300, 0, 1e6);
  RooRealVar NB("NB","yieldB",  500, 0, 1e7);
  RooAddPdf model("model","sig1+sig2+ps", RooArgList(sig1,sig2,ps), RooArgList(N1,N2,NB));

  // Fit (quiet, with SumW2 errors)
  model.fitTo(data, PrintLevel(-1), SumW2Error(kTRUE));

  // ------------------- Plot with CMS-like style -------------------
  gStyle->SetOptStat(0);
  RooPlot* fr = dm.frame(Bins(nbins));
  fr->SetTitle("");
  fr->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}K^{+}K^{-}} - m_{#mu^{+}#mu^{-}}  [GeV]");
  fr->GetYaxis()->SetTitle("N / 20 MeV");
  fr->GetXaxis()->SetTitleSize(0.045);
  fr->GetYaxis()->SetTitleSize(0.045);
  fr->GetXaxis()->SetLabelSize(0.040);
  fr->GetYaxis()->SetLabelSize(0.040);

  // Data points
  data.plotOn(fr, DataError(RooAbsData::SumW2),
              Name("dataPoints"),
              MarkerStyle(kFullCircle), MarkerSize(1.0),
              LineColor(kBlack));

  // Total fit (thick red line)
  model.plotOn(fr, LineColor(kRed), LineWidth(3), Name("totCurve"));
  // Three-body PS (blue dotted)
  model.plotOn(fr, Components(ps.GetName()),
               LineColor(kBlue+1), LineStyle(kDotted), LineWidth(3),
               Name("psCurve"));

  TCanvas c("c_dm2","Two signal-hypothesis (relBW ⊗ G)", 900, 720);
  c.SetLeftMargin(0.12); c.SetRightMargin(0.04);
  c.SetTopMargin(0.08);  c.SetBottomMargin(0.12);
  gPad->SetTicks(1,1);

  fr->SetMinimum(0);
  fr->SetMaximum(fr->GetMaximum()*1.15);
  fr->Draw();

  // CMS-style header
  TLatex head; head.SetNDC(); head.SetTextSize(0.035);
  head.DrawLatex(0.14,0.93,"CMS, #sqrt{s}=7 TeV, L = 217.8 pb^{-1}");

  // Legend (similar to the paper)
  TLegend leg(0.53,0.73,0.88,0.90);
  leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.032);
  leg.AddEntry(fr->findObject("dataPoints"), "Data", "lep");
  leg.AddEntry(fr->findObject("psCurve"),    "three-body PS", "l");
  leg.AddEntry(fr->findObject("totCurve"),   "2 BWs + three-body PS", "l");
  leg.Draw();

  // Caption-style label
  TLatex cap; cap.SetNDC(); cap.SetTextSize(0.030);
  cap.DrawLatex(0.38,0.02,"Two signal-hypothesis fit");

  // Save plots
  c.SaveAs("dm_twoRelBW_ps.png");
  c.SaveAs("dm_twoRelBW_ps.pdf");

  // (Optional) save rebinned Δm histogram used for the plot
  TFile fout("filtered_dm.root","RECREATE");
  hdata.Write();
  fout.Close();

  // Console summary
  printf("\n[Two signal-hypothesis: S-wave relBW ⊗ Gauss]\n");
  printf(" Peak1: m=%.4f±%.4f GeV | Γ0=%.4f±%.4f GeV | N1=%.1f±%.1f\n",
         m01.getVal(),m01.getError(), g01.getVal(),g01.getError(), N1.getVal(),N1.getError());
  printf(" Peak2: m=%.4f±%.4f GeV | Γ0=%.4f±%.4f GeV | N2=%.1f±%.1f\n",
         m02.getVal(),m02.getError(), g02.getVal(),g02.getError(), N2.getVal(),N2.getError());
  printf(" Resolution: σ=%.4f±%.4f GeV | NB=%.1f±%.1f\n",
         sigma.getVal(),sigma.getError(), NB.getVal(),NB.getError());

  f->Close();
}

