// fit_Bplus_mass_dm_roofit.C
//
// Purpose:
//   Unbinned (NON-EXTENDED) RooFit of the B⁺ candidate mass (cand_mass_fit)
//   under a Δm = m(μμKK) − m(μμ) selection window, using the final
//   JPsiTrkTrkTrk merged ntuple.
//
// Model:
//   - Signal: single Gaussian (mean, sigma).
//   - Background: 2nd-order Chebyshev polynomial.
//   - Two fit scenarios:
//       (i) mean and sigma floating;
//       (ii) mean fixed to the PDG B⁺ mass, sigma floating.
//
// Inputs (defaults):
//   - infile   : "muonia_all_jpsitrktrktrk_merged.root"
//   - treePath : "myJPsiTrkTrkTrk/JPsiTrkTrkTrkTree"
//   - dm_low   : 1.008 GeV (lower edge of Δm window, m(μμKK) − m(μμ))
//   - dm_high  : 1.568 GeV (upper edge of Δm window)
//   - PDG_Bp   : 5.2793 GeV (PDG B⁺ mass used as reference)
//   - XMIN,XMAX: [5.15, 5.45] GeV (explicit mass window for cand_mass_fit)
//
// Output:
//   - PNG/PDF plots: fit_floatMean_floatSigma.{png,pdf},
//                    fit_fixedMean_floatSigma.{png,pdf}.
//   - Console summary of fSig, nSig(est.), mean and sigma for both scenarios.
//
// Requirements:
//   - ROOT >= 6 with RooFit enabled.

#include <vector>
#include <cmath>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLorentzVector.h>
#include <TString.h>

#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"

using namespace RooFit;

static const double M_K = 0.493677; // GeV

static TLorentzVector p4(float pt,float eta,float phi,double m){
  TLorentzVector v; const double px=pt*cos(phi), py=pt*sin(phi), pz=pt*sinh(eta);
  v.SetPxPyPzE(px,py,pz, std::sqrt(px*px+py*py+pz*pz+m*m)); return v;
}

// Select the K+K− pair corresponding to mKK_min
static void chooseKK(double m12,double m13,double m23,double mmin,int& i,int& j){
  const double tol=0.005; i=j=-1;
  if (fabs(m12-mmin)<tol){i=1;j=2;return;}
  if (fabs(m13-mmin)<tol){i=1;j=3;return;}
  if (fabs(m23-mmin)<tol){i=2;j=3;return;}
  const double d12=fabs(m12-mmin), d13=fabs(m13-mmin), d23=fabs(m23-mmin);
  if (d12<=d13 && d12<=d23){i=1;j=2;} else if (d13<=d12 && d13<=d23){i=1;j=3;} else {i=2;j=3;}
}

// Plot text block: labels Bp_peak/Bp_width/c1/c2/nBckPol/nSigBp
// - No box, bold font
// - Text anchored to the frame area (independent of pad margins)
static void stampNoExt(const RooRealVar& mean,const RooRealVar& sigma,
                       const RooRealVar& c1,const RooRealVar& c2,
                       const RooRealVar& fSig, int Ntot,
                       TPad* pad)
{
  const double nSig  = fSig.getVal()*Ntot;
  const double enSig = fSig.getError()*Ntot;
  const double nBkg  = (1.0 - fSig.getVal())*Ntot;
  const double enBkg = fSig.getError()*Ntot;

  // Pad/canvas margins
  const double xmin = pad->GetLeftMargin();
  const double ymin = pad->GetBottomMargin();
  const double xmax = 1.0 - pad->GetRightMargin();
  const double ymax = 1.0 - pad->GetTopMargin();

  // Relative position inside the usable frame area
  const double fx_center = 0.30;   // 0=left edge, 1=right edge (0.50 = centered)
  const double fy_top    = 0.30;   // distance from bottom edge

  // Convert to absolute NDC coordinates
  const double x = xmin + fx_center * (xmax - xmin);
  const double y = ymin + fy_top    * (ymax - ymin);
  const double dy = 0.04 * (ymax - ymin);  // line spacing, scaled with frame height

  TLatex t;
  t.SetNDC();
  t.SetTextFont(62);   // bold
  t.SetTextSize(0.032);
  t.SetTextColor(kBlack);
  t.SetTextAlign(23);  // horizontally centered at x, y anchored at top

  t.DrawLatex(x, y,            Form("Bp_peak  = %.4f #pm %.6f", mean.getVal(),  mean.getError()));
  t.DrawLatex(x, y - 1*dy,     Form("Bp_width = %.6f #pm %.6f", sigma.getVal(), sigma.getError()));
  t.DrawLatex(x, y - 2*dy,     Form("c1 = %.6f #pm %.6f", c1.getVal(), c1.getError()));
  t.DrawLatex(x, y - 3*dy,     Form("c2 = %.6f #pm %.6f", c2.getVal(), c2.getError()));
  t.DrawLatex(x, y - 4*dy,     Form("nBckPol = %.1f #pm %.2f", nBkg, enBkg));
  t.DrawLatex(x, y - 5*dy,     Form("nSigBp  = %.1f #pm %.2f", nSig, enSig));
}

// --------------------------------------------------------------------
// CMS-style header, anchored to the top margin of the canvas/pad
// --------------------------------------------------------------------
static void drawCmsHeader(TPad *pad,
                          const char *text = "CMS, #sqrt{s} = 7 TeV, L = 5.2 fb^{-1}")
{
  // Pad margins
  const double xmin = pad->GetLeftMargin();
  const double xmax = 1.0 - pad->GetRightMargin();
  const double top  = pad->GetTopMargin();

  // Horizontal center of the frame
  const double x = 0.5*(xmin + xmax);
  // Point inside the top margin (roughly mid-margin)
  const double y = 1.0 - 0.5*top;

  TLatex t;
  t.SetNDC();
  t.SetTextFont(62);   // bold
  t.SetTextSize(0.045);
  t.SetTextAlign(22);  // centered at (x,y)

  t.DrawLatex(x, y, text);
}

// Single NON-EXTENDED fit (fSig fraction), returns basic summary values
static void fitOnce(RooDataSet& data,
                    RooRealVar& mB,
                    const char* tag,
                    bool fixMean, double PDG,
                    int nbins,
                    // summarized outputs:
                    double& fSig_val, double& fSig_err,
                    double& mean_val,double& mean_err,
                    double& sigma_val,double& sigma_err)
{
  // Signal
  RooRealVar mean("mean","B^{+} mass", PDG, PDG-0.010, PDG+0.010);
  if (fixMean){ mean.setVal(PDG); mean.setConstant(kTRUE); }
  RooRealVar sigma("sigma","resolution", 0.010, 0.004, 0.020);
  RooGaussian sig("sig","signal", mB, mean, sigma);

  // Background: Chebyshev(2)
  RooRealVar c1("c1","c1", -0.25, -2.0, 2.0);
  RooRealVar c2("c2","c2",  0.00, -2.0, 2.0);
  RooChebychev bkg("bkg","bkg", mB, RooArgList(c1,c2));

  // NON-EXTENDED model (signal fraction)
  RooRealVar fSig("fSig","signal fraction", 0.2, 0.0, 1.0);
  RooAddPdf  model("model","sig+bkg", RooArgList(sig,bkg), RooArgList(fSig));

  // Fit
  model.fitTo(data, PrintLevel(-1));

  // Plot
  RooPlot* fr = mB.frame(Bins(nbins));
  data.plotOn(fr, DataError(RooAbsData::SumW2), Name("data_points"));
  model.plotOn(fr, Name("model_total"), LineColor(kBlue));                         // solid blue
  model.plotOn(fr, Components("bkg"), LineStyle(kDashed), LineColor(kBlue), Name("bkg_only")); // dashed blue

  // Canvas/style
  TCanvas c(Form("c_%s",tag), "", 900, 700);
  gStyle->SetOptStat(0);

  // Axes (X as m_{J/psi phiK}; Y keeps "Events / bin")
  fr->GetXaxis()->SetTitle("m_{J/#psi#phiK^{+}} (GeV)");
  fr->GetYaxis()->SetTitle("Events / bin");

  fr->GetXaxis()->SetTitleSize(0.045);
  fr->GetYaxis()->SetTitleSize(0.045);
  fr->GetXaxis()->SetLabelSize(0.04);
  fr->GetYaxis()->SetLabelSize(0.04);
  
  fr->SetTitle("");
  gStyle->SetOptTitle(0);

  // Margins
  c.SetLeftMargin(0.15);
  c.SetBottomMargin(0.15);
  c.SetRightMargin(0.10);
  c.SetTopMargin(0.10);

  fr->Draw();

  // CMS-style header at the top (you can update the text once you quote your effective L)
  drawCmsHeader(&c, "CMS, #sqrt{s} = 7 TeV, L = 5.2 fb^{-1}");

  // Bold text with fit parameters, without a box (using current canvas/pad)
  stampNoExt(mean, sigma, c1, c2, fSig, data.numEntries(), &c);

  c.SaveAs(Form("fit_%s.png", tag));
  c.SaveAs(Form("fit_%s.pdf", tag));

  // Outputs
  fSig_val = fSig.getVal(); fSig_err = fSig.getError();
  mean_val = mean.getVal(); mean_err = mean.getError();
  sigma_val= sigma.getVal(); sigma_err= sigma.getError();
}

void fit_Bplus_dm_roofit_textAligned(const char* infile="muonia_all_jpsitrktrktrk_merged.root",
                                     const char* treePath="myJPsiTrkTrkTrk/JPsiTrkTrkTrkTree",
                                     double dm_low=1.008,      // GeV (note)
                                     double dm_high=1.568,     // GeV (note)
                                     double PDG_Bp=5.2793,     // GeV
                                     double XMIN=5.15, double XMAX=5.45,
                                     int nbins=35)
{
  // Open input file and TTree
  TFile* f=TFile::Open(infile);
  if(!f||f->IsZombie()){printf("ERROR: cannot open %s\n", infile); return;}
  TTree* t=(TTree*)f->Get(treePath);
  if(!t){printf("ERROR: cannot find %s\n", treePath); return;}

  // Helper to set vector<float>* with alternative branch names (*_phl / *_phi)
  auto setVec = [&](const char* a, const char* b, std::vector<float>** p)->bool{
    if (t->GetBranch(a)) { t->SetBranchAddress(a,p); return true; }
    if (b && t->GetBranch(b)) { t->SetBranchAddress(b,p); return true; }
    printf("ERROR: missing branch '%s'%s%s\n", a, b?" or '":"", b?b:"");
    return false;
  };

  // Per-candidate vector branches
  std::vector<float> *dimu_pt_raw=0, *dimu_eta_raw=0, *dimu_phi_raw=0, *dimu_mass_raw=0;
  std::vector<float> *trk1_pt=0, *trk1_eta=0, *trk1_phi=0;
  std::vector<float> *trk2_pt=0, *trk2_eta=0, *trk2_phi=0;
  std::vector<float> *trk3_pt=0, *trk3_eta=0, *trk3_phi=0;
  std::vector<float> *m12=0, *m13=0, *m23=0, *mmin=0;
  std::vector<float> *mBfit=0;

  // Branch addresses
  t->SetBranchAddress("dimu_pt_raw",   &dimu_pt_raw);
  t->SetBranchAddress("dimu_eta_raw",  &dimu_eta_raw);
  if (!setVec("dimu_phl_raw","dimu_phi_raw",&dimu_phi_raw)) return;
  t->SetBranchAddress("dimu_mass_raw", &dimu_mass_raw);

  t->SetBranchAddress("trk1_pt",&trk1_pt); t->SetBranchAddress("trk1_eta",&trk1_eta); if(!setVec("trk1_phl","trk1_phi",&trk1_phi)) return;
  t->SetBranchAddress("trk2_pt",&trk2_pt); t->SetBranchAddress("trk2_eta",&trk2_eta); if(!setVec("trk2_phl","trk2_phi",&trk2_phi)) return;
  t->SetBranchAddress("trk3_pt",&trk3_pt); t->SetBranchAddress("trk3_eta",&trk3_eta); if(!setVec("trk3_phl","trk3_phi",&trk3_phi)) return;

  t->SetBranchAddress("mKK_12",&m12);
  t->SetBranchAddress("mKK_13",&m13);
  t->SetBranchAddress("mKK_23",&m23);
  t->SetBranchAddress("mKK_min",&mmin);
  t->SetBranchAddress("cand_mass_fit",&mBfit);

  // RooFit observable and unbinned dataset
  RooRealVar mB("cand_mass_fit","cand_mass_fit", XMIN, XMAX, "GeV");
  RooDataSet data("data","data", RooArgSet(mB));

  const Long64_t nEv = t->GetEntries();
  Long64_t nCandPass = 0, nCandTot = 0;

  for (Long64_t ie=0; ie<nEv; ++ie){
    t->GetEntry(ie);

    // Ensure consistent vector sizes for this event/candidate collection
    size_t n = mBfit->size();
    auto shrink = [&](const std::vector<float>* v){ n = std::min(n, v->size()); };
    shrink(dimu_pt_raw); shrink(dimu_eta_raw); shrink(dimu_phi_raw); shrink(dimu_mass_raw);
    shrink(trk1_pt); shrink(trk1_eta); shrink(trk1_phi);
    shrink(trk2_pt); shrink(trk2_eta); shrink(trk2_phi);
    shrink(trk3_pt); shrink(trk3_eta); shrink(trk3_phi);
    shrink(m12); shrink(m13); shrink(m23); shrink(mmin);

    nCandTot += n;

    for (size_t ic=0; ic<n; ++ic){
      // 4-vectors
      TLorentzVector pdimu = p4(dimu_pt_raw->at(ic), dimu_eta_raw->at(ic), dimu_phi_raw->at(ic), dimu_mass_raw->at(ic));
      TLorentzVector pk1   = p4(trk1_pt->at(ic), trk1_eta->at(ic), trk1_phi->at(ic), M_K);
      TLorentzVector pk2   = p4(trk2_pt->at(ic), trk2_eta->at(ic), trk2_phi->at(ic), M_K);
      TLorentzVector pk3   = p4(trk3_pt->at(ic), trk3_eta->at(ic), trk3_phi->at(ic), M_K);

      // Select the lowest-mass K+K− pair
      int a=-1,b=-1; chooseKK(m12->at(ic), m13->at(ic), m23->at(ic), mmin->at(ic), a, b);
      TLorentzVector pphi = (a==1&&b==2)? pk1+pk2 : (a==1&&b==3? pk1+pk3 : pk2+pk3);

      const double dm    = (pdimu+pphi).M() - pdimu.M();
      const double mBval = mBfit->at(ic);

      // Δm window + B-mass range
      if (dm>dm_low && dm<dm_high && mBval>XMIN && mBval<XMAX){
        mB.setVal(mBval);
        data.add(RooArgSet(mB));
        ++nCandPass;
      }
    }
  }

  printf("Eventos: %lld | Candidatos totales: %lld | Aceptados tras Δm y rango mB: %lld\n",
         nEv, nCandTot, nCandPass);
  if (nCandPass==0){ printf("Aviso: 0 candidatos tras cortes; revisa ramas y rangos.\n"); }

  // --- Fits (two scenarios) ---
  double f1=0,ef1=0, m1=0,em1=0, s1=0,es1=0;
  double f2=0,ef2=0, m2=0,em2=0, s2=0,es2=0;

  fitOnce(data, mB, "floatMean_floatSigma", false, PDG_Bp, nbins, f1,ef1, m1,em1, s1,es1);
  fitOnce(data, mB, "fixedMean_floatSigma", true,  PDG_Bp, nbins, f2,ef2, m2,em2, s2,es2);

  const int N = data.numEntries();
  printf("\n[Resumen NO-EXTENDIDO]\n");
  printf("  float mean & sigma : fSig = %.3f ± %.3f | nSig(est.) = %.1f | mean = %.6f ± %.6f | sigma = %.6f ± %.6f\n",
         f1, ef1, f1*N, m1, em1, s1, es1);
  printf("  fixed mean, float σ : fSig = %.3f ± %.3f | nSig(est.) = %.1f | mean = %.6f (fixed) | sigma = %.6f ± %.6f\n",
         f2, ef2, f2*N, PDG_Bp, s2, es2);

  f->Close();
}
