# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

"""
JPsiTrkTrkTrk_cfi
-----------------
Default configuration of the EDAnalyzer "JPsiTrkTrkTrk" (LITE version with FULL SELECTION).

Inputs:
  - muonLabel: muon collection (their globalTrack() will be used).
  - trackLabel: track collection (e.g. generalTracks).
  - vertexLabel: primary vertex collection.
  - triggerResultsLabel: HLT results (usually from the "HLT" process).
  - hltPaths: list of HLT path names (they can be EXACT or include a single '*'
              wildcard, i.e. prefix*suffix; the .cc code already implements
              that matching).

Preselection (reminder):
  - RAW dimuon in [mJpsiMin, mJpsiMax], fitted with a constraint to m(J/ψ)
    (fallback to PSEUDO if the constraint fails).
  - Triplets of kaon candidates: pT and ΔR with respect to J/ψ_used, 2 OS pairs
    and |q1+q2+q3|=1 (if enabled), mKK_min < mKKminCut, and a B-mass window
    on the “used/fit” candidate.

Selection (enabled with applyFinalSelection):
  1) Kaons: pT>1.0 GeV, |η|≤2.4, ≥1 pixel hit and ≥3 strip hits.
  2) Muons: |η|≤2.4; |m(μμ)_raw − mPDG| ≤ 0.150 GeV; pT(J/ψ)>7 GeV if the event
     passed the HLT.
  3) φ: choose an OS kaon pair in [1.008, 1.035] GeV (keep the closest to m_φ).
  4) Kinematic refit of the B candidate with (μ,K) mass hypotheses.
  5) Vertex probabilities: P(J/ψ)>10% and P(B)>1% (constrained J/ψ if available).
  6) Significance: Lxy(J/ψ)/σ ≥ 3 (with respect to the PV, using covariances).

Compatibility:
  - CMSSW_5_3_X (OpenData 2011/2012).
"""

JPsiTrkTrkTrk = cms.EDAnalyzer(
    "JPsiTrkTrkTrk",

    # ----------------------
    # Input collections
    # ----------------------
    muonLabel           = cms.InputTag("muons"),
    trackLabel          = cms.InputTag("generalTracks"),
    vertexLabel         = cms.InputTag("offlinePrimaryVertices"),
    # If your HLT process has a different name (e.g. "HLT2"), change the third argument:
    triggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),

    # ----------------------
    # HLT paths (EXACT or with '*')
    # ----------------------
    hltPaths = cms.vstring(
        "HLT_DoubleMu3_Jpsi_v1",
        "HLT_DoubleMu3_Jpsi_v2",
        "HLT_Dimuon6p5_Jpsi_Displaced_v1",
        "HLT_Dimuon7_Jpsi_Displaced_v2",
        "HLT_Dimuon7_Jpsi_Displaced_v3"
        # Examples with wildcard (if needed):
        # "HLT_DoubleMu3_Jpsi_v*",
        # "HLT_Dimuon7_Jpsi_Displaced_v*",
    ),

    # ============================
    # UNTRACKED parameters (cuts)
    # ============================
    enforceHLT          = cms.untracked.bool(True),

    # ----- PRESELECTION -----
    # Muons (globalTrack → hits in the tracker component)
    minMuonPixelHits    = cms.untracked.int32(1),
    minMuonSiliconHits  = cms.untracked.int32(8),   # pixel+strip

    # J/ψ (RAW dimuon) mass window
    mJpsiMin            = cms.untracked.double(2.7),
    mJpsiMax            = cms.untracked.double(3.4),

    # Additional tracks (ΔR vs J/ψ_used)
    minTrkPt            = cms.untracked.double(0.5),
    maxDrJpsiTrk        = cms.untracked.double(1.5),

    # Charge combinatorics / OS pairs
    requireTwoOSPairs   = cms.untracked.bool(True),   # exactly 2 OS pairs
    requireQsumAbs1     = cms.untracked.bool(True),   # |q1+q2+q3| == 1

    # Mass cuts in preselection
    mKKminCut           = cms.untracked.double(1.06), # reject if mKK_min >= 1.06
    bMassMin            = cms.untracked.double(5.0),
    bMassMax            = cms.untracked.double(5.6),

    # Kinematic fit for J/ψ (constraint)
    muMassSigma         = cms.untracked.double(1e-6),

    # ----- FINAL SELECTION (after preselection) -----
    applyFinalSelection     = cms.untracked.bool(True),

    # Kaons
    selKaonPtMin            = cms.untracked.double(1.0),  # GeV
    selKaonAbsEtaMax        = cms.untracked.double(2.4),
    selKaonMinPixelHits     = cms.untracked.int32(1),
    selKaonMinStripHits     = cms.untracked.int32(3),

    # Muons
    selMuonAbsEtaMax        = cms.untracked.double(2.4),

    # J/ψ RAW mass ± window
    selJpsiRawMassWin       = cms.untracked.double(0.150),  # GeV

    # φ: mass window
    selPhiMassMin           = cms.untracked.double(1.008),
    selPhiMassMax           = cms.untracked.double(1.035),

    # pT(J/ψ) > 7 GeV if the event passed the HLT (no explicit run dependence)
    selRequireHLTForJpsiPt7 = cms.untracked.bool(True),
    selJpsiPtMinWithHLT     = cms.untracked.double(7.0),    # GeV

    # Vertex probabilities
    selJpsiVtxProbMin       = cms.untracked.double(0.10),   # 10%
    selBVtxProbMin          = cms.untracked.double(0.01),   # 1%
    selUseConstrVtxForJpsi  = cms.untracked.bool(True),     # use constrained P(vtx) if available

    # Lxy(J/ψ) significance
    selJpsiLxySigMin        = cms.untracked.double(3.0)
)

# -------------------------------------------------------------------------
# Reference clones (optional) — uncomment/tune in your cfg if needed.
# -------------------------------------------------------------------------

# 1) Without HLT requirement (tests without TriggerResults.accept)
# JPsiTrkTrkTrk_NoHLT = JPsiTrkTrkTrk.clone(
#     enforceHLT = cms.untracked.bool(False)
# )

# 2) Wider B-mass window (exploration/systematics)
# JPsiTrkTrkTrk_WideB = JPsiTrkTrkTrk.clone(
#     bMassMin = cms.untracked.double(4.8),
#     bMassMax = cms.untracked.double(5.8)
# )

# 3) Relaxed selection for quick testing (toggle categories as needed)
# JPsiTrkTrkTrk_Relaxed = JPsiTrkTrkTrk.clone(
#     applyFinalSelection     = cms.untracked.bool(True),
#     selKaonPtMin            = cms.untracked.double(0.8),
#     selKaonMinStripHits     = cms.untracked.int32(2),
#     selPhiMassMin           = cms.untracked.double(1.000),
#     selPhiMassMax           = cms.untracked.double(1.040),
#     selJpsiVtxProbMin       = cms.untracked.double(0.05),
#     selBVtxProbMin          = cms.untracked.double(0.005),
#     selJpsiLxySigMin        = cms.untracked.double(2.0)
# )
