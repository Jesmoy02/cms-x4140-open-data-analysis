# -*- coding: utf-8 -*-
import os
import sys
import FWCore.ParameterSet.Config as cms
from FWCore.PythonUtilities import LumiList  # apply JSON of certified good runs and lumis

# =============================================================================
# 1) DYNAMIC PARAMETERS (CLI)
# -----------------------------------------------------------------------------
# Usage:
#   cmsRun poet_cfg.py input.root output.root
# If no output is given, JPsiTrkTrkTrkTree_output.root will be used.
# =============================================================================
if len(sys.argv) > 2:
    input_file = sys.argv[2]
else:
    raise RuntimeError("Debes proporcionar el archivo ROOT de entrada como tercer argumento.\n"
                       "Ejemplo: cmsRun poet_cfg.py input.root output.root")

if len(sys.argv) > 3:
    output_file = sys.argv[3]
else:
    output_file = "JPsiTrkTrkTrkTree_output.root"

print("[CFG] Archivo de entrada: %s" % input_file)
print("[CFG] Archivo de salida : %s" % output_file)

# Set to False when running on MC
isData = True  # MuOnia 2011/2012 → datos

# =============================================================================
# 2) PROCESS, SERVICES AND OPTIONS
# =============================================================================
process = cms.Process("JPsiTrkTrkTrkAna")

# MessageLogger (compact output with final summary)
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = "WARNING"
process.options = cms.untracked.PSet(
    wantSummary=cms.untracked.bool(True)
)

# Number of events (-1 = all)
process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(-1)
)

# =============================================================================
# 2.1) CONDITIONS / ES: GlobalTag, Geometry, B field and TTBuilder
# =============================================================================
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# ============================
# Helpers for ASCII-safe prints (Py2.6)
# ============================
def _plain(x):
    try:
        return x.value() if hasattr(x, 'value') else x
    except Exception:
        return x

def _fmt(x, nd=6):
    try:
        s = ("%." + str(nd) + "f") % float(x)
        s = s.rstrip('0').rstrip('.')
        return s
    except Exception:
        return str(x)

def _p(label, value, width=28):
    # Alignment using %, compatible with Py2.6
    print(("%-" + str(width) + "s %s") % (label, _plain(value)))

# ============================
# GlobalTag (DATA/MC)
# ============================
if isData:
    sqlite_path = '/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_data_stripped.db'
    if os.path.exists(sqlite_path):
        process.GlobalTag.connect   = cms.string('sqlite_file:' + sqlite_path)
        process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'
        _p("[GlobalTag] (DATA) source :", "SQLite OpenData")
        _p("[GlobalTag] (DATA) tag    :", process.GlobalTag.globaltag)
        _p("[GlobalTag] (DATA) connect:", process.GlobalTag.connect)
    else:
        process.GlobalTag.globaltag = 'GR_R_53_LV6::All'
        _p("[GlobalTag] (DATA) source :", "Frontier")
        _p("[GlobalTag] (DATA) tag    :", process.GlobalTag.globaltag)
else:
    sqlite_path = '/cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1_MC_stripped.db'
    if os.path.exists(sqlite_path):
        process.GlobalTag.connect   = cms.string('sqlite_file:' + sqlite_path)
        process.GlobalTag.globaltag = 'START53_LV6A1::All'
        _p("[GlobalTag] (MC) source   :", "SQLite OpenData")
        _p("[GlobalTag] (MC) tag      :", process.GlobalTag.globaltag)
        _p("[GlobalTag] (MC) connect  :", process.GlobalTag.connect)
    else:
        process.GlobalTag.globaltag = 'START53_LV6A1::All'
        _p("[GlobalTag] (MC) source   :", "Frontier")
        _p("[GlobalTag] (MC) tag      :", process.GlobalTag.globaltag)

# Ideal geometry (IdealGeometryRecord)
process.load('Configuration.Geometry.GeometryIdeal_cff')

# Static magnetic field 3.8T
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

# ESProducer for TransientTrackBuilder (TransientTrackRecord)
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')

# =============================================================================
# 2.2) PRE-HLT FILTER (Option A: enabled by default)
# -----------------------------------------------------------------------------
# - Filters by HLT BEFORE the analyzer.
# - Keep enforceHLT=True in the EDAnalyzer so that 'passedTriggers' is also stored.
# - andOr=True  ⇒ event passes if ANY of the paths is accepted.
# - throw=False ⇒ do not throw if a path does not exist in the menu.
# - EXACT PATH NAMES are used (no wildcards), in line with the reference note.
# =============================================================================
process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltFilterForJPsi = process.hltHighLevel.clone(
    HLTPaths = [
        "HLT_DoubleMu3_Jpsi_v1",
        "HLT_DoubleMu3_Jpsi_v2",
        "HLT_Dimuon6p5_Jpsi_Displaced_v1",
        "HLT_Dimuon7_Jpsi_Displaced_v2",
        "HLT_Dimuon7_Jpsi_Displaced_v3"
    ],
    andOr = True,
    throw = False
)

# =============================================================================
# 3) INPUT
# =============================================================================
# Supports local files and XRootD (root://)
if input_file.startswith("root://"):
    file_name = input_file
else:
    # Soft check of local file
    if not os.path.exists(input_file):
        raise RuntimeError("El archivo local no existe: %s" % input_file)
    file_name = "file:" + input_file

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(file_name)
)

# (Optional) Apply good-lumi JSON if found on disk
# ============================
# Good-lumi JSON (DATA)
# ============================
if isData:
    json_candidates = [
        "/code/CMSSW_5_3_32/src/PhysObjectExtractorTool/PhysObjectExtractor/python/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt",
        "data/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt",
    ]
    for _p_json in json_candidates:
        if os.path.exists(_p_json):
            print("[JSON] Usando: %s" % _p_json)
            myLumis = LumiList.LumiList(filename=_p_json).getCMSSWString().split(',')
            process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
            process.source.lumisToProcess.extend(myLumis)
            break

# =============================================================================
# 4) ROOT OUTPUT (TFileService)
# =============================================================================
process.TFileService = cms.Service(
    "TFileService",
    fileName=cms.string(output_file)
)

# =============================================================================
# 5) LOAD ANALYZER FROM THE _CFI
# -----------------------------------------------------------------------------
# Typical path for the _cfi within the package:
# PhysObjectExtractorTool/PhysObjectExtractor/python/JPsiTrkTrkTrk_cfi.py
# =============================================================================
process.load("PhysObjectExtractorTool.PhysObjectExtractor.JPsiTrkTrkTrk_cfi")

# =============================================================================
# 6) ANALYZER CONFIGURATION (clone from base _cfi)
# -----------------------------------------------------------------------------
# - Keeps:
#   * Exact HLT list from the _cfi (the analyzer also stores 'passedTriggers').
#   * enforceHLT = True.
#   * Default cuts consistent with .h/.cc implementation.
# - Below are commented examples for tests/debugging (without changing the base logic).
# =============================================================================
process.myJPsiTrkTrkTrk = process.JPsiTrkTrkTrk.clone(
    # --- (A) For debugging WITHOUT pre-HLT filter (use Path Option B below):
    # enforceHLT = cms.untracked.bool(False),

    # --- (B) Temporary “loose mode” to inspect triplets/candidates:
    # minMuonPixelHits    = cms.untracked.int32(0),
    # minMuonSiliconHits  = cms.untracked.int32(0),
    # minTrkPt            = cms.untracked.double(0.2),
    # maxDrJpsiTrk        = cms.untracked.double(3.0),
    # requireTwoOSPairs   = cms.untracked.bool(False),
    # requireQsumAbs1     = cms.untracked.bool(False),
    # mKKminCut           = cms.untracked.double(99.0),
    # bMassMin            = cms.untracked.double(4.0),
    # bMassMax            = cms.untracked.double(6.5),

    # --- (C) J/ψ mass window (parametrizable; default 2.7–3.4 GeV)
    # mJpsiMin = cms.untracked.double(2.7),
    # mJpsiMax = cms.untracked.double(3.4),

    # --- (D) Change the HLT list here if you add new specific versions:
    # hltPaths = cms.vstring(
    #     "HLT_DoubleMu3_Jpsi_v3",
    #     "HLT_Dimuon7_Jpsi_Displaced_v4",
    # ),

    # --- (E) Final selection enabled (1–6)
    # applyFinalSelection     = cms.untracked.bool(True),

    # # Kaons
    # selKaonPtMin            = cms.untracked.double(1.0),
    # selKaonAbsEtaMax        = cms.untracked.double(2.4),
    # selKaonMinPixelHits     = cms.untracked.int32(1),
    # selKaonMinStripHits     = cms.untracked.int32(3),

    # # Muons
    # selMuonAbsEtaMax        = cms.untracked.double(2.4),

    # # J/ψ RAW mass ± window
    # selJpsiRawMassWin       = cms.untracked.double(0.150),

    # # φ: mass window
    # selPhiMassMin           = cms.untracked.double(1.008),
    # selPhiMassMax           = cms.untracked.double(1.035),

    # # pT(J/ψ) > 7 GeV if it passed HLT
    # selRequireHLTForJpsiPt7 = cms.untracked.bool(True),
    # selJpsiPtMinWithHLT     = cms.untracked.double(7.0),

    # # Vertex probabilities
    # selJpsiVtxProbMin       = cms.untracked.double(0.10),
    # selBVtxProbMin          = cms.untracked.double(0.01),
    # selUseConstrVtxForJpsi  = cms.untracked.bool(True),

    # # Lxy(J/ψ) significance
    # selJpsiLxySigMin        = cms.untracked.double(3.0)
)

# =============================================================================
# PRINT BLOCK
# =============================================================================

def _fmt(x, nd=6):
    """Format floats without weird notation or trailing zeros (ASCII-safe)."""
    try:
        s = ("%." + str(nd) + "f") % float(x)
        s = s.rstrip('0').rstrip('.')
        return s
    except Exception:
        return str(x)

print("========== JPsiTrkTrkTrk CONFIG ==========")
print("[CFG] enforceHLT            = %s" % process.myJPsiTrkTrkTrk.enforceHLT.value())
print("[CFG] HLT paths (analyser)  = %s" % ", ".join(list(process.myJPsiTrkTrkTrk.hltPaths)))

print("--- PRESELECCIÓN ---")
print("[CFG] minMuonPixelHits      = %s" % process.myJPsiTrkTrkTrk.minMuonPixelHits.value())
print("[CFG] minMuonSiliconHits    = %s" % process.myJPsiTrkTrkTrk.minMuonSiliconHits.value())
print("[CFG] mJpsi window (GeV)    = [%s, %s]" % (
    _fmt(process.myJPsiTrkTrkTrk.mJpsiMin.value()),
    _fmt(process.myJPsiTrkTrkTrk.mJpsiMax.value())
))
print("[CFG] minTrkPt (GeV)        = %s" % _fmt(process.myJPsiTrkTrkTrk.minTrkPt.value()))
print("[CFG] maxDrJpsiTrk          = %s" % _fmt(process.myJPsiTrkTrkTrk.maxDrJpsiTrk.value()))
print("[CFG] requireTwoOSPairs     = %s" % process.myJPsiTrkTrkTrk.requireTwoOSPairs.value())
print("[CFG] requireQsumAbs1       = %s" % process.myJPsiTrkTrkTrk.requireQsumAbs1.value())
print("[CFG] mKKminCut (GeV)       = %s" % _fmt(process.myJPsiTrkTrkTrk.mKKminCut.value()))
print("[CFG] bMass window (GeV)    = [%s, %s]" % (
    _fmt(process.myJPsiTrkTrkTrk.bMassMin.value()),
    _fmt(process.myJPsiTrkTrkTrk.bMassMax.value())
))
print("[CFG] muMassSigma (GeV)     = %s" % _fmt(process.myJPsiTrkTrkTrk.muMassSigma.value()))

print("--- SELECCIÓN (1-6) ---")
print("[CFG] applyFinalSelection   = %s" % process.myJPsiTrkTrkTrk.applyFinalSelection.value())
print("[CFG] K pT min (GeV)        = %s" % _fmt(process.myJPsiTrkTrkTrk.selKaonPtMin.value()))
print("[CFG] K |eta| max           = %s" % _fmt(process.myJPsiTrkTrkTrk.selKaonAbsEtaMax.value()))
print("[CFG] K min pixel hits      = %s" % process.myJPsiTrkTrkTrk.selKaonMinPixelHits.value())
print("[CFG] K min strip hits      = %s" % process.myJPsiTrkTrkTrk.selKaonMinStripHits.value())
print("[CFG] mu |eta| max          = %s" % _fmt(process.myJPsiTrkTrkTrk.selMuonAbsEtaMax.value()))
print("[CFG] |m(mumu)_raw - mPDG| win  = %s" % _fmt(process.myJPsiTrkTrkTrk.selJpsiRawMassWin.value()))
print("[CFG] phi mass window (GeV) = [%s, %s]" % (
    _fmt(process.myJPsiTrkTrkTrk.selPhiMassMin.value()),
    _fmt(process.myJPsiTrkTrkTrk.selPhiMassMax.value())
))
print("[CFG] requireHLT for pT(J/psi) > 7 = %s" % process.myJPsiTrkTrkTrk.selRequireHLTForJpsiPt7.value())
print("[CFG] pT(J/psi) min with HLT (GeV) = %s" % _fmt(process.myJPsiTrkTrkTrk.selJpsiPtMinWithHLT.value()))
print("[CFG] VtxProb(J/psi) min    = %s" % _fmt(process.myJPsiTrkTrkTrk.selJpsiVtxProbMin.value()))
print("[CFG] VtxProb(B) min        = %s" % _fmt(process.myJPsiTrkTrkTrk.selBVtxProbMin.value()))
print("[CFG] use constrained Vtx(J/psi) = %s" % process.myJPsiTrkTrkTrk.selUseConstrVtxForJpsi.value())
print("[CFG] Lxy(J/psi)/sigma min  = %s" % _fmt(process.myJPsiTrkTrkTrk.selJpsiLxySigMin.value()))
print("==========================================")

# =============================================================================
# 7) PATH(s)
# -----------------------------------------------------------------------------
# Option A (ACTIVE): pre-HLT filter + analyzer (recommended for data).
# Option B (debug alternative): without pre-HLT filter, rely only on enforceHLT.
# =============================================================================
process.p = cms.Path(process.hltFilterForJPsi * process.myJPsiTrkTrkTrk)
# For debugging without pre-filter (and e.g. enforceHLT=False):
# process.p = cms.Path(process.myJPsiTrkTrkTrk)

# =============================================================================
# 8) QUICK NOTES
# -----------------------------------------------------------------------------
# - The analyzer uses a "constrained" J/ψ if the constraint converges; otherwise it
#   falls back to a PSEUDO candidate.
# - ΔR(J/ψ, trk) and the B mass use J/ψ_used (fit if available; otherwise PSEUDO).
# - The current code **avoids reusing** the muon tracks inside the triplet (internal veto).
# - To run on MC, set isData=False above and review the MC GlobalTag.
# =============================================================================
