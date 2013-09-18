import FWCore.ParameterSet.Config as cms

correctHits = False
useHLTFilter = False
process = cms.Process("PI0DUMPERGUN")

if correctHits:
    print "CORRECTING HITS"
    
if useHLTFilter:
    print "FILTERING PI0 EVENTS"

import RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi

from Geometry.CaloEventSetup.CaloTopology_cfi import *
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
process.load("RecoTracker.GeometryESProducer.TrackerRecoGeometryESProducer_cfi")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#load transparency loss
process.GlobalTag.globaltag = 'GR_R_42_V21B::All'

process.options = cms.untracked.PSet( 
    wantSummary = cms.untracked.bool(False),
    SkipEvent   = cms.untracked.vstring('ProductNotFound')
) 


if useHLTFilter:
     import copy
     from HLTrigger.HLTfilters.hltHighLevel_cfi import *
     process.AlcaP0Filter = copy.deepcopy(hltHighLevel)
     process.AlcaP0Filter.throw = cms.bool(False)
     process.AlcaP0Filter.HLTPaths = ["AlCa_EcalPi0_*"]


process.MessageLogger.cerr.FwkReport.reportEvery = 100

#process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMenuConfig_cff')
#process.load('L1TriggerConfig.L1GtConfigProducers.L1GtStableParametersConfig_cff')

if correctHits:
    process.ecalPi0ReCorrected =  RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi.ecalRecHit.clone(
        doEnergyScale = cms.bool(True),
        doIntercalib = cms.bool(True),
        doLaserCorrections = cms.bool(True),
        EBRecHitCollection = cms.InputTag('ecalRecHit','EcalRecHitsEB', 'RECO'),
        EERecHitCollection = cms.InputTag('ecalRecHit','EcalRecHitsEE', 'RECO'),
        EBRecalibRecHitCollection = cms.string('EcalRecHitsEB'),
        EERecalibRecHitCollection = cms.string('EcalRecHitsEE')
    )
#NewPi0Dumper_Gun
process.load("Analysis.Modules.NewPi0Dumper_Gun_cfi")
##process.newPi0Dumper_Gun.OutputFile = cms.untracked.string('/tmp/LocalPio0Gun.root')
process.newPi0Dumper_Gun.OutputFile = cms.untracked.string('LocalPio0Gun_GEOLOCAL.root')
#process.newPi0Dumper_Gun.pathExtGeo = cms.untracked.string('/afs/cern.ch/user/l/lpernie/scratch1/pi0Calib/pi0/CMSSW_4_2_4/src/CalibCode/submit/common/caloGeometry.root')
process.newPi0Dumper_Gun.ExternalGeometry = cms.untracked.string('/afs/cern.ch/user/l/lpernie/scratch1/pi0Calib/pi0/CMSSW_4_2_4/src/CalibCode/submit/common/caloGeometry.root')
#Cuts for AlcaPi0
process.newPi0Dumper_Gun.ptpi0Cut = 0.7
process.newPi0Dumper_Gun.s1CluCutEE = 0.5
process.newPi0Dumper_Gun.s1CluCut = 0.35
process.newPi0Dumper_Gun.ptCluCut = 0.35
process.newPi0Dumper_Gun.s4s9CluCut = 0.85
process.newPi0Dumper_Gun.DoOffGeom = cms.untracked.bool(False)
process.newPi0Dumper_Gun.StoreTransparencyCorrection = cms.untracked.bool(True)
### choosing proper input tag (recalibration module changes the collection names)
if correctHits:
    process.newPi0Dumper_Gun.EBRecHitCollectionTag = cms.untracked.InputTag("ecalPi0ReCorrected","pi0EcalRecHitsEB")
    process.newPi0Dumper_Gun.EERecHitCollectionTag = cms.untracked.InputTag("ecalPi0ReCorrected","pi0EcalRecHitsEE")
else:
    process.newPi0Dumper_Gun.EBRecHitCollectionTag = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEB", "RECO")
    process.newPi0Dumper_Gun.EERecHitCollectionTag = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEE", "RECO")
process.newPi0Dumper_Gun.ESRecHitCollectionTag = cms.untracked.InputTag("ecalPreshowerRecHit","EcalRecHitsES","RECO")
process.newPi0Dumper_Gun.PFRecHitCollectionTag = cms.untracked.InputTag("particleFlowClusterECAL","","RECO")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.p = cms.Path()

#use HLT
if useHLTFilter:
    process.p *= process.AlcaP0Filter

if correctHits:
    print "ADDING RECALIB RECHIT MODULE WITH PARAMETERS"
    print "ENERGY SCALE "+str(process.ecalPi0ReCorrected.doEnergyScale)
    print "INTERCALIBRATION "+str(process.ecalPi0ReCorrected.doIntercalib)
    print "LASER "+str(process.ecalPi0ReCorrected.doLaserCorrections)
    process.p *= process.ecalPi0ReCorrected

#use beam spot
if (process.newPi0Dumper_Gun.useBeamSpotPosition == True):
    process.p *= process.offlineBeamSpot

#build ntuple
process.p *= process.newPi0Dumper_Gun

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:/afs/cern.ch/user/l/lpernie/scratch1/pi0Calib/pi0/CMSSW_4_2_4/src/Pi0Gun/SinglePi0E10_1000Ev.root'
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_93_1_C40.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_940_1_DCz.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_941_1_o4O.root',            
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_942_1_iXr.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_943_1_WMC.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_944_1_OIt.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_945_1_8Xq.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_885_1_F2N.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_983_1_MCD.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_946_1_Zx4.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_92_1_TAS.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_999_1_SAu.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_989_1_mbP.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_974_1_bfb.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_957_1_wKy.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_960_1_hkw.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_952_1_OFG.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_998_1_MCI.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_95_1_nWm.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_949_1_t0e.root',
        'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_754_1_hBT.root'
     )
)
