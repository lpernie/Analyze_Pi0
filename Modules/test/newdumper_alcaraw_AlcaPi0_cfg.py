import FWCore.ParameterSet.Config as cms

isGun = False

if isGun:
   correctHits = False
   useHLTFilter = False
   process = cms.Process('PI0DUMPERGUN')
else:
   correctHits = True
   useHLTFilter = True
   process = cms.Process('PI0DUMPER')

if correctHits:
    print 'CORRECTING HITS'
    
if useHLTFilter:
    print 'FILTERING PI0 EVENTS'

if useHLTFilter:
    import copy
    from HLTrigger.HLTfilters.hltHighLevel_cfi import *
    process.AlcaP0Filter = copy.deepcopy(hltHighLevel)
    process.AlcaP0Filter.throw = cms.bool(False)
    process.AlcaP0Filter.HLTPaths = ['AlCa_EcalPi0_*']

import RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi

from Geometry.CaloEventSetup.CaloTopology_cfi import *
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Geometry.TrackerGeometryBuilder.trackerGeometry_cfi')
process.load('RecoTracker.GeometryESProducer.TrackerRecoGeometryESProducer_cfi')
process.load('RecoVertex.BeamSpotProducer.BeamSpot_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#load transparency loss
process.GlobalTag.globaltag = 'GR_R_42_V21B::All'

process.options = cms.untracked.PSet( 
    wantSummary = cms.untracked.bool(False),
    SkipEvent   = cms.untracked.vstring('ProductNotFound')
) 

process.MessageLogger.cerr.FwkReport.reportEvery = 100

#process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMenuConfig_cff')
#process.load('L1TriggerConfig.L1GtConfigProducers.L1GtStableParametersConfig_cff')
if correctHits:
    process.ecalPi0ReCorrected =  RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi.ecalRecHit.clone(
        doEnergyScale = cms.bool(False),
        doIntercalib = cms.bool(False),
        doLaserCorrections = cms.bool(True),
        EBRecHitCollection = cms.InputTag('hltAlCaPi0RecHitsFilter','pi0EcalRecHitsEB', 'HLT'),
        EERecHitCollection = cms.InputTag('hltAlCaPi0RecHitsFilter','pi0EcalRecHitsEE', 'HLT'),
        EBRecalibRecHitCollection = cms.string('pi0EcalRecHitsEB'),
        EERecalibRecHitCollection = cms.string('pi0EcalRecHitsEE')
    )

#NewPi0Dumper_Gun
process.load('Analysis.Modules.NewPi0Dumper_Gun_cfi')
#process.newPi0Dumper_Gun.OutputFile = cms.untracked.string('/tmp/LocalPio0Gun.root')
if isGun:
   process.newPi0Dumper_Gun.OutputFile = cms.untracked.string('LocalPi0_Gun.root')
else:   
   process.newPi0Dumper_Gun.OutputFile = cms.untracked.string('LocalPi0_AlcaPi0.root')
process.newPi0Dumper_Gun.ExternalGeometry = cms.untracked.string('/afs/cern.ch/user/l/lpernie/scratch1/pi0Calib/pi0/CMSSW_4_2_4/src/CalibCode/submit/common/caloGeometry.root')
#Cuts for AlcaPi0
if isGun:
   process.newPi0Dumper_Gun.StoreMCTruth = cms.untracked.bool(True)
   process.newPi0Dumper_Gun.useES = cms.untracked.bool(False)
else:
   process.newPi0Dumper_Gun.StoreMCTruth = cms.untracked.bool(False)
   process.newPi0Dumper_Gun.useES = cms.untracked.bool(True)
process.newPi0Dumper_Gun.ptpi0Cut = 0.7
process.newPi0Dumper_Gun.s1CluCutEE = 0.5
process.newPi0Dumper_Gun.s1CluCut = 0.35
process.newPi0Dumper_Gun.ptCluCut = 0.35
process.newPi0Dumper_Gun.s4s9CluCut = 0.85
process.newPi0Dumper_Gun.DoOffGeom = cms.untracked.bool(False)
process.newPi0Dumper_Gun.StoreTransparencyCorrection = cms.untracked.bool(True)
### choosing proper input tag (recalibration module changes the collection names)
if( not(isGun) and correctHits) :
    process.newPi0Dumper_Gun.EBRecHitCollectionTag = cms.untracked.InputTag('ecalPi0ReCorrected','pi0EcalRecHitsEB')
    process.newPi0Dumper_Gun.EERecHitCollectionTag = cms.untracked.InputTag('ecalPi0ReCorrected','pi0EcalRecHitsEE')
    process.newPi0Dumper_Gun.ESRecHitCollectionTag = cms.untracked.InputTag('hltAlCaPi0RecHitsFilter','pi0EcalRecHitsES', 'HLT')
if isGun:
    process.newPi0Dumper_Gun.EBRecHitCollectionTag = cms.untracked.InputTag('ecalRecHit','EcalRecHitsEB', 'RECO')
    process.newPi0Dumper_Gun.EERecHitCollectionTag = cms.untracked.InputTag('ecalRecHit','EcalRecHitsEE', 'RECO')
    process.newPi0Dumper_Gun.ESRecHitCollectionTag = cms.untracked.InputTag('ecalPreshowerRecHit','EcalRecHitsES','RECO')
    #process.newPi0Dumper_Gun.endcapRawSuperClusterCollection = cms.untracked.InputTag('multi5x5SuperClusters','multi5x5EndcapSuperClusters')
process.newPi0Dumper_Gun.PFRecHitCollectionTag = cms.untracked.InputTag('particleFlowClusterECAL','','RECO')

if isGun:
    process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
    )
else:
    process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(700000)
    )

process.p = cms.Path()

#use HLT
if useHLTFilter:
    process.p *= process.AlcaP0Filter

if correctHits :
    print 'ADDING RECALIB RECHIT MODULE WITH PARAMETERS'
    print 'ENERGY SCALE '+str(process.ecalPi0ReCorrected.doEnergyScale)
    print 'INTERCALIBRATION '+str(process.ecalPi0ReCorrected.doIntercalib)
    print 'LASER '+str(process.ecalPi0ReCorrected.doLaserCorrections)
    process.p *= process.ecalPi0ReCorrected

#use beam spot
if (process.newPi0Dumper_Gun.useBeamSpotPosition == True):
    process.p *= process.offlineBeamSpot

#build ntuple
process.p *= process.newPi0Dumper_Gun

if isGun:
   process.source = cms.Source('PoolSource',
       fileNames = cms.untracked.vstring(
           #'file:/afs/cern.ch/user/l/lpernie/scratch1/pi0Calib/pi0/CMSSW_4_2_4/src/Pi0Gun/SinglePi0E10_1000Ev.root'
           'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_727_1_Wvj.root',
           'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_729_1_8CW.root',
           'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_72_1_ZbK.root',
           'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_730_1_SFX.root',
           'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_731_1_qc9.root',
           'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_732_1_t9s.root',
           'root://eoscms//eos/cms/store/caf/user/lpernie/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/Pi0_Gun_Pt1To15_CMSSW424_STARTV12_RECOSIM/0012ec1f5a1886e394599653b32b0ae8/SinglePi0E10_cfi_py_GEN_SIM_DIGI_L1_DIGI2RAW_RAW2DIGI_L1Reco_RECO_733_1_BgZ.root'
        )
   )
else:
   process.source = cms.Source('PoolSource',
       fileNames = cms.untracked.vstring(
#            'root://eoscms//eos/cms/store/data/Run2011B/AlCaP0/RAW/v1/000/177/789/DCFDCD2D-40EE-E011-A7DC-003048F118E0.root'
           'root://eoscms//eos/cms/store/data/Run2011B/AlCaP0/RAW/v1/000/175/906/D4FC42BE-9DDA-E011-AB80-001D09F232B9.root'
#           'root://eoscms//eos/cms/store/data/Run2011B/AlCaP0/RAW/v1/000/176/933/E8477639-4EE5-E011-884F-0030487C7828.root'
#           'root://eoscms//eos/cms/store/data/Run2011B/AlCaP0/RAW/v1/000/177/789/DCFDCD2D-40EE-E011-A7DC-003048F118E0.root'
        )
   )
