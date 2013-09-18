import FWCore.ParameterSet.Config as cms

process = cms.Process("PI0DUMPER")

correctHits = True

if correctHits == True:
    print "CORRECTING HITS"
    
from Geometry.CaloEventSetup.CaloTopology_cfi import *

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MixingNoPileUp_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
process.load("RecoTracker.GeometryESProducer.TrackerRecoGeometryESProducer_cfi")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'MC_3XY_V18::All'
# Higgs sample
process.GlobalTag.globaltag = 'GR_R_39X_V6::All'

if correctHits: 
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string("EcalIntercalibConstantsRcd"),
                 tag = cms.string("EcalIntercalibConstants_2010_V3_forAlCaRawTest"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_ECAL") #intercalib Ratios V03 Calib / HLT
                 )
        )
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string("EcalADCToGeVConstantRcd"),
                 tag = cms.string("EcalADCToGeVConstant_forAlCaRawTest"),
                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_ECAL") #ADCToGeV ratios V03 Calib/ HLT
                 )
        )
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRcd"),
                 tag = cms.string("EcalLaserAPDPNRatios_test_sat0"),
                 connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL") # laser correction V3
                 )
        )


# Minbias Data commissiong 10 jetmet 
#process.GlobalTag.globaltag = 'GR_R_36X_V12A::All'
process.MessageLogger.cerr.FwkReport.reportEvery = 100


#process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMenuConfig_cff')
#process.load('L1TriggerConfig.L1GtConfigProducers.L1GtStableParametersConfig_cff')

import RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi

if correctHits:
    process.ecalPi0ReCorrected =  RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi.ecalRecHit.clone(
        doEnergyScale = cms.bool(True),
        doIntercalib = cms.bool(True),
        doLaserCorrections = cms.bool(True),
        EBRecHitCollection = cms.InputTag("ecalPi0Corrected","pi0EcalRecHitsEB"),
        EERecHitCollection = cms.InputTag("ecalPi0Corrected","pi0EcalRecHitsEE"),
        EBRecalibRecHitCollection = cms.string('pi0EcalRecHitsEB'),
        EERecalibRecHitCollection = cms.string('pi0EcalRecHitsEE')
        )


process.load("Analysis.Modules.NewPi0Dumper_cfi")
process.newPi0Dumper.OutputFile = cms.untracked.string('NewPi0Tuple.root')
    
#Cuts for AlcaPi0
process.newPi0Dumper.ptpi0Cut = 1.2
process.newPi0Dumper.ptCluCut = 0.5
process.newPi0Dumper.s4s9CluCut = 0.8
process.newPi0Dumper.StoreTransparencyCorrection = cms.untracked.bool(True)
if correctHits:
    process.newPi0Dumper.EBRecHitCollectionTag = cms.untracked.InputTag("ecalPi0ReCorrected","pi0EcalRecHitsEB")
    process.newPi0Dumper.EERecHitCollectionTag = cms.untracked.InputTag("ecalPi0ReCorrected","pi0EcalRecHitsEE")
else:
    process.newPi0Dumper.EBRecHitCollectionTag = cms.untracked.InputTag("ecalPi0Corrected","pi0EcalRecHitsEB")
    process.newPi0Dumper.EERecHitCollectionTag = cms.untracked.InputTag("ecalPi0Corrected","pi0EcalRecHitsEE")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

process.p = cms.Path()

# if (process.newPi0Dumper.useBeamSpotPosition == True):
#     p *= process.offlineBeamSpot
    
if correctHits:
    print "ADDING RECALIB RECHIT MODULE WITH PARAMETERS"
    print "ENERGY SCALE "+str(process.ecalPi0ReCorrected.doEnergyScale)
    print "INTERCALIBRATION "+str(process.ecalPi0ReCorrected.doIntercalib)
    print "LASER "+str(process.ecalPi0ReCorrected.doLaserCorrections)
    process.p *= process.ecalPi0ReCorrected

process.p *= process.offlineBeamSpot*process.newPi0Dumper

# input goes here
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #'/store/mc/Spring10/Hgg_wbf_M120/GEN-SIM-RECO/START3X_V26_S09-v1/0028/A6883EFD-3948-DF11-AB50-0025B31E3D3C.root'
      #'/store/data/BeamCommissioning09/MinimumBias/RAW-RECO/May8thSkim_GOODCOLL-v1/0120/08E61D3F-725C-DF11-AEC9-00151796D7C0.root'
      #'/store/data/Commissioning10/MinimumBias/RECO/SD_JetMETTauMonitor-Jun14thSkim_v1/0119/1264B42A-547F-DF11-B8B5-002618943821.root'
      #'file:/cmsrm/pc23_2/meridian/AlcaP0.root'
#	'/store/data/Run2010B/AlCaP0/RAW/v1/000/148/864/EC2B7078-2EE0-DF11-AF09-000423D33970.root'
	'/store/data/Run2010B/AlCaP0/ALCARECO/ALCARECOEcalCalPi0Calib-v2/000/148/864/807DB5C4-FCE0-DF11-BBCC-001617C3B76A.root'
     )
)
