import FWCore.ParameterSet.Config as cms

process = cms.Process("PI0DUMPER")

correctHits = True
useHLTFilter = True

if correctHits:
    print "CORRECTING HITS"
    
if useHLTFilter:
    print "FILTERING PI0 EVENTS"

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
process.GlobalTag.globaltag = 'GR_P_V17::All'
process.GlobalTag.toGet = cms.VPSet(
     cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRcd"),
                  tag = cms.string("EcalLaserAPDPNRatios_v6_noVPT_online"),
                  connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
             ),
     cms.PSet(record = cms.string("EcalLaserAlphasRcd"),
                  tag = cms.string("EcalLaserAlphas_v1_prompt"),
                  connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_31X_ECAL")
             )
     )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) ) 


if useHLTFilter:
     import copy
     from HLTrigger.HLTfilters.hltHighLevel_cfi import *
     process.AlcaP0Filter = copy.deepcopy(hltHighLevel)
     process.AlcaP0Filter.throw = cms.bool(False)
     process.AlcaP0Filter.HLTPaths = ["AlCa_EcalPi0_*"]


process.MessageLogger.cerr.FwkReport.reportEvery = 100

#process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMenuConfig_cff')
#process.load('L1TriggerConfig.L1GtConfigProducers.L1GtStableParametersConfig_cff')

import RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi

if correctHits:
    process.ecalPi0ReCorrected =  RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi.ecalRecHit.clone(
        doEnergyScale = cms.bool(False),
        doIntercalib = cms.bool(False),
        doLaserCorrections = cms.bool(True),
        EBRecHitCollection = cms.InputTag("hltAlCaPi0RecHitsFilter","pi0EcalRecHitsEB","HLT"),
        EERecHitCollection = cms.InputTag("hltAlCaPi0RecHitsFilter","pi0EcalRecHitsEE","HLT"),
        EBRecalibRecHitCollection = cms.string('pi0EcalRecHitsEB'),
        EERecalibRecHitCollection = cms.string('pi0EcalRecHitsEE')
        )


#NewPi0Dumper
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
    process.newPi0Dumper.EBRecHitCollectionTag = cms.untracked.InputTag("hltAlCaPi0RecHitsFilter","pi0EcalRecHitsEB","HLT")
    process.newPi0Dumper.EERecHitCollectionTag = cms.untracked.InputTag("hltAlCaPi0RecHitsFilter","pi0EcalRecHitsEE","HLT")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
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
if (process.newPi0Dumper.useBeamSpotPosition == True):
    process.p *= process.offlineBeamSpot

#build ntuple
process.p *= process.newPi0Dumper

# input goes here
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      #'/store/data/Run2010B/AlCaP0/RAW/v1/000/148/864/EC2B7078-2EE0-DF11-AF09-000423D33970.root'
      #'file:/cmsrm/pc24/mgrassi/80E17F48-5F71-E011-8CC0-001D09F25217.root'
      #'file:/cmsrm/pc24/mgrassi/8EE5BA1D-2471-E011-8BCA-003048F01E88.root'
      '/store/data/Run2011A/AlCaP0/RAW/v1/000/160/939/F8ED531B-1153-E011-9D2A-0030487CAEAC.root'
     )
)
