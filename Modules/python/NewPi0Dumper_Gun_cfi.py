import FWCore.ParameterSet.Config as cms

newPi0Dumper_Gun = cms.EDAnalyzer('NewPi0Dumper_Gun',
                              OutputFile = cms.untracked.string('NewPi0Tuple.root'),
                              pathExtGeo = cms.untracked.string('/afs/cern.ch/user/l/lpernie/scratch1/pi0Calib/pi0/CMSSW_4_2_4/src/CalibCode/submit/common/caloGeometry.root'),
                              MCTruthCollection = cms.untracked.InputTag("source"),
                              
                              EBRecHitCollectionTag = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEB"),
                              EERecHitCollectionTag = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEE"),
                              ESRecHitCollectionTag = cms.untracked.InputTag("ecalPreshowerRecHit","EcalRecHitsES"),
                              PFRecHitCollectionTag = cms.untracked.InputTag("particleFlowClusterECAL"),
                              conversionsTag = cms.untracked.InputTag("trackerOnlyConversions"),

                              useES            = cms.untracked.bool(False),
                              StoreMCTruth     = cms.untracked.bool(True),
                              StoreConversions = cms.untracked.bool(False),
                              StoreTransparencyCorrection = cms.untracked.bool(True),
                              useTracks = cms.untracked.bool(False),
                              useHCAL = cms.untracked.bool(False),
                              goodCollSelection = cms.untracked.bool(False),
                              useBeamSpotPosition = cms.untracked.bool(True),
                              
                              masspi0Cut = cms.untracked.double(1.),
                              ptpi0Cut = cms.untracked.double(0.7),
                              s1CluCut = cms.untracked.double(0.350),
                              s1CluCutEE = cms.untracked.double(0.5),
                              ptCluCut = cms.untracked.double(0.300),
                              s4s9CluCut = cms.untracked.double(0.6),
                              posCalcParameters = cms.PSet(
                                   T0_barl      = cms.double(7.4),
                                   T0_endc      = cms.double(3.1),
                                   T0_endcPresh = cms.double(1.2),
                                   LogWeighted  = cms.bool(True),
                                   W0           = cms.double(4.2),
                                   X0           = cms.double(0.89)
                                                          )
                              )
