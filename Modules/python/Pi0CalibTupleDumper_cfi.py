import FWCore.ParameterSet.Config as cms

EcalPi0Mon = cms.EDFilter("Pi0CalibTupleDumper",
    prescaleFactor = cms.untracked.int32(1),
    FolderName = cms.untracked.string('AlCaReco/EcalPi0'),

    AlCaStreamEBpi0Tag = cms.untracked.InputTag("hltAlCaPi0RegRecHits","pi0EcalRecHitsEB"),
    AlCaStreamEEpi0Tag = cms.untracked.InputTag("hltAlCaPi0RegRecHits","pi0EcalRecHitsEE"),
    AlCaStreamESpi0Tag = cms.untracked.InputTag("hltAlCaPi0RegRecHits","pi0EcalRecHitsES"),

    AlCaStreamEBetaTag = cms.untracked.InputTag("hltAlCaEtaRegRecHits","etaEcalRecHitsEB"),
    AlCaStreamEEetaTag = cms.untracked.InputTag("hltAlCaEtaRegRecHits","etaEcalRecHitsEE"),
    AlCaStreamESetaTag = cms.untracked.InputTag("hltAlCaEtaRegRecHits","etaEcalRecHitsES"),

    isMonEEpi0 = cms.untracked.bool(True),
    isMonEBpi0 = cms.untracked.bool(True),
    isMonESpi0 = cms.untracked.bool(True),

    isMonEEeta = cms.untracked.bool(True),
    isMonEBeta = cms.untracked.bool(True),
    isMonESeta = cms.untracked.bool(True),

    isMaskEB = cms.untracked.bool(True),
    isMaskEE = cms.untracked.bool(False),
    isMaskES = cms.untracked.bool(False),

    maskEBFile = cms.untracked.string('maskEB.txt'),
    maskEEFile = cms.untracked.string('maskEE.txt'),
    maskESFile = cms.untracked.string('maskES.txt'),

    useRecoFlag = cms.untracked.bool(False),

    SaveToFile = cms.untracked.bool(False),
    FileName = cms.untracked.string('MonitorAlCaEcalPi0.root'),

    OutputFile = cms.untracked.string('AlCaPi0EBTree.root'),
    OutputFileEta = cms.untracked.string('AlCaEtaEBTree.root'),

    OutputFileEE = cms.untracked.string('AlCaPi0EETree.root'),
    OutputFileEEEta = cms.untracked.string('AlCaEtaEETree.root'),

    clusSeedThr = cms.double( 0.5 ),
    clusSeedThrEndCap = cms.double( 1.0 ),
    clusEtaSize = cms.int32( 3 ),
    clusPhiSize = cms.int32( 3 ),
    seleXtalMinEnergy = cms.double( -0.15 ),
    seleXtalMinEnergyEndCap = cms.double( -0.75 ),
    selePtGamma = cms.double(1 ),
    selePtPi0 = cms.double( 2. ),
    seleMinvMaxPi0 = cms.double( 0.22 ),
    seleMinvMinPi0 = cms.double( 0.06 ),
    seleS4S9Gamma = cms.double( 0.83 ),
    selePi0Iso = cms.double( 0.5 ),
    ptMinForIsolation = cms.double( 1 ),
    selePi0BeltDR = cms.double( 0.2 ),
    selePi0BeltDeta = cms.double( 0.05 ),

    selePtGammaEndCap = cms.double( 0.8 ),
    selePtPi0EndCap = cms.double( 3.0 ),
    seleS4S9GammaEndCap = cms.double( 0.9 ),
    seleMinvMaxPi0EndCap = cms.double( 0.3 ),
    seleMinvMinPi0EndCap = cms.double( 0.05 ),
    ptMinForIsolationEndCap = cms.double( 0.5 ),
    selePi0IsoEndCap = cms.double( 0.5 ),
    selePi0BeltDREndCap  = cms.double( 0.2 ),
    selePi0BeltDetaEndCap  = cms.double( 0.05 ),

region1_Pi0EndCap = cms.double(2),
selePtGammaPi0EndCap_region1 = cms.double(0.7),
selePtPi0EndCap_region1 = cms.double(3),

region2_Pi0EndCap = cms.double(2.5),
selePtGammaPi0EndCap_region2 = cms.double(0.5),
selePtPi0EndCap_region2 = cms.double(2),

selePtGammaPi0EndCap_region3 = cms.double(0.3),
selePtPi0EndCap_region3 = cms.double(1.2),



    selePtGammaEta = cms.double(1.2),
    selePtEta = cms.double(4.0),
    seleS4S9GammaEta  = cms.double(0.9),
    seleS9S25GammaEta  = cms.double(0.8),
    seleMinvMaxEta = cms.double(0.8),
    seleMinvMinEta = cms.double(0.3),
    ptMinForIsolationEta = cms.double(1.0),
    seleEtaIso = cms.double(0.5),
    seleEtaBeltDR = cms.double(0.3),
    seleEtaBeltDeta = cms.double(0.1),
    massLowPi0Cand = cms.double(0.104),
    massHighPi0Cand = cms.double(0.163),

    selePtGammaEtaEndCap = cms.double(1.5),
    selePtEtaEndCap = cms.double(5),
    seleS4S9GammaEtaEndCap  = cms.double(0.9),
    seleS9S25GammaEtaEndCap  = cms.double(0.85),
    seleMinvMaxEtaEndCap = cms.double(0.8),
    seleMinvMinEtaEndCap = cms.double(0.3),
    ptMinForIsolationEtaEndCap = cms.double(0.5),
    seleEtaIsoEndCap = cms.double(0.5),
    seleEtaBeltDREndCap = cms.double(0.3),
    seleEtaBeltDetaEndCap = cms.double(0.1),

region1_EtaEndCap = cms.double(2),
selePtGammaEtaEndCap_region1 = cms.double(1.3),
selePtEtaEndCap_region1 = cms.double(5),

region2_EtaEndCap = cms.double(2.5),
selePtGammaEtaEndCap_region2 = cms.double(1),
selePtEtaEndCap_region2 = cms.double(4),

selePtGammaEtaEndCap_region3 = cms.double(0.7),
selePtEtaEndCap_region3 = cms.double(3),

    preshStripEnergyCut      = cms.double(0.0),
    preshClusterEnergyCut    = cms.double(0.0),
    preshSeededNstrip        = cms.int32(15),
    preshNclust              = cms.int32(4),
    preshCalibPlaneX  = cms.double(1.0),
    preshCalibPlaneY  = cms.double(0.70),
    preshCalibGamma   = cms.double(0.024),
    preshCalibMIP     = cms.double(9.0E-5),
    debugLevelES      = cms.string(" "),
                              
    ParameterLogWeighted = cms.bool( True ),
    ParameterX0 = cms.double( 0.89 ),
    ParameterT0_barl = cms.double( 5.7 ),
    ParameterT0_endc = cms.double( 3.1 ),
    ParameterT0_endcPresh = cms.double( 1.2 ),
    ParameterW0 = cms.double( 4.2 )

)

                          
