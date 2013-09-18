#ifndef ANALYSIS_MODULE_Pi0CalibTupleDumper_H
#define ANALYSIS_MODULE_Pi0CalibTupleDumper_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

// Geometry
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"

// ES stuff
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "RecoEcal/EgammaClusterAlgos/interface/PreshowerClusterAlgo.h"

// ROOT stuff
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"


typedef std::map<DetId, EcalRecHit> RecHitsMap;
// Less than operator for sorting EcalRecHits according to energy.
class ecalRecHitLess : public std::binary_function<EcalRecHit, EcalRecHit, bool> 
{
public:
  bool operator()(EcalRecHit x, EcalRecHit y) 
  { 
    return (x.energy() > y.energy()); 
  }
};


class TFile;
class TTree;


//class DQMStore;
//class MonitorElement;

class Pi0CalibTupleDumper : public edm::EDAnalyzer {

public:

  Pi0CalibTupleDumper( const edm::ParameterSet& );
  ~Pi0CalibTupleDumper();

protected:
   
  void beginJob();

  void beginRun(const edm::Run& r, const edm::EventSetup& c);

  void analyze(const edm::Event& e, const edm::EventSetup& c) ;

  void beginLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
                            const edm::EventSetup& context) ;

  void endLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
                          const edm::EventSetup& c);

  void endRun(const edm::Run& r, const edm::EventSetup& c);

  void endJob();

  void convxtalid(int & , int &);
  int diff_neta_s(int,int);
  int diff_nphi_s(int,int);



private:
 

  //DQMStore*   dbe_;  
  int eventCounter_;      
  PositionCalc posCalculator_ ;                        

  /// Pi0 invariant mass in EB
  ////MonitorElement * hMinvPi0EB_;

  /// Pi0 invariant mass in EE
  //MonitorElement * hMinvPi0EE_;

  /// Eta invariant mass in EB
  //MonitorElement * hMinvEtaEB_;

  /// Eta invariant mass in EE
  //MonitorElement * hMinvEtaEE_;

  /// Pt of the 1st most energetic Pi0 photon in EB 
  //MonitorElement *hPt1Pi0EB_;

  /// Pt of the 1st most energetic Pi0 photon in EE
  //MonitorElement *hPt1Pi0EE_;

  /// Pt of the 1st most energetic Eta photon in EB
  //MonitorElement *hPt1EtaEB_;

  /// Pt of the 1st most energetic Eta photon in EE
  //MonitorElement *hPt1EtaEE_;

  
  /// Pt of the 2nd most energetic Pi0 photon in EB
  //MonitorElement *hPt2Pi0EB_;

  /// Pt of the 2nd most energetic Pi0 photon in EE
  //MonitorElement *hPt2Pi0EE_;

  /// Pt of the 2nd most energetic Eta photon in EB
  //MonitorElement *hPt2EtaEB_;

  /// Pt of the 2nd most energetic Eta photon in EE
  //MonitorElement *hPt2EtaEE_;

  
  /// Pi0 Pt in EB
  //MonitorElement * hPtPi0EB_;

  /// Pi0 Pt in EE
  //MonitorElement * hPtPi0EE_;

  /// Eta Pt in EB
  //MonitorElement * hPtEtaEB_;

  /// Eta Pt in EE
  //MonitorElement * hPtEtaEE_;


  /// EB: eta Pi0, after selection
  //MonitorElement * eta_output_EBpi0_; 
  /// EB: phi Pi0, after selection
  //MonitorElement * phi_output_EBpi0_; 

  /// EE: eta Pi0, after selection
  //MonitorElement * eta_output_EEpi0_; 
  /// EE: phi Pi0, after selection
  //MonitorElement * phi_output_EEpi0_; 


  /// EB: eta Eta, after selection
  //MonitorElement * eta_output_EBeta_; 
  /// EB: phi Eta, after selection
  //MonitorElement * phi_output_EBeta_; 

  /// EE: eta Eta, after selection
  //MonitorElement * eta_output_EEeta_; 
  /// EE: phi Eta, after selection
  //MonitorElement * phi_output_EEeta_; 


  /// object to monitor
  edm::InputTag productMonitoredEBpi0_;
  edm::InputTag productMonitoredEBeta_;

 /// object to monitor
  edm::InputTag productMonitoredEEpi0_;
  edm::InputTag productMonitoredEEeta_;

 /// object to monitor
  edm::InputTag productMonitoredESpi0_;
  edm::InputTag productMonitoredESeta_;

      int gammaCandEtaSize_;
      int gammaCandPhiSize_;

      double seleXtalMinEnergy_;
      double seleXtalMinEnergyEndCap_;

  double clusSeedThr_;
  int clusEtaSize_;
  int clusPhiSize_;

  double clusSeedThrEndCap_;

      //// for pi0->gg barrel 
      double selePtGamma_;
      double selePtPi0_;
      double seleMinvMaxPi0_;
      double seleMinvMinPi0_;
      double seleS4S9Gamma_;
      double selePi0BeltDR_;
      double selePi0BeltDeta_;
      double selePi0Iso_;
      double ptMinForIsolation_; 

      ///for pi0->gg endcap
      double selePtGammaEndCap_;
      double selePtPi0EndCap_;
      double seleMinvMaxPi0EndCap_;
      double seleMinvMinPi0EndCap_;
      double seleS4S9GammaEndCap_;
      double selePi0IsoEndCap_;
      double selePi0BeltDREndCap_;
      double selePi0BeltDetaEndCap_;
      double ptMinForIsolationEndCap_; 

  double region1_Pi0EndCap_;
  double selePtGammaPi0EndCap_region1_; 
  double selePtPi0EndCap_region1_;
  double region2_Pi0EndCap_;
  double selePtGammaPi0EndCap_region2_; 
  double selePtPi0EndCap_region2_;
  double selePtGammaPi0EndCap_region3_; 
  double selePtPi0EndCap_region3_;



      ///for eta->gg barrel
      double selePtGammaEta_;
      double selePtEta_;
      double seleS4S9GammaEta_; 
      double seleS9S25GammaEta_; 
      double seleMinvMaxEta_; 
      double seleMinvMinEta_; 
      double ptMinForIsolationEta_; 
      double seleEtaIso_; 
      double seleEtaBeltDR_; 
      double seleEtaBeltDeta_; 

      ///for eta->gg endcap
      double selePtGammaEtaEndCap_;
      double seleS4S9GammaEtaEndCap_;
      double seleS9S25GammaEtaEndCap_;
      double selePtEtaEndCap_;
      double seleMinvMaxEtaEndCap_;
      double seleMinvMinEtaEndCap_;
      double ptMinForIsolationEtaEndCap_;
      double seleEtaIsoEndCap_;
      double seleEtaBeltDREndCap_;
      double seleEtaBeltDetaEndCap_;

  double region1_EtaEndCap_;
  double selePtGammaEtaEndCap_region1_; 
  double selePtEtaEndCap_region1_;
  double region2_EtaEndCap_;
  double selePtGammaEtaEndCap_region2_; 
  double selePtEtaEndCap_region2_;
  double selePtGammaEtaEndCap_region3_; 
  double selePtEtaEndCap_region3_;

  ///ES 

  int preshNclust_;
  float preshClustECut;
  double etThresh_;
  double calib_planeX_;
  double calib_planeY_;
  double mip_;
  double gamma_;
  PreshowerClusterAlgo * presh_algo; // algorithm doing the real work
  PreshowerClusterAlgo::DebugLevel debugL;  



  bool ParameterLogWeighted_;
  double ParameterX0_;
  double ParameterT0_barl_;
  double ParameterT0_endc_;
  double ParameterT0_endcPresh_;
  double ParameterW0_;



  std::vector<EBDetId> detIdEBRecHits; 
  std::vector<EcalRecHit> EBRecHits; 
 
  
  std::vector<EEDetId> detIdEERecHits; 
  std::vector<EcalRecHit> EERecHits; 



  /// Monitor every prescaleFactor_ events
  unsigned int prescaleFactor_;
  
  /// DQM folder name
  std::string folderName_; 
 
  /// Write to file 
  bool saveToFile_;

  /// which subdet will be monitored
  bool isMonEBpi0_;
  bool isMonEBeta_;
  bool isMonEEpi0_;
  bool isMonEEeta_;
  bool isMonESpi0_;
  bool isMonESeta_;

  /// Masks for EB EE(unused) and ES(unused) 
  bool isMaskEB_;
  bool isMaskEE_;
  bool isMaskES_;

  /// Use Reco Flag from RH
  bool useRecoFlag_;

  /// Files with xtals masked
  std::string maskEBFile_;
  std::string maskEEFile_;
  std::string maskESFile_;


  /// Output file name if required
  std::string fileName_;

  /// Store TTree 


  /// TTree for calibration Pi0 EB
  int             runn;
  int             eventn;
  float           iso;

  int nxt9[2];
  int nevpair;
  
  float xtalE9[2][9];
  int xtalEta9[2][9];
  int xtalPhi9[2][9];
  float xtalEtaf9[2][9];
  float xtalPhif9[2][9];
  
  /// TTree for calibration Eta EB
  int             runn_eta;
  int             eventn_eta;
  float           iso_eta;

  int nxt9_eta[2];
  int nevpair_eta;
  
  float xtalE9_eta[2][25];
  int xtalEta9_eta[2][25];
  int xtalPhi9_eta[2][25];
  float xtalEtaf9_eta[2][25];
  float xtalPhif9_eta[2][25];
  

  /// TTree for calibration Pi0 EE
  int             runn_ee;
  int             eventn_ee;
  float           iso_ee;

  int nxt9_ee[2];
  int nevpair_ee;
  
  float xtalE9_ee[2][9];
  int xtalEta9_ee[2][9];
  int xtalPhi9_ee[2][9];
  float xtalEtaf9_ee[2][9];
  float xtalPhif9_ee[2][9];
  
  /// TTree for calibration Eta EE
  int             runn_ee_eta;
  int             eventn_ee_eta;
  float           iso_ee_eta;

  int nxt9_ee_eta[2];
  int nevpair_ee_eta;
  
  float xtalE9_ee_eta[2][25];
  int xtalEta9_ee_eta[2][25];
  int xtalPhi9_ee_eta[2][25];
  float xtalEtaf9_ee_eta[2][25];
  float xtalPhif9_ee_eta[2][25];
  

  /// file and ttree objects : Pi0 EB
  TFile* m_file;
  TTree* m_tree;
  std::string m_outfilename;

  /// file and ttree objects : Pi0 EE
  TFile* m_file_ee;
  TTree* m_tree_ee;
  std::string m_outfilename_ee;

  /// file and ttree objects : Eta EB
  TFile* m_file_eta;
  TTree* m_tree_eta;
  std::string m_outfilename_eta;

  /// file and ttree objects : Eta EE
  TFile* m_file_ee_eta;
  TTree* m_tree_ee_eta;
  std::string m_outfilename_ee_eta;


};

#endif

