#ifndef EcalCalibAlgo_H
#define EcalCalibAlgo_H

//
// $Id: EcalCalibAlgo.h,v 1.1 2013/02/04 13:23:02 lpernie Exp $
// Shahram Rahatlou, Sapienza Universita` di Roma & INFN
// May 2010

#include <vector>
#include <iostream>
#include <math.h>
#include <stdexcept>
#include <string>
#include <iomanip>

using namespace std;

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TProfile2D.h"

#include "PhysicsTools/FWLite/interface/EventContainer.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

#include "RecoCaloTools/Navigation/interface/EcalPreshowerNavigator.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "Geometry/CaloTopology/interface/EcalBarrelHardcodedTopology.h"
#include "Geometry/CaloTopology/interface/EcalEndcapHardcodedTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerHardcodedTopology.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerHardcodedGeometry.h"
#include "Geometry/EcalAlgo/interface/ECALGeometry.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "Analysis/Pi0Calib/interface/EcalCalibMap.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "Analysis/Pi0Calib/interface/PosCalcParams.h"

#define PI0MASS 0.135

//this class simply stores all the basic informations of a preshowercluster
//is used as output of the clustering algo 


class PreshowerCluster;

class RegionWeight {
 public:
  uint32_t   iRegion;
  float value;
}; 
typedef std::vector<RegionWeight> RegionWeightVector;

typedef std::pair<DetId, float> EnergyFraction;
typedef std::vector< EnergyFraction > EnergyFractionVector;
//used in the preshower clustering algo
typedef std::map<DetId, EcalRecHit> RecHitsMap;
typedef std::set<DetId> HitsID;

//----------------------------------------
template<class Type,int NMaxIter=10> class EcalCalibAlgo {
  //----------------------------------------
 public:
 //EcalCalibAlgo(fwlite::EventContainer* evt);
 EcalCalibAlgo(optutl::CommandLineParser* p);
 
 ~EcalCalibAlgo() {
   delete ebtopology;
   delete eetopology;
   delete estopology;
   delete geom;
   delete hardcodedPreshowerGeometry;
   //delete[] corrFunc; // causes crash at the end. why??
 }

 void iterate();
 void analyze(int it=0);
 float DeltaPhi(float phi1, float phi2);
 float GetDeltaR(float eta1, float eta2, float phi1, float phi2);
 void findESRoad(int stripwindow, ESDetId strip, EcalPreshowerNavigator theESNav, int plane);
 bool goodStrip(RecHitsMap::iterator candidate_it);
 void loadCalibMapFromFile(const char* cfile);
 void loadCorrectionsFromFile(const char* cfile);
 static const int MaxEnergyBins = 50;
 void fillResultTree(); // store coeff and DetId info for each xtal
 PreshowerCluster makeOnePreshowerCluster(int stripwindow,ESDetId *strip,HitsID *used_strips, RecHitsMap *the_rechitsMap_p,CaloSubdetectorTopology* topology_p);

 private:
 
 // compute regional weights for a given cluster
 RegionWeightVector getWeights(const reco::CaloCluster* clus) const;
 // { return RegionWeightVector(); } 
 
 void printIterSummary(int iter);// print eps,calib per region for iteration iter
 void resetParams(); // reset parameters
 void bookHistograms(); // book histograms per iteration - probably just the last one
 void makeSummaryHisto(); // book summary histograms
 
 //making vector of strips of the road reconstructed
 vector<ESDetId> esroad_2d;

 // The set of used DetID's for a given event used in the preshower clustering algo:
 HitsID  used_strips;

// creates the map of rechits used in the preshower clustering algo
 RecHitsMap  rechits_map;

 optutl::CommandLineParser* parser;
 fwlite::EventContainer* eventCont; // event container to access collections and 
 // outut histograms
 CaloTopology *ebtopology;         // hardcoded topology
 CaloTopology *eetopology;
 CaloSubdetectorTopology *estopology;
 ECALGeometry* geom;             // hardcoded geometry from root file
 EcalPreshowerHardcodedGeometry* hardcodedPreshowerGeometry;
 EcalCalibMap<Type>*   calibMap; // calibration coeffcients
 
 // these are needed for position calculation
 std::map<std::string,double> params_;
 PosCalcParams          PCparams_;
 
 // energy corrections
 TF1* corrFunc[MaxEnergyBins];
 bool noCorrections_;            // cache to decide whether we have corrections to apply
 float lowEnergyEdge[MaxEnergyBins+1];
 void fillLowEnergyEdge();
 float energyCorrection(float energy, float eta);
 
 
 float eps[Type::nRegions];   // eps_j for each region
 float calibCoeff[Type::nRegions][NMaxIter];   // 1/[1+eps_j] for each region
 float tmpcalib[Type::nRegions];   // 1/[1+eps_j] for each region
 float entries[Type::nRegions][NMaxIter];   // weighted number of entries
 
 float Numer[Type::nRegions]; //           sum_k <eps>^k W_j^k
                              // eps_j =  --------------------
 float Denom[Type::nRegions]; //               sum_k W_j^k
 
 bool lastIter_; // cache variable to book&fill histograms ONLY in the last iteration... for now

 //PRESHOWER CONSTANTS
 
 double mip_;
 double gamma_;
 double calib_planeX_;
 double calib_planeY_;


 //DEBUG COUNTERS & VAR

 int clusterwindowsize;

};


//========================================================================================
template<class Type,int NMaxIter>
  //EcalCalibAlgo<Type>::EcalCalibAlgo(fwlite::EventContainer* evt) {
  EcalCalibAlgo<Type,NMaxIter>::EcalCalibAlgo(optutl::CommandLineParser* p) {
  //========================================================================================
  //if(!evt) throw std::invalid_argument("null pointer for fwlite::EventContainer. No data. no party!");
  cout << "building EcalCalibAlgo ..." ;
  parser = p;
  //eventCont = evt;

   //DEBUG COUNTERS & VAR
   
   clusterwindowsize = 15;


  //PRESHOWER CONSTANTS

    mip_ = 81.1e-06;
    gamma_ = 0.024;
    calib_planeX_ = 1.0;
    calib_planeY_ = 0.7;

  //Topology and Geometry initialization

  ebtopology = new CaloTopology();  
  EcalBarrelHardcodedTopology* ebTopology=new EcalBarrelHardcodedTopology();
  ebtopology->setSubdetTopology(DetId::Ecal,EcalBarrel,ebTopology);
  
  eetopology = new CaloTopology();  
  EcalEndcapHardcodedTopology* eeTopology=new EcalEndcapHardcodedTopology();
  eetopology->setSubdetTopology(DetId::Ecal,EcalEndcap,eeTopology);
 
  estopology = new EcalPreshowerHardcodedTopology();
 
  calibMap = new EcalCalibMap<Type>();
 
  //Initializing properly the ECAL Geometry 
  TFile* f = TFile::Open("caloGeometry.root");
  geom = ECALGeometry::getGeometry(f);
  
  hardcodedPreshowerGeometry = new EcalPreshowerHardcodedGeometry(geom);
  hardcodedPreshowerGeometry->initializeParms();


  //std::map<std::string,double> providedParameters;
  //providedParameters.insert(std::make_pair("LogWeighted",true));
  //providedParameters.insert(std::make_pair("T0_barl",5.7));
  //providedParameters.insert(std::make_pair("T0_endc",3.1));
  //providedParameters.insert(std::make_pair("T0_endcPresh",1.2));
  //providedParameters.insert(std::make_pair("W0",4.2));
  //providedParameters.insert(std::make_pair("X0",0.89));
  PCparams_.param_LogWeighted_ = true;
  PCparams_.param_T0_barl_     = 5.7;
  PCparams_.param_T0_endc_     = 3.1;
  PCparams_.param_T0_endcES_   = 1.2;
  PCparams_.param_W0_          = 4.2;
  PCparams_.param_X0_          = 0.89;
  
  // initialize all parameter vectors
  resetParams();
  // the following should be reset ONLY at the very beginning
  for(uint32_t i=0; i<Type::nRegions;++i) {
    tmpcalib[i] = 1.;
    for(uint32_t j=0; j<NMaxIter; ++j) {
      entries[i][j] = 0.;
      calibCoeff[i][j] = 0.;
    }
  }
  
  // No correction fucntions by default
  for(int i=0; i<MaxEnergyBins; ++i) {
    corrFunc[i] = 0;
  }
  fillLowEnergyEdge();
  noCorrections_ = true;
  
  lastIter_ = false;
  
  cout << "... done. " << endl;
}


//============================================
template<class Type,int NMaxIter>
void EcalCalibAlgo<Type,NMaxIter>::printIterSummary(int iIter) {
//============================================
   for(uint32_t j=0; j<Type::nRegions; ++j)  {
     cout << "Region: " << Type::detIdFromRegion(j)
          << "  #cands: " << entries[j][iIter]
          << "\teps: " << eps[j]
          << "  1/(1+eps): " << 1./(1.+eps[j])
          << "  calib coeff: " << tmpcalib[j]
          << endl; 
   }
}

//============================================
template<class Type,int NMaxIter>
void EcalCalibAlgo<Type,NMaxIter>::resetParams() {
//============================================
  
  
  for(uint32_t i = 0; i<Type::nRegions; ++i) {
    eps[i] = 0.;
    Numer[i] = 0.;
    Denom[i] = 0.;
  }
}


// eta on x and phi on y
TH2F* makeEBXtalTH2(const char* name, const char* title) {
  TH2F* h2 = new TH2F(name,title,2*EBDetId::MAX_IETA+1,-EBDetId::MAX_IETA-0.5,EBDetId::MAX_IETA+0.5
		      ,EBDetId::MAX_IPHI, EBDetId::MIN_IPHI+0.5, EBDetId::MAX_IPHI+0.5 );
  h2->GetXaxis()->SetTitle("i#eta");
  h2->GetYaxis()->SetTitle("i#phi");
  return h2;
}


// eta on y and phi on x
TH2F* makeEBXtalTH2Std(const char* name, const char* title) {
  TH2F* h2 = new TH2F(name,title,
		      EBDetId::MAX_IPHI, EBDetId::MIN_IPHI+0.5, EBDetId::MAX_IPHI+0.5, 
		      2*EBDetId::MAX_IETA+1,-EBDetId::MAX_IETA-0.5,EBDetId::MAX_IETA+0.5
		      );
  h2->GetYaxis()->SetTitle("i#eta");
  h2->GetXaxis()->SetTitle("i#phi");
  return h2;
}

//============================================
template<class Type,int NMaxIter>
void EcalCalibAlgo<Type,NMaxIter>::bookHistograms() {
//============================================
  if(!lastIter_) return;
  
  // pi0 mass with no cut
  eventCont->add( 
		 new TH2F("Pi0MassvsetaEB","#gamma #gamma invariant mass vs. #eta",170,-1.479,1.479,120,0.050, 0.200)
		  );
  eventCont->hist("Pi0MassvsetaEB")->GetXaxis()->SetTitle("#eta");
  eventCont->hist("Pi0MassvsetaEB")->GetYaxis()->SetTitle("#gamma#gamma invariant mass");
  eventCont->hist("Pi0MassvsetaEB")->SetStats(false);
  
  int histo_bins;
  double histo_min_range;
  double histo_max_range;
  double bin_size;
  
  histo_bins = 120;
  histo_min_range = 0.050 ;
  histo_max_range = 0.200;
  bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins)*1000; //Mev
  
  
  ostringstream strs;
  strs<<setprecision(3)<<bin_size;
  string s_amplitude = strs.str();
  strs.clear();
  strs.str("");
  
  eventCont->add( new TH1F("Pi0MassEB","#gamma #gamma invariant mass for |#eta|<1.5",histo_bins ,histo_min_range,histo_max_range) );
  eventCont->hist("Pi0MassEB")->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");
  eventCont->hist("Pi0MassEB")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c^{2}").c_str());
  eventCont->hist("Pi0MassEB")->SetFillColor(kOrange);
  
  histo_bins = 120;
  histo_min_range = 0.050 ;
  histo_max_range = 0.200;
  bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins)*1000; //Mev
  
  strs<<setprecision(3)<<bin_size;
  s_amplitude = strs.str();
  strs.clear();
  strs.str("");
  
  eventCont->add( new TH1F("Pi0MassEE","#gamma #gamma invariant mass for 1.566<|#eta|<2.5",histo_bins ,histo_min_range, histo_max_range) );
  eventCont->hist("Pi0MassEE")->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");
  eventCont->hist("Pi0MassEE")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c^{2}").c_str());
  eventCont->hist("Pi0MassEE")->SetFillColor(kOrange);
  
  eventCont->add( new TH1F("Pi0MassEEetacut","#gamma #gamma invariant mass for 1.479<|#eta|<2.2",histo_bins ,histo_min_range, histo_max_range) );
  eventCont->hist("Pi0MassEEetacut")->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");
  eventCont->hist("Pi0MassEEetacut")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c^{2}").c_str());
  eventCont->hist("Pi0MassEEetacut")->SetFillColor(kOrange);
 
  eventCont->add( new TH1F("Pi0MassESEE","#gamma #gamma invariant mass corrected with Preshower energy for 1.479<|#eta|<2.5",histo_bins ,histo_min_range, histo_max_range) );
  eventCont->hist("Pi0MassESEE")->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");
  eventCont->hist("Pi0MassESEE")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c^{2}").c_str());
  eventCont->hist("Pi0MassESEE")->SetFillColor(kOrange);

  eventCont->add( new TH1F("Pi0MassESEEetacut","#gamma #gamma invariant mass corrected with Preshower energy for 1.479<|#eta|<1.22",histo_bins ,histo_min_range, histo_max_range) );
  eventCont->hist("Pi0MassESEEetacut")->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");
  eventCont->hist("Pi0MassESEEetacut")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c^{2}").c_str());
  eventCont->hist("Pi0MassESEEetacut")->SetFillColor(kOrange);

  eventCont->add( new TH1F("Pi0MassESEEetacutr1","#gamma #gamma invariant mass corrected with Preshower energy for 1.479<|#eta|<1.7",histo_bins ,histo_min_range, histo_max_range) );
  eventCont->hist("Pi0MassESEEetacutr1")->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");
  eventCont->hist("Pi0MassESEEetacutr1")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c^{2}").c_str());
  eventCont->hist("Pi0MassESEEetacutr1")->SetFillColor(kOrange);


  eventCont->add( new TH1F("Pi0MassESEEetacutr2","#gamma #gamma invariant mass corrected with Preshower energy for 1.7<|#eta|<2.0",histo_bins ,histo_min_range, histo_max_range) );
  eventCont->hist("Pi0MassESEEetacutr2")->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");
  eventCont->hist("Pi0MassESEEetacutr2")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c^{2}").c_str());
  eventCont->hist("Pi0MassESEEetacutr2")->SetFillColor(kOrange);


  eventCont->add( new TH1F("Pi0MassESEEetacutr3","#gamma #gamma invariant mass corrected with Preshower energy for 2.0<|#eta|<2.5",histo_bins ,histo_min_range, histo_max_range) );
  eventCont->hist("Pi0MassESEEetacutr3")->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");
  eventCont->hist("Pi0MassESEEetacutr3")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c^{2}").c_str());
  eventCont->hist("Pi0MassESEEetacutr3")->SetFillColor(kOrange);


  eventCont->add( new TH1F("Pi0MassESEEetacutr4","#gamma #gamma invariant mass corrected with Preshower energy for 2.2<|#eta|<2.8",histo_bins ,histo_min_range, histo_max_range) );
  eventCont->hist("Pi0MassESEEetacutr4")->GetXaxis()->SetTitle("#gamma#gamma invariant mass (GeV/c^{2})");
  eventCont->hist("Pi0MassESEEetacutr4")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c^{2}").c_str());
  eventCont->hist("Pi0MassESEEetacutr4")->SetFillColor(kOrange);
 
  
  histo_bins = 120;
  histo_min_range = 0. ;
  histo_max_range = 10.;
  bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins)*1000; //Mev
 
  strs<<setprecision(2)<<bin_size;
  s_amplitude = strs.str();
  strs.clear();
  strs.str("");
  
  eventCont->add(new TH1F("PtPi0EB", "P_{t} distribution of Pi0 for |#eta|<1.5", histo_bins ,histo_min_range, histo_max_range ) );
  eventCont->hist("PtPi0EB")->GetXaxis()->SetTitle("P_{t} of Pi0 (GeV/c)");
  eventCont->hist("PtPi0EB")->GetYaxis()->SetTitle(("entries/"+ s_amplitude+" MeV/c").c_str());
  eventCont->hist("PtPi0EB")->SetFillColor(kOrange);
  eventCont->hist("PtPi0EB")->SetStats(false);
  
  histo_bins = 120;
  histo_min_range = 0. ;
  histo_max_range = 10.;
  bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins)*1000; //Mev
  
  strs<<setprecision(2)<<bin_size;
  s_amplitude = strs.str();
  strs.clear();
  strs.str("");
  
  eventCont->add(new TH1F("PtPi0EE", "P_{t} distribution of Pi0 for 1.566<|#eta|<2.5",histo_bins, histo_min_range,  histo_max_range) );
  eventCont->hist("PtPi0EE")->GetXaxis()->SetTitle("P_{t} of Pi0 (GeV/c)");
  eventCont->hist("PtPi0EE")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c").c_str());
  eventCont->hist("PtPi0EE")->SetFillColor(kOrange);
  eventCont->hist("PtPi0EE")->SetStats(false);
  
  histo_bins = 100;
  histo_min_range = -2.5 ;
  histo_max_range = 2.5;
  bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
  
  strs<<setprecision(1)<<bin_size;
  s_amplitude = strs.str();
  strs.clear();
  strs.str("");
  
  eventCont->add(new TH1F("etaPi0EB", "#eta distribution of Pi0 |#eta|<1.5", histo_bins, histo_min_range,histo_max_range) );
  eventCont->hist("etaPi0EB")->GetXaxis()->SetTitle("#eta of Pi0");
  eventCont->hist("etaPi0EB")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
  eventCont->hist("etaPi0EB")->SetFillColor(kOrange);
  eventCont->hist("etaPi0EB")->SetStats(false);
  
  histo_bins = 100;
  histo_min_range = -3. ;
  histo_max_range = 3.;
  bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
  
  strs<<setprecision(1)<<bin_size;
  s_amplitude = strs.str();
  strs.clear();
  strs.str("");
  
  eventCont->add(new TH1F("etaPi0EE", "#eta distribution of Pi0 for 1.479<|#eta|<3.0", histo_bins, histo_min_range,histo_max_range) );
  eventCont->hist("etaPi0EE")->GetXaxis()->SetTitle("#eta of Pi0");
  eventCont->hist("etaPi0EE")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
  eventCont->hist("etaPi0EE")->SetFillColor(kOrange);
  eventCont->hist("etaPi0EE")->SetStats(false);
  
  histo_bins = 100;
  histo_min_range = -3.14 ;
  histo_max_range = 3.14;
  bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
  
  strs<<setprecision(1)<<bin_size;
  s_amplitude = strs.str();
  strs.clear();
  strs.str("");
  
  eventCont->add(new TH1F("phiPi0EB", "#phi distribution of Pi0 |#eta|<1.5", histo_bins, histo_min_range,histo_max_range) );
  eventCont->hist("phiPi0EB")->GetXaxis()->SetTitle("#phi of Pi0");
  eventCont->hist("phiPi0EB")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
  eventCont->hist("phiPi0EB")->SetFillColor(kOrange);
  eventCont->hist("phiPi0EB")->SetStats(false);
  
  histo_bins = 100;
  histo_min_range = -3.14 ;
  histo_max_range = 3.14;
  bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
  
  strs<<setprecision(1)<<bin_size;
  s_amplitude = strs.str();
  strs.clear();
  strs.str("");
  
  eventCont->add(new TH1F("phiPi0EE", "#phi distribution of Pi0 1.566<|#eta|<2.5", histo_bins, histo_min_range,histo_max_range) );
  eventCont->hist("phiPi0EE")->GetXaxis()->SetTitle("#phi of Pi0");
  eventCont->hist("phiPi0EE")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
  eventCont->hist("phiPi0EE")->SetFillColor(kOrange);
  eventCont->hist("phiPi0EE")->SetStats(false);
  
  histo_bins = 120;
  histo_min_range = 0. ;
  histo_max_range = 9.5;
  bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins)*1000; //Mev
  
  strs<<setprecision(2)<<bin_size;
  s_amplitude = strs.str();
  strs.clear();
  strs.str("");
  
  eventCont->add(new TH1F("Ptgamma1EB", "P_{t} distribution of first #gamma|#eta|<1.5 ",histo_bins, histo_min_range,histo_max_range) );
  eventCont->hist("Ptgamma1EB")->GetXaxis()->SetTitle("P_{t} of first #gamma (GeV/c)");
  eventCont->hist("Ptgamma1EB")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c").c_str());
  eventCont->hist("Ptgamma1EB")->SetFillColor(kRed);
  eventCont->hist("Ptgamma1EB")->SetStats(false);
  
  histo_bins = 100;
  histo_min_range = -1.6 ;
  histo_max_range = 1.6;
  bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
  
  strs<<setprecision(1)<<bin_size;
  s_amplitude = strs.str();
  strs.clear();
  strs.str("");
  
  eventCont->add(new TH1F("etagamma1EB", "#eta distribution of first #gamma |#eta|<1.5 ",histo_bins, histo_min_range,histo_max_range) );
  eventCont->hist("etagamma1EB")->GetXaxis()->SetTitle("#eta of first #gamma");
  eventCont->hist("etagamma1EB")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
  eventCont->hist("etagamma1EB")->SetFillColor(kRed);
  eventCont->hist("etagamma1EB")->SetStats(false);
  
  histo_bins = 100;
  histo_min_range = -3.14 ;
  histo_max_range = 3.14;
  bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
  
  strs<<setprecision(1)<<bin_size;
  s_amplitude = strs.str();
  strs.clear();
  strs.str("");
  
  eventCont->add(new TH1F("phigamma1EB", "#phi distribution of first #gamma |#eta|<1.5",histo_bins, histo_min_range,histo_max_range) );
  eventCont->hist("phigamma1EB")->GetXaxis()->SetTitle("#phi of first #gamma");
  eventCont->hist("phigamma1EB")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
  eventCont->hist("phigamma1EB")->SetFillColor(kRed);
  eventCont->hist("phigamma1EB")->SetStats(false);
  
  histo_bins = 120;
  histo_min_range = 0. ;
  histo_max_range = 9.5;
  bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins)*1000; //Mev
  
   strs<<setprecision(2)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");
   
   eventCont->add(new TH1F("Ptgamma2EB", "P_{t} distribution of second #gamma |#eta|<1.5", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("Ptgamma2EB")->GetXaxis()->SetTitle("P_{t} of second #gamma (GeV/c)");
   eventCont->hist("Ptgamma2EB")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c").c_str());
   eventCont->hist("Ptgamma2EB")->SetFillColor(kYellow);
   eventCont->hist("Ptgamma2EB")->SetStats(false);
   
   histo_bins = 100;
   histo_min_range = -1.6 ;
   histo_max_range = 1.6;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
   
   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");
   
   eventCont->add(new TH1F("etagamma2EB", "#eta distribution of second #gamma |#eta|<1.5", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("etagamma2EB")->GetXaxis()->SetTitle("#eta of second #gamma");
   eventCont->hist("etagamma2EB")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("etagamma2EB")->SetFillColor(kYellow);
   eventCont->hist("etagamma2EB")->SetStats(false);
   
   histo_bins = 100;
   histo_min_range = -3.14 ;
   histo_max_range = 3.14;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
   
   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");
   
   eventCont->add(new TH1F("phigamma2EB", "#phi distribution of second #gamma|#eta|<1.5 ", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("phigamma2EB")->GetXaxis()->SetTitle("#phi of second #gamma");
   eventCont->hist("phigamma2EB")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("phigamma2EB")->SetFillColor(kYellow);
   eventCont->hist("phigamma2EB")->SetStats(false);
   
   histo_bins = 120;
   histo_min_range = 0. ;
   histo_max_range = 9.5;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins)*1000; //MeV
   
   strs<<setprecision(2)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");
   
   eventCont->add(new TH1F("Ptgamma1EE", "P_{t} distribution of first #gamma 1.479<|#eta|<2.5", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("Ptgamma1EE")->GetXaxis()->SetTitle("P_{t} of first #gamma (GeV/c)");
   eventCont->hist("Ptgamma1EE")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c").c_str());
   eventCont->hist("Ptgamma1EE")->SetFillColor(kRed);
   eventCont->hist("Ptgamma1EE")->SetStats(false);
   
   histo_bins = 120;
   histo_min_range = 0. ;
   histo_max_range = 9.5;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins)*1000; //MeV
   
   strs<<setprecision(2)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");
   
   eventCont->add(new TH1F("Ptgamma1EEr1", "P_{t} distribution of first #gamma 1.479<|#eta|<2.0", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("Ptgamma1EEr1")->GetXaxis()->SetTitle("P_{t} of first #gamma (GeV/c)");
   eventCont->hist("Ptgamma1EEr1")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c").c_str());
   eventCont->hist("Ptgamma1EEr1")->SetFillColor(kRed);
   eventCont->hist("Ptgamma1EEr1")->SetStats(false);
   
   histo_bins = 120;
   histo_min_range = 0. ;
   histo_max_range = 9.5;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins)*1000; //MeV
   
   strs<<setprecision(2)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");
   
   eventCont->add(new TH1F("Ptgamma1EEr2", "P_{t} distribution of first #gamma 2.0<|#eta|<2.5",  histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("Ptgamma1EEr2")->GetXaxis()->SetTitle("P_{t} of first #gamma (GeV/c)");
   eventCont->hist("Ptgamma1EEr2")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c").c_str());
   eventCont->hist("Ptgamma1EEr2")->SetFillColor(kRed);
   eventCont->hist("Ptgamma1EEr2")->SetStats(false);
   
   histo_bins = 120;
   histo_min_range = 0. ;
   histo_max_range = 9.5;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins)*1000; //MeV
   
   strs<<setprecision(2)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");
   
   eventCont->add(new TH1F("Ptgamma1EEr3", "P_{t} distribution of first #gamma 2.5<|#eta|<3.0", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("Ptgamma1EEr3")->GetXaxis()->SetTitle("P_{t} of first #gamma (GeV/c)");
   eventCont->hist("Ptgamma1EEr3")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c").c_str());
   eventCont->hist("Ptgamma1EEr3")->SetFillColor(kRed);
   eventCont->hist("Ptgamma1EEr3")->SetStats(false);
   
   histo_bins = 100;
   histo_min_range = -3. ;
   histo_max_range = 3. ;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
   
   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");
   
   eventCont->add(new TH1F("etagamma1EE", "#eta distribution of first #gamma 1.566<|#eta|<3.0", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("etagamma1EE")->GetXaxis()->SetTitle("#eta of first #gamma");
   eventCont->hist("etagamma1EE")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("etagamma1EE")->SetFillColor(kRed);
   eventCont->hist("etagamma1EE")->SetStats(false);
   
   histo_bins = 100;
   histo_min_range = -3.14 ;
   histo_max_range = 3.14 ;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
   
   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");
   
   eventCont->add(new TH1F("phigamma1EE", "#phi distribution of first #gamma 1.566<|#eta|<2.5 ", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("phigamma1EE")->GetXaxis()->SetTitle("#phi of first #gamma");
   eventCont->hist("phigamma1EE")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("phigamma1EE")->SetFillColor(kRed);
   eventCont->hist("phigamma1EE")->SetStats(false);
   
   histo_bins = 120;
   histo_min_range = 0. ;
   histo_max_range = 9.5;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins)*1000; //MeV
   
   strs<<setprecision(2)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");
   
   eventCont->add(new TH1F("Ptgamma2EE", "P_{t} distribution of second #gamma 1.566<|#eta|<2.5", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("Ptgamma2EE")->GetXaxis()->SetTitle("P_{t} of second #gamma (GeV/c)");
   eventCont->hist("Ptgamma2EE")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c").c_str());
   eventCont->hist("Ptgamma2EE")->SetFillColor(kYellow);
   eventCont->hist("Ptgamma2EE")->SetStats(false);
   
   eventCont->add(new TH1F("Ptgamma2EE+", "P_{t} distribution of second #gamma 1.566< #eta <2.5", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("Ptgamma2EE+")->GetXaxis()->SetTitle("P_{t} of second #gamma (GeV/c)");
   eventCont->hist("Ptgamma2EE+")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c").c_str());
   eventCont->hist("Ptgamma2EE+")->SetFillColor(kYellow);
   eventCont->hist("Ptgamma2EE+")->SetStats(false);
   
   eventCont->add(new TH1F("Ptgamma2EE-", "P_{t} distribution of second #gamma -1.566< #eta < -2.5", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("Ptgamma2EE-")->GetXaxis()->SetTitle("P_{t} of second #gamma (GeV/c)");
   eventCont->hist("Ptgamma2EE-")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c").c_str());
   eventCont->hist("Ptgamma2EE-")->SetFillColor(kYellow);
   eventCont->hist("Ptgamma2EE-")->SetStats(false);
   
   histo_bins = 120;
   histo_min_range = 0. ;
   histo_max_range = 9.5;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins)*1000; //MeV
   
   strs<<setprecision(2)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");

   eventCont->add(new TH1F("Ptgamma2EEr1", "P_{t} distribution of second #gamma 1.479<|#eta|<2.0",  histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("Ptgamma2EEr1")->GetXaxis()->SetTitle("P_{t} of first #gamma (GeV/c)");
   eventCont->hist("Ptgamma2EEr1")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c").c_str());
   eventCont->hist("Ptgamma2EEr1")->SetFillColor(kYellow);
   eventCont->hist("Ptgamma2EEr1")->SetStats(false);

   histo_bins = 120;
   histo_min_range = 0.5 ;
   histo_max_range = 9.5;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins)*1000; //MeV
   
   strs<<setprecision(2)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");

   eventCont->add(new TH1F("Ptgamma2EEr2", "P_{t} distribution of second #gamma 2.0<|#eta|<2.5", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("Ptgamma2EEr2")->GetXaxis()->SetTitle("P_{t} of first #gamma (GeV/c)");
   eventCont->hist("Ptgamma2EEr2")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c").c_str());
   eventCont->hist("Ptgamma2EEr2")->SetFillColor(kYellow);
   eventCont->hist("Ptgamma2EEr2")->SetStats(false);

   histo_bins = 120;
   histo_min_range = 0. ;
   histo_max_range = 9.5;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins)*1000; //MeV
   
   strs<<setprecision(2)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");

   eventCont->add(new TH1F("Ptgamma2EEr3", "P_{t} distribution of second #gamma 2.5<|#eta|<3.0", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("Ptgamma2EEr3")->GetXaxis()->SetTitle("P_{t} of first #gamma (GeV/c)");
   eventCont->hist("Ptgamma2EEr3")->GetYaxis()->SetTitle(("entries/"+s_amplitude+" MeV/c").c_str());
   eventCont->hist("Ptgamma2EEr3")->SetFillColor(kYellow);
   eventCont->hist("Ptgamma2EEr3")->SetStats(false);

   histo_bins = 100;
   histo_min_range = -3. ;
   histo_max_range = 3. ;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
   
   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");

   eventCont->add(new TH1F("etagamma2EE", "#eta distribution of second #gamma 1.566<|#eta|<2.5", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("etagamma2EE")->GetXaxis()->SetTitle("#eta of second #gamma");
   eventCont->hist("etagamma2EE")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("etagamma2EE")->SetFillColor(kYellow);
   eventCont->hist("etagamma2EE")->SetStats(false);

   histo_bins = 100;
   histo_min_range = -3.14 ;
   histo_max_range = 3.14 ;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
   
   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");

   eventCont->add(new TH1F("phigamma2EE", "#phi distribution of second #gamma 1.566<|#eta|<2.5", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("phigamma2EE")->GetXaxis()->SetTitle("#phi of second #gamma");
   eventCont->hist("phigamma2EE")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("phigamma2EE")->SetFillColor(kYellow);
   eventCont->hist("phigamma2EE")->SetStats(false);

   histo_bins = 50;
   histo_min_range = 0.8 ;
   histo_max_range = 1.15 ;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
   
   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");

   eventCont->add(new TH1F("s4s9S_EB", "s4/s9 in the signal region (0.95 - 0.171 Gev/c^{2}) |#eta|<1.5", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("s4s9S_EB")->GetXaxis()->SetTitle("s4/s9");
   eventCont->hist("s4s9S_EB")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("s4s9S_EB")->SetFillColor(kOrange);
   eventCont->hist("s4s9S_EB")->SetStats(false);

   histo_bins = 50;
   histo_min_range = 0.8 ;
   histo_max_range = 1.15 ;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
   
   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");

   eventCont->add(new TH1F("s4s9SB_EB", "s4/s9 in the sideband region (0.2 - 0.3 Gev/c^{2}) |#eta|<1.5", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("s4s9SB_EB")->GetXaxis()->SetTitle("s4/s9");
   eventCont->hist("s4s9SB_EB")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("s4s9SB_EB")->SetFillColor(kViolet);
   eventCont->hist("s4s9SB_EB")->SetStats(false);

   histo_bins = 50;
   histo_min_range = 0.9 ;
   histo_max_range = 1.15 ;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
   
   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");

   eventCont->add(new TH1F("s4s9S_EE", "s4/s9 in the signal region (0.95 - 0.171 Gev/c^{2}) 1.566<|#eta|<2.5", histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("s4s9S_EE")->GetXaxis()->SetTitle("s4/s9");
   eventCont->hist("s4s9S_EE")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("s4s9S_EE")->SetFillColor(kOrange);
   eventCont->hist("s4s9S_EE")->SetStats(false);

   histo_bins = 50;
   histo_min_range = 0.9 ;
   histo_max_range = 1.15 ;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
   
   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");

   eventCont->add(new TH1F("s4s9SB_EE", "s4/s9 in the sideband region (0.2 - 0.3 Gev/c^{2}) 1.566<|#eta|<2.5",  histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("s4s9SB_EE")->GetXaxis()->SetTitle("s4/s9");
   eventCont->hist("s4s9SB_EE")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("s4s9SB_EE")->SetFillColor(kViolet);
   eventCont->hist("s4s9SB_EE")->SetStats(false);

   histo_bins = 360;
   histo_min_range = 0. ;
   histo_max_range = 7. ;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);
   
   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");

   eventCont->add(new TH1F("PreshowerClusterEnergy", "Energy gathered from valid matched Preshower Clusters",  histo_bins, histo_min_range,histo_max_range) );
   eventCont->hist("PreshowerClusterEnergy")->GetXaxis()->SetTitle("Energy [GeV]");
   eventCont->hist("PreshowerClusterEnergy")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("PreshowerClusterEnergy")->SetFillColor(kViolet);
   eventCont->hist("PreshowerClusterEnergy")->SetStats(false);

   eventCont->add(
    new TH2F("BlindstripPositionZpos","Position of the strips not matched for Z+",1200,-130.,130.,1200,-130., 130.)
   );
   eventCont->hist("BlindstripPositionZpos")->GetXaxis()->SetTitle("X Axis Position [cm]");
   eventCont->hist("BlindstripPositionZpos")->GetYaxis()->SetTitle("Y Axis Position [cm]");
   eventCont->hist("BlindstripPositionZpos")->SetStats(false);
   eventCont->hist("BlindstripPositionZpos")->SetFillColor(kWhite);

   eventCont->add(
    new TH2F("BlindstripPositionZneg","Position of the strips not matched for Z-",1200,-130.,130.,1200,-130., 130.)
   );
   eventCont->hist("BlindstripPositionZneg")->GetXaxis()->SetTitle("X Axis Position [cm]");
   eventCont->hist("BlindstripPositionZneg")->GetYaxis()->SetTitle("Y Axis Position [cm]");
   eventCont->hist("BlindstripPositionZneg")->SetStats(false);
   eventCont->hist("BlindstripPositionZneg")->SetFillColor(kWhite);

   eventCont->add(
    new TH2F("ESclusterPositionZpos","Position of ESclusters for Z+",1200,-130.,130.,1200,-130., 130.)
   );
   eventCont->hist("ESclusterPositionZpos")->GetXaxis()->SetTitle("X Axis Position [cm]");
   eventCont->hist("ESclusterPositionZpos")->GetYaxis()->SetTitle("Y Axis Position [cm]");
   eventCont->hist("ESclusterPositionZpos")->SetStats(false);
   eventCont->hist("ESclusterPositionZpos")->SetFillColor(kWhite);

   eventCont->add(
    new TH2F("ESclusterPositionZneg","Position of ESclusters for Z-",1200,-130.,130.,1200,-130., 130.)
   );
   eventCont->hist("ESclusterPositionZneg")->GetXaxis()->SetTitle("X Axis Position [cm]");
   eventCont->hist("ESclusterPositionZneg")->GetYaxis()->SetTitle("Y Axis Position [cm]");
   eventCont->hist("ESclusterPositionZneg")->SetStats(false);
   eventCont->hist("ESclusterPositionZneg")->SetFillColor(kWhite);


   eventCont->add(
    new TH2F("nonBlindstripPositionZpos","Position of the strips not matched for Z+",1200,-130.,130.,1200,-130., 130.)
   );
   eventCont->hist("nonBlindstripPositionZpos")->GetXaxis()->SetTitle("X Axis Position [cm]");
   eventCont->hist("nonBlindstripPositionZpos")->GetYaxis()->SetTitle("Y Axis Position [cm]");
   eventCont->hist("nonBlindstripPositionZpos")->SetStats(false);
   eventCont->hist("nonBlindstripPositionZpos")->SetFillColor(kWhite);

   eventCont->add(
    new TH2F("nonBlindstripPositionZneg","Position of the strips not matched for Z-",1200,-130.,130.,1200,-130., 130.)
   );
   eventCont->hist("nonBlindstripPositionZneg")->GetXaxis()->SetTitle("X Axis Position [cm]");
   eventCont->hist("nonBlindstripPositionZneg")->GetYaxis()->SetTitle("Y Axis Position [cm]");
   eventCont->hist("nonBlindstripPositionZneg")->SetStats(false);
   eventCont->hist("nonBlindstripPositionZneg")->SetFillColor(kWhite);


   eventCont->add(
    new TH2F("PreshowerEnergy_vs_Etapos","Preshower Energy vs #eta(+) Preshower Cluster",360,1.5,2.4,360,0., 7.)
   );
   eventCont->hist("PreshowerEnergy_vs_Etapos")->GetXaxis()->SetTitle("#eta Preshower Cluster");
   eventCont->hist("PreshowerEnergy_vs_Etapos")->GetYaxis()->SetTitle("Preshower Cluster Energy [GeV]");
   eventCont->hist("PreshowerEnergy_vs_Etapos")->SetStats(false);
   eventCont->hist("PreshowerEnergy_vs_Etapos")->SetFillColor(kWhite);

   eventCont->add(
    new TH2F("PreshowerEnergy_vs_Etaneg","Preshower Energy vs #eta(-) Preshower Cluster",360,-2.4,-1.5,360,0., 7.)
   );
   eventCont->hist("PreshowerEnergy_vs_Etaneg")->GetXaxis()->SetTitle("#eta Preshower Cluster");
   eventCont->hist("PreshowerEnergy_vs_Etaneg")->GetYaxis()->SetTitle("Preshower Cluster Energy [GeV]");
   eventCont->hist("PreshowerEnergy_vs_Etaneg")->SetStats(false);
   eventCont->hist("PreshowerEnergy_vs_Etaneg")->SetFillColor(kWhite);


   eventCont->add(
    new TH2F("PreshowerEnergy_vs_EEenergy","Preshower Energy vs Endcap Energy",360,0.,20.,360,0., 7.)
   );
   eventCont->hist("PreshowerEnergy_vs_EEenergy")->GetXaxis()->SetTitle("Endcap Cluster Energy [GeV]");
   eventCont->hist("PreshowerEnergy_vs_EEenergy")->GetYaxis()->SetTitle("Preshower Cluster Energy [GeV]");
   eventCont->hist("PreshowerEnergy_vs_EEenergy")->SetStats(false);
   eventCont->hist("PreshowerEnergy_vs_EEenergy")->SetFillColor(kWhite);

   eventCont->add(
    new TH2F("PreshowerEnergy_vs_TotEnergy","Preshower Cluster Energy vs (EE+ES) Cluster Energy",360,0.,30.,360,0., 7.)
   );
   eventCont->hist("PreshowerEnergy_vs_TotEnergy")->GetXaxis()->SetTitle("Preshower Cluster + Endcap Cluster Energy [GeV]");
   eventCont->hist("PreshowerEnergy_vs_TotEnergy")->GetYaxis()->SetTitle("Preshower Cluster Energy [GeV]");
   eventCont->hist("PreshowerEnergy_vs_TotEnergy")->SetStats(false);
   eventCont->hist("PreshowerEnergy_vs_TotEnergy")->SetFillColor(kWhite);

  eventCont->add(
    new TH2F("EEcluster_ESstrip_distance_vs_EEclusterEta_2D_Etapos","|EEclusterPosition - EsMatchedStripPosition|(2D) vs (+)#eta EEcluster",360,1.6,2.1,360,0., 15.)
   );
   eventCont->hist("EEcluster_ESstrip_distance_vs_EEclusterEta_2D_Etapos")->GetXaxis()->SetTitle("#eta EEcluster");
   eventCont->hist("EEcluster_ESstrip_distance_vs_EEclusterEta_2D_Etapos")->GetYaxis()->SetTitle("|EEclusterPosition - EsMatchedStripPosition| 2D [cm]");
   eventCont->hist("EEcluster_ESstrip_distance_vs_EEclusterEta_2D_Etapos")->SetStats(false);
   eventCont->hist("EEcluster_ESstrip_distance_vs_EEclusterEta_2D_Etapos")->SetFillColor(kWhite);

  eventCont->add(
    new TH2F("EEcluster_ESstrip_distance_vs_EEclusterEta_2D_Etaneg","|EEclusterPosition - EsMatchedStripPosition|(2D) vs (-)#eta EEcluster",360,-2.1,-1.6,360,0., 15.)
   );
   eventCont->hist("EEcluster_ESstrip_distance_vs_EEclusterEta_2D_Etaneg")->GetXaxis()->SetTitle("#eta EEcluster");
   eventCont->hist("EEcluster_ESstrip_distance_vs_EEclusterEta_2D_Etaneg")->GetYaxis()->SetTitle("|EEclusterPosition - EsMatchedStripPosition| 2D [cm]");
   eventCont->hist("EEcluster_ESstrip_distance_vs_EEclusterEta_2D_Etaneg")->SetStats(false);
   eventCont->hist("EEcluster_ESstrip_distance_vs_EEclusterEta_2D_Etaneg")->SetFillColor(kWhite);

   eventCont->add(
    new TH2F("EEcluster_ESstrip_distance_vs_EEclusterEta_3D_Etapos","|EEclusterPosition - EsMatchedStripPosition|(3D) vs (+)#eta EEcluster",360,1.6,2.1,360,0., 30.)
   );
   eventCont->hist("EEcluster_ESstrip_distance_vs_EEclusterEta_3D_Etapos")->GetXaxis()->SetTitle("#eta EEcluster");
   eventCont->hist("EEcluster_ESstrip_distance_vs_EEclusterEta_3D_Etapos")->GetYaxis()->SetTitle("|EEclusterPosition - EsMatchedStripPosition| 3D [cm]");
   eventCont->hist("EEcluster_ESstrip_distance_vs_EEclusterEta_3D_Etapos")->SetStats(false);
   eventCont->hist("EEcluster_ESstrip_distance_vs_EEclusterEta_3D_Etapos")->SetFillColor(kWhite);

  eventCont->add(
    new TH2F("EEcluster_ESstrip_distance_vs_EEclusterEta_3D_Etaneg","|EEclusterPosition - EsMatchedStripPosition|(3D) vs (-)#eta EEcluster",360,-2.1,-1.6,360,0., 30.)

   );
   eventCont->hist("EEcluster_ESstrip_distance_vs_EEclusterEta_3D_Etaneg")->GetXaxis()->SetTitle("#eta EEcluster");
   eventCont->hist("EEcluster_ESstrip_distance_vs_EEclusterEta_3D_Etaneg")->GetYaxis()->SetTitle("|EEclusterPosition - EsMatchedStripPosition| 3D [cm]");
   eventCont->hist("EEcluster_ESstrip_distance_vs_EEclusterEta_3D_Etaneg")->SetStats(false);
   eventCont->hist("EEcluster_ESstrip_distance_vs_EEclusterEta_3D_Etaneg")->SetFillColor(kWhite);


   histo_bins = 720;
   histo_min_range = 0. ;
   histo_max_range = 0.25 ;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);

   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");

   eventCont->add(
    new TH1F("ES_EE_DeltaR","#DeltaR between the ES Cluster and th EE Cluster",histo_bins,histo_min_range,histo_max_range)
   );
   eventCont->hist("ES_EE_DeltaR")->GetXaxis()->SetTitle("#DeltaR [rad]");
   eventCont->hist("ES_EE_DeltaR")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("ES_EE_DeltaR")->SetFillColor(kViolet);
   eventCont->hist("ES_EE_DeltaR")->SetStats(false);
   eventCont->hist("ES_EE_DeltaR")->SetFillColor(kWhite);

   histo_bins = 720;
   histo_min_range = 0. ;
   histo_max_range = 15. ;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);

   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");

   eventCont->add(
    new TH1F("ES_EE_Deltar_2D","Distance between the ES Cluster and th EE Cluster inside X-Y plane",histo_bins,histo_min_range,histo_max_range)
   );
   eventCont->hist("ES_EE_Deltar_2D")->GetXaxis()->SetTitle("#Deltar [cm]");
   eventCont->hist("ES_EE_Deltar_2D")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("ES_EE_Deltar_2D")->SetFillColor(kViolet);
   eventCont->hist("ES_EE_Deltar_2D")->SetStats(false);
   eventCont->hist("ES_EE_Deltar_2D")->SetFillColor(kWhite);

   histo_bins = 720;
   histo_min_range = 0. ;
   histo_max_range = 30. ;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);

   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");

   eventCont->add(
    new TH1F("ES_EE_Deltar_3D","Distance between the ES Cluster and th EE Cluster (3D)",histo_bins,histo_min_range,histo_max_range)
   );
   eventCont->hist("ES_EE_Deltar_3D")->GetXaxis()->SetTitle("#Deltar (3D) [cm]");
   eventCont->hist("ES_EE_Deltar_3D")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("ES_EE_Deltar_3D")->SetFillColor(kViolet);
   eventCont->hist("ES_EE_Deltar_3D")->SetStats(false);
   eventCont->hist("ES_EE_Deltar_3D")->SetFillColor(kWhite);


   histo_bins = 720;
   histo_min_range = 0. ;
   histo_max_range = 15. ;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);

   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");

   eventCont->add(
    new TH1F("ESseed_EE_Deltar_2D","Distance between the ES Window seed and th EE Cluster inside X-Y plane",histo_bins,histo_min_range,histo_max_range)
   );
   eventCont->hist("ESseed_EE_Deltar_2D")->GetXaxis()->SetTitle("#Deltar [cm]");
   eventCont->hist("ESseed_EE_Deltar_2D")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("ESseed_EE_Deltar_2D")->SetFillColor(kViolet);
   eventCont->hist("ESseed_EE_Deltar_2D")->SetStats(false);
   eventCont->hist("ESseed_EE_Deltar_2D")->SetFillColor(kWhite);

   histo_bins = 360;
   histo_min_range = 8. ;
   histo_max_range = 30. ;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);

   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");

   eventCont->add(
    new TH1F("deltaZcrosscheck","#DeltaZ just to see 2 spikes",histo_bins,histo_min_range,histo_max_range)
   );
   eventCont->hist("deltaZcrosscheck")->GetXaxis()->SetTitle("#Deltaz [cm]");
   eventCont->hist("deltaZcrosscheck")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("deltaZcrosscheck")->SetFillColor(kViolet);
   eventCont->hist("deltaZcrosscheck")->SetStats(false);
   eventCont->hist("deltaZcrosscheck")->SetFillColor(kWhite);


   histo_bins = 720;
   histo_min_range = 0. ;
   histo_max_range = 30. ;
   bin_size = (double)(( histo_max_range - histo_min_range)/  histo_bins);

   strs<<setprecision(1)<<bin_size;
   s_amplitude = strs.str();
   strs.clear();
   strs.str("");

   eventCont->add(
    new TH1F("ESseed_EE_Deltar_3D","#Distance between the ES Window seed and th EE Cluster",histo_bins,histo_min_range,histo_max_range)
   );
   eventCont->hist("ESseed_EE_Deltar_3D")->GetXaxis()->SetTitle("#Deltar [cm]");
   eventCont->hist("ESseed_EE_Deltar_3D")->GetYaxis()->SetTitle(("entries/"+s_amplitude).c_str());
   eventCont->hist("ESseed_EE_Deltar_3D")->SetFillColor(kViolet);
   eventCont->hist("ESseed_EE_Deltar_3D")->SetStats(false);
   eventCont->hist("ESseed_EE_Deltar_3D")->SetFillColor(kWhite);


   eventCont->add(
    new TProfile2D("ESclus_Avg_Energy_Zpos","ESCluster average energy map for Z+",1200,-130.,130.,1200,-130., 130.)
   );
   eventCont->hist("ESclus_Avg_Energy_Zpos")->GetXaxis()->SetTitle("X Axis Position [cm]");
   eventCont->hist("ESclus_Avg_Energy_Zpos")->GetYaxis()->SetTitle("Y Axis Position [cm]");
   eventCont->hist("ESclus_Avg_Energy_Zpos")->SetStats(false);
   eventCont->hist("ESclus_Avg_Energy_Zpos")->SetFillColor(kWhite);

   eventCont->add(
    new TProfile2D("ESclus_Avg_Energy_Zneg","ESCluster average energy map for Z-",1200,-130.,130.,1200,-130., 130.)
   );
   eventCont->hist("ESclus_Avg_Energy_Zneg")->GetXaxis()->SetTitle("X Axis Position [cm]");
   eventCont->hist("ESclus_Avg_Energy_Zneg")->GetYaxis()->SetTitle("Y Axis Position [cm]");
   eventCont->hist("ESclus_Avg_Energy_Zneg")->SetStats(false);
   eventCont->hist("ESclus_Avg_Energy_Zneg")->SetFillColor(kWhite);


   eventCont->add( 
    new TH2F("Pi0MassvsphiposEB","#gamma #gamma invariant mass vs. #phi in EB+ |#eta|<1",72,-3.14,3.14,120,0.050, 0.200)
   );
   eventCont->add( 
    new TH2F("Pi0MassvsphinegEB","#gamma #gamma invariant mass vs. #phi in EB- |#eta|<1",72,-3.14,3.14,120,0.050, 0.200)
   );
   eventCont->hist("Pi0MassvsphiposEB")->GetXaxis()->SetTitle("#phi");
   eventCont->hist("Pi0MassvsphiposEB")->GetYaxis()->SetTitle("#gamma#gamma invariant mass");
   eventCont->hist("Pi0MassvsphiposEB")->SetStats(false);
   eventCont->hist("Pi0MassvsphiposEB")->SetFillColor(kWhite);
   eventCont->hist("Pi0MassvsphinegEB")->GetXaxis()->SetTitle("#phi");
   eventCont->hist("Pi0MassvsphinegEB")->GetYaxis()->SetTitle("#gamma#gamma invariant mass");
   eventCont->hist("Pi0MassvsphinegEB")->SetStats(false);
   eventCont->hist("Pi0MassvsphinegEB")->SetFillColor(kWhite);

   const int ptBins = 68;
   double ptEdge[ptBins] = {0.0,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
                         1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                         2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
                         3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9,
                         4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9,
                         5.0, 5.2, 5.4, 5.6, 5.8,
                         6.0, 6.5, 7.0, 7.5, 8.0,
                         10., 12., 14., 16., 18., 
                         20., 25., 30., 35., 100.}; // 43 filled 
   eventCont->add( 
    new TH2F("Pi0Massvsptpi0EB","#gamma #gamma invariant mass vs. p_{T}(#pi^{0}) for |#eta|<1",ptBins-1,ptEdge,120,0.050, 0.200)
   );
   eventCont->hist("Pi0Massvsptpi0EB")->GetXaxis()->SetTitle("p_{T}(#pi^{0}) (GeV/c)");
   eventCont->hist("Pi0Massvsptpi0EB")->GetYaxis()->SetTitle("#gamma#gamma invariant mass");
   eventCont->hist("Pi0Massvsptpi0EB")->SetStats(false);
   eventCont->hist("Pi0Massvsptpi0EB")->SetFillColor(kWhite);


   // ===== begin : too many for xtal calibration - data stored in TTree
   /**
   char name[100], title[100];
   for(uint32_t iR=0; iR<Type::nRegions; ++iR) { 

      // pi0 mass in each region after cut
      sprintf(name,"Pi0Mass_%u",iR);
      sprintf(title,"#gamma #gamma invariant mass in calibration region %u",iR);
      eventCont->add( new TH1F(name,title,120,0.050, 0.200) );

      // weighted eps for each region
      sprintf(name,"eps_%u",iR);
      sprintf(title,"<#epsilon> #times w in calibration region %u",iR);
      eventCont->add( new TH1F(name,title,360,-0.4,0.4) );
   } // loop over regions
   **/
   // ===== end : too many for xtal calibration - data stored in TTree

}


   
//============================================
template<class Type,int NMaxIter>
void EcalCalibAlgo<Type,NMaxIter>::makeSummaryHisto() {
//============================================

   // reopen the file to add summary histo
   TFile* f = TFile::Open( parser->stringValue("outputFile").c_str(), "UPDATE" );

   // 2D calib map in the barrel
   TH2F* hmap =  makeEBXtalTH2("calibMap","EB calib coefficients: #eta on x, #phi on y");
   TH2F* hent =  makeEBXtalTH2("entriesMap","EB weighted entries: #eta on x, #phi on y");

   // eta on y and phi on x as usually shown at the calibration meeting
   TH2F* hmapstd =  makeEBXtalTH2Std("calibMapStd","EB calib coefficients: #phi on x, #eta on y");
   TH2F* hentstd =  makeEBXtalTH2Std("entriesMapStd","EB weighted entries: #phi on x, #eta on y");

   for(uint32_t j=0; j<Type::nRegions; ++j)  {
      std::vector<DetId> ids = Type::allDetIdsInRegion(j);
      //cout << "Region: " << j << " # xtals: " << ids.size() << endl;
      for(std::vector<DetId>::const_iterator iid = ids.begin();
                   			     iid != ids.end(); ++iid) {
         EBDetId ebid(*iid);
         //int ix = ebid.ieta()< 0 ? ebid.ieta()+EBDetId::MAX_IETA+1 : ebid.ieta()+EBDetId::MAX_IETA+2;
         int ix = ebid.ieta()+EBDetId::MAX_IETA+1;
         //cout << "iR: " << j << " id: " << ebid << " ix: " << ix
         //     << " entries: " << Denom[j] << " coeff: " << tmpcalib[j] << endl;

         hmap->SetBinContent( ix, ebid.iphi(), tmpcalib[j] );
         hent->SetBinContent( ix, ebid.iphi(), Denom[j] );
         hmapstd->SetBinContent( ebid.iphi(), ix, tmpcalib[j] );
         hentstd->SetBinContent( ebid.iphi(), ix, Denom[j] );
      } // loop over DetId in regions
   }
   hmap->SetMinimum(0.9);
   hmap->SetStats(false);
   hmap->Write();

   hmapstd->SetMinimum(0.9);
   hmapstd->SetStats(false);
   hmapstd->Write();

   hent->SetStats(false);
   hent->Write();

   hentstd->SetStats(false);
   hentstd->Write();

   f->mkdir("summary","plots vs iteration");
   f->cd("summary");
   char name[100],title[100];

      // calib coeff vs iter
// ========  begin: way too many for xtal based calibration - data stored in TTree
/**
   for(uint32_t iR=0; iR<Type::nRegions; ++iR) {
      sprintf(name,"coeff_%u",iR);
      sprintf(title,"calib coeff vs iteration in region %u",iR);
      TH1F* h = new TH1F(name,title,NMaxIter,0.5,float(NMaxIter)+0.5);
      h->GetXaxis()->SetTitle("# iteration");
      h->GetYaxis()->SetTitle("calibration coefficient");

      TH1F* h2 = new TH1F(*h);
      sprintf(name,"entries_%u",iR);
      sprintf(title,"weighted entries vs iteration in region %u",iR);
      h2->SetName(name);
      h2->SetTitle(title);
      // loop over iterations and fill each iter bin
      for(int j=0; j<NMaxIter; ++j) {
        if(entries[iR][j] == 0. ) continue; // skip iterations we haven't done
        h->SetBinContent(j+1, calibCoeff[iR][j]);
        h2->SetBinContent(j+1, entries[iR][j]);
      } 
      h->SetMinimum(0.8);
      h->Write();
      h2->Write();
   }
**/
// ========  end: way too many for xtal based calibration - data stored in TTree

   // calibration precision
   sprintf(name,"coeffmu_vs_iter");
   sprintf(title,"mean calib coeff vs iteration");
   TH1F* hm = new TH1F(name,title,NMaxIter,0.5,float(NMaxIter)+0.5);
   hm->GetXaxis()->SetTitle("# iteration");
   hm->GetYaxis()->SetTitle("mean calibration coefficient");

   TH1F* hmlim = new TH1F(*hm);
   hmlim->SetName("coeffmu_vs_iter_lim");
   hmlim->SetTitle("mean calib coeff vs iteration |i#eta|<=50");

   sprintf(name,"coeffrms_vs_iter");
   sprintf(title,"RMS of calib coeff vs iteration");
   TH1F* hs = new TH1F(name,title,NMaxIter,0.5,float(NMaxIter)+0.5);
   hs->GetXaxis()->SetTitle("# iteration");
   hs->GetYaxis()->SetTitle("RMS calibration coefficient");

   TH1F* hslim = new TH1F(*hm);
   hslim->SetName("coeffrms_vs_iter_lim");
   hslim->SetTitle("RMS calib coeff vs iteration |i#eta|<=50");

   sprintf(name,"coeff_precision_vs_iter");
   sprintf(title,"Precision of calib coeff vs iteration");
   TH1F* hp = new TH1F(name,title,NMaxIter,0.5,float(NMaxIter)+0.5);
   hp->GetXaxis()->SetTitle("# iteration");
   hp->GetYaxis()->SetTitle("RMS/#mu calibration coefficient");

   TH1F* hplim = new TH1F(*hm);
   hplim->SetName("coeff_precision_vs_iter_lim");
   hplim->SetTitle("Precision of calib coeff vs iteration |i#eta|<=50");


   // distribution of calb coeff after each iteration
   for(int it=0; it<NMaxIter; ++it) {
      // only fill plots for iterations we have run
      float ncand(0.);
      for(uint32_t iR=0; iR<Type::nRegions; ++iR) ncand += entries[iR][it];
      if(ncand == 0.) continue;

      sprintf(name,"coeff_iter_%u",it+1);
      sprintf(title,"Distribution of calib coeff after iteration %u",it+1);

      TH1F* h = new TH1F(name,title,360,0.75,1.15);
      h->GetXaxis()->SetTitle("# iteration");
      h->GetYaxis()->SetTitle("a.u.");

      sprintf(name,"coeff_lim_iter_%u",it+1);
      sprintf(title,"Distribution of calib coeff with |i#eta|<50 after iteration %u",it+1);

      TH1F* hlim = new TH1F(*h);
      hlim->SetName(name);
      hlim->SetTitle(title);

      for(uint32_t iR=0; iR<Type::nRegions; ++iR) {
         h->Fill( calibCoeff[iR][it] );
         if(Type::detIdFromRegion(iR).ietaAbs()<=50) hlim->Fill(calibCoeff[iR][it] );
      }
      h->Write();
      hlim->Write();

      hm->SetBinContent(it+1, h->GetMean());
      hs->SetBinContent(it+1, h->GetRMS());
      hp->SetBinContent(it+1, h->GetRMS()/h->GetMean());

      hmlim->SetBinContent(it+1, hlim->GetMean());
      hslim->SetBinContent(it+1, hlim->GetRMS());
      hplim->SetBinContent(it+1, hlim->GetRMS()/hlim->GetMean());
   }
   hm->Write();
   hs->Write();
   hp->Write();
   hmlim->Write();
   hslim->Write();
   hplim->Write();

   f->Close();
}

//============================================
template<class Type,int NMaxIter>
void EcalCalibAlgo<Type,NMaxIter>::loadCalibMapFromFile(const char* cfile) {
//============================================

  TFile* f = TFile::Open(cfile);

  TH2F* hmap = (TH2F*) f->Get("calibMap");
  if(!hmap) throw std::runtime_error("cannot find TH2F::calibMap in the file provided");

  cout << "loading constants from TH2F::calibMap in <" << cfile << "> ..." << endl;

  for(int ix=1; ix<= hmap->GetXaxis()->GetNbins(); ++ix) {
     if(ix==86) continue;
     for(int iphi=1; iphi<= hmap->GetYaxis()->GetNbins(); ++iphi) {
        EBDetId id(ix-1-85,iphi);
        //cout << "id: " << id << " calib coeff: " << hmap->GetBinContent(ix,iphi) << endl;
        calibMap->coeff(id) =  hmap->GetBinContent(ix,iphi);
     }
  }
  cout <<  "done." << endl;   
}


//============================================
template<class Type,int NMaxIter>
void EcalCalibAlgo<Type,NMaxIter>::loadCorrectionsFromFile(const char* cfile) {
//============================================

  TFile* f = TFile::Open(cfile);
  cout << "loading corrections from <" << cfile << "> ..." << endl;

  TF1* fn = (TF1*) f->Get("NEnergyBins");
  if(!fn) throw std::runtime_error("cannot find TF1::NEnergyBins in ths file");

  cout << "Energy correction functions vs eta to be read for "
       << int(fn->Eval(0.)) << " bins of energy" << endl;

  char name[100];
  for(int i=0; i< int(fn->Eval(0.)); ++i) {
    sprintf(name,"f_corr_vs_eta_en_%d",i);
    TF1* ftemp = (TF1*)f->Get(name);
    if(!ftemp) {
      cout << "cannot find correction function " << name << endl; 
      continue;
    }
    corrFunc[i] = new TF1(*(TF1*)f->Get(name));
    //corrFunc[i]->Print();
  }
  f->Close();

  noCorrections_ = false; // we do have corrections to apply
  cout << "done." << endl;
}


//============================================
template<class Type,int NMaxIter>
float EcalCalibAlgo<Type,NMaxIter>::energyCorrection(float energy, float eta) {
//============================================
   if(energy<0.5 || energy>6.) return 1.;
   if(fabs(eta)>1.1) return 1.;

   int ien;
   for(ien=5; ien< 40; ++ien) { // BAAAAD. remove hardcoded value
     if(energy<= lowEnergyEdge[ien+1]) break;
   }
   return 1./corrFunc[ien]->Eval(eta);

}

//============================================
template<class Type,int NMaxIter>
void EcalCalibAlgo<Type,NMaxIter>::fillLowEnergyEdge() {
//============================================
   lowEnergyEdge[0] = 0. ;
   lowEnergyEdge[1] = 0.1 ;
   lowEnergyEdge[2] = 0.2 ;
   lowEnergyEdge[3] = 0.3 ;
   lowEnergyEdge[4] = 0.4;
   lowEnergyEdge[5] = 0.5;
   lowEnergyEdge[6] = 0.6;
   lowEnergyEdge[7] = 0.7;
   lowEnergyEdge[8] = 0.8;
   lowEnergyEdge[9] = 0.9;
   lowEnergyEdge[10] = 1.0;
   lowEnergyEdge[11] = 1.1;
   lowEnergyEdge[12] = 1.2;
   lowEnergyEdge[13] = 1.3;
   lowEnergyEdge[14] = 1.4;
   lowEnergyEdge[15] = 1.5;
   lowEnergyEdge[16] = 1.6;
   lowEnergyEdge[17] = 1.7;
   lowEnergyEdge[18] = 1.8;
   lowEnergyEdge[19] = 1.9;
   lowEnergyEdge[20] = 2.0;
   lowEnergyEdge[21] = 2.1;
   lowEnergyEdge[22] = 2.2;
   lowEnergyEdge[23] = 2.3;
   lowEnergyEdge[24] = 2.4;
   lowEnergyEdge[25] = 2.5;
   lowEnergyEdge[26] = 2.6;
   lowEnergyEdge[27] = 2.7;
   lowEnergyEdge[28] = 2.8;
   lowEnergyEdge[29] = 2.9;
   lowEnergyEdge[30] = 3.0;
   lowEnergyEdge[31] = 3.1;
   lowEnergyEdge[32] = 3.2;
   lowEnergyEdge[33] = 3.3;
   lowEnergyEdge[34] = 3.4;
   lowEnergyEdge[35] = 3.5;
   lowEnergyEdge[36] = 4.0;
   lowEnergyEdge[37] = 5.0;
   lowEnergyEdge[38] = 6.0;
   lowEnergyEdge[39] = 7.0;
   lowEnergyEdge[40] = 8.0;
   lowEnergyEdge[41] = 9.0;
   lowEnergyEdge[42] = 10.0;
   lowEnergyEdge[43] = 12.0;
   lowEnergyEdge[44] = 16.0;
   lowEnergyEdge[45] = 24.0;
   lowEnergyEdge[46] = 40.0;
}


//============================================
template<class Type,int NMaxIter>
void EcalCalibAlgo<Type,NMaxIter>::fillResultTree() {
//============================================

   TFile* f = TFile::Open( parser->stringValue("outputFile").c_str(), "UPDATE" );

   TTree* tree = new TTree("calib","Tree of EB Inter-calibration constants");

   uint32_t   rawId;
   int        hashedIndex;
   int        ieta;
   int        iphi;
   int        iSM;
   int        iMod;
   int        iTT;
   int        iTTeta;
   int        iTTphi;
   int        iter;
   float      regCoeff[NMaxIter];
   float      regEntries[NMaxIter];


   // GENERAL block branches
   tree->Branch("rawId",&rawId,"rawId/i");
   tree->Branch("hashedIndex",&hashedIndex,"hashedIndex/I");
   tree->Branch("ieta",&ieta,"ieta/I");
   tree->Branch("iphi",&iphi,"iphi/I");
   tree->Branch("iSM",&iSM,"iSM/I");
   tree->Branch("iMod",&iMod,"iMod/I");
   tree->Branch("iTT",&iTT,"iTT/I");
   tree->Branch("iTTeta",&iTTeta,"iTTeta/I");
   tree->Branch("iTTphi",&iTTphi,"iTTphi/I");
   tree->Branch("iter",&iter,"iter/I");
   tree->Branch("coeff",&regCoeff,"coeff[iter]/F");
   tree->Branch("entries",&regEntries,"entries[iter]/F");

   for(uint32_t iR=0; iR<Type::nRegions; ++iR)  {
      std::vector<DetId> ids = Type::allDetIdsInRegion(iR);
      for(std::vector<DetId>::const_iterator iid = ids.begin();
                   			     iid != ids.end(); ++iid) {
         for(int k=0; k<NMaxIter; ++k) { regCoeff[k]=0.; regEntries[k] = 0.; } 
         EBDetId ebid(*iid);
         hashedIndex = ebid.hashedIndex();
         ieta = ebid.ieta();
         iphi = ebid.iphi();
         iSM = ebid.ism();
         iMod = ebid.im();
         iTT  = ebid.tower().hashedIndex();
         iTTeta = ebid.tower_ieta();
         iTTphi = ebid.tower_iphi();
         iter = parser->integerValue("iterations");
         if(iter==1) { // just dump the map. No calibration!
            //cout << "id: " << ebid << " coeff: " << calibMap->coeff(*iid) << endl;
            regCoeff[0] = calibMap->coeff(*iid);
            regEntries[0] = 0.;
            //cout << "id: " << ebid << " coeff: " << regCoeff[0] << endl;
         } else { 
            for(int j=0; j<iter; ++j) {
              regCoeff[j] = calibCoeff[iR][j];
              regEntries[j] = entries[iR][j];
            } // loop over iterations
         }
         tree->Fill();
      } // loop over DetId in regions
   } // loop over regions

   tree->Write();
   f->Close();

}

class PreshowerCluster {
 public: 
  
  PreshowerCluster(){
    x = 0.;
    y = 0.;
    z = 0.;
    energy = 0.;
    plane = 0.;
    goodcluster = false;
  };
  
  ~PreshowerCluster(){};

  double get_x(){return x;}
  double get_y(){return y;}
  double get_z(){return z;}
  double get_energy(){return energy;}
  int get_plane(){return plane;}
  bool get_goodcluster(){return goodcluster;}

  void set_x(double a) { x = a;}
  void set_y(double a) { y = a;}
  void set_z(double a) { z = a;}
  void set_energy(double a) { energy = a;}
  void set_plane(int a){plane = a;}
  void set_goodcluster(bool a){goodcluster = a;}

  

 private:

  double x;
  double y;
  double z;

  bool goodcluster;
  double energy;
  int plane;
  
  };

// nasty trick to keep the long function implementation elsewhere
#include "Analysis/Pi0Calib/src/EcalCalibAlgo.icc"

#endif
