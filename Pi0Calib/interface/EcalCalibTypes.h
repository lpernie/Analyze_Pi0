#ifndef EcalCalibTypes_H
#define EcalCalibTypes_H

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include<vector>
//#include<iostream>
//using namespace std;

namespace EcalCalibType {

  class Xtal {
    public:
      typedef EBDetId ID;
      static const uint32_t nRegions = EBDetId::MAX_SM*EBDetId::kCrystalsInPhi*EBDetId::kCrystalsInEta;
      //                    indices from 0 -> 61200-1
      static uint32_t iRegion(const DetId& id) { return EBDetId(id).hashedIndex(); }
      static ID detIdFromRegion(uint32_t iR) {
         return EBDetId::detIdFromDenseIndex( iR );
      }

      static std::vector<DetId> allDetIdsInRegion(uint32_t iR) {
         std::vector<DetId> ids;
         ids.push_back( detIdFromRegion(iR) );
         return ids;
      }

    private:
  };

  class EtaRing {
    public:
      typedef EBDetId ID;
      static const uint32_t nRegions = 2*EBDetId::kCrystalsInEta;
      static uint32_t iRegion(const DetId& id) { 
        EBDetId ebid(id); 
        int ieta = ebid.ieta()+EBDetId::kCrystalsInEta;
        //                   EB- : 0->84   EB+ : 85 -> 169 
        return (ebid.zside() < 0) ? ieta : ieta-1 ;
      }

      static ID detIdFromRegion(uint32_t iR) {
         int ieta = iR - EBDetId::kCrystalsInEta;
         ieta = (iR<uint32_t(EBDetId::kCrystalsInEta)) ?  ieta : ieta+1;
         return EBDetId( ieta, 1 );
      }

      static std::vector<DetId> allDetIdsInRegion(uint32_t iR) {
         int ieta = int(iR) - EBDetId::kCrystalsInEta;
         ieta = (iR<uint32_t(EBDetId::kCrystalsInEta)) ?  ieta : ieta+1;
         std::vector<DetId> ids;
         ids.reserve(EBDetId::MAX_IPHI);
         for(int iphi=1; iphi<=EBDetId::MAX_IPHI; ++iphi) 
            ids.push_back( EBDetId( ieta, iphi ) );
         //cout << "allDetIdsInRegion == iR " << iR << " xtals: " 
         //     << ids.size() << endl;
         return ids;
      }

    private:
  };

  class TrigTower {
    public:
      typedef EcalTrigTowerDetId ID;
      static const uint32_t nRegions = EcalTrigTowerDetId::kEBTowersPerSM*EBDetId::MAX_SM;
      static uint32_t iRegion(const DetId& id) { 
        EBDetId ebid(id); 
        //       0 -> 2448-1
        return (ebid.tower().hashedIndex());
      }
      static ID detIdFromRegion(uint32_t iR) {
         return EcalTrigTowerDetId::detIdFromDenseIndex( iR ); // was +1 but it's a mistake! 6/6/10
      }

      static std::vector<DetId> allDetIdsInRegion(uint32_t iR) {
         EcalTrigTowerDetId iTT = EcalTrigTowerDetId::detIdFromDenseIndex( iR );
         std::vector<DetId> ids;
         //cout << "iR: " << iR  << " iTT: " << iTT << endl;
         int ieta;
         for(int i=1; i <=5; ++i) {
            ieta = ( (iTT.ietaAbs()-1)*5 + i ) * iTT.zside();
            int phixtalMin=((iTT.iphi()-1)*5+11)%360;
            int phixtalMax=((iTT.iphi())*5+10)%360;
            if(phixtalMax<=0) phixtalMax+=360;
            for(int iphi=phixtalMin; iphi<=phixtalMax; ++iphi) {
               //cout << "---- ieta: " << ieta << "  iphi: " << iphi  << endl;
               ids.push_back( EBDetId( ieta,iphi, EBDetId::ETAPHIMODE ) );
            }
         }
         return ids;
      }

    private:
  };

} // end of namespace

#endif
