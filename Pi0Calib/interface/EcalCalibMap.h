#ifndef EcalCalibMap_H
#define EcalCalibMap_H

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

//EcalCalibType :: { Xtal, EtaRing, TrigTower };

template<typename Type>
class EcalCalibMap {
  public:
    EcalCalibMap() {
     for(uint32_t i=0; i<Type::nRegions; ++i) { map[i] = 1.; }
    };

    static float& coeff(const DetId& id){ return map[Type::iRegion(id)]; }
    static uint32_t nCoeff() { return nRegions; }

    float& operator[](const DetId& id) { return map[Type::iRegion(id)]; }
    const float& operator[](const DetId& id) const { return map[Type::iRegion(id)]; }


  private:
    static float map[Type::nRegions];
    static const uint32_t nRegions = Type::nRegions;

};
#endif
