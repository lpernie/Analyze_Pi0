#ifndef EcalPreshowerHardcodedGeometry_h
#define EcalPreshowerHardcodedGeometry_h

#include "Geometry/EcalCommonData/interface/EcalPreshowerNumberingScheme.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/EcalPreshowerGeometryRecord.h"
#include "CondFormats/AlignmentRecord/interface/ESAlignmentRcd.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "Geometry/Records/interface/PEcalPreshowerRcd.h"
#include <vector>
#include "Analysis/CalibTools/interface/ECALGeometry.h"

class EcalPreshowerHardcodedGeometry : public CaloSubdetectorGeometry
{
   public:

      typedef IdealGeometryRecord         IdealRecord   ;
      typedef EcalPreshowerGeometryRecord AlignedRecord ;
      typedef ESAlignmentRcd              AlignmentRecord ;
      typedef PEcalPreshowerRcd           PGeometryRecord ;

      typedef EcalPreshowerNumberingScheme NumberingScheme ;
      typedef CaloSubdetectorGeometry::ParVec ParVec ;
      typedef CaloSubdetectorGeometry::ParVecVec ParVecVec ;
      typedef ESDetId DetIdType ;

      enum { k_NumberOfCellsForCorners = ESDetId::kSizeForDenseIndexing } ;

      enum { k_NumberOfShapes = 2 } ;

      enum { k_NumberOfParametersPerShape = 3 } ;

      static std::string dbString() { return "PEcalPreshowerRcd" ; }

      virtual unsigned int numberOfShapes() const { return k_NumberOfShapes ; }
      virtual unsigned int numberOfParametersPerShape() const { return k_NumberOfParametersPerShape ; }

      EcalPreshowerHardcodedGeometry(ECALGeometry *geom) ;
  
      /// The EcalPreshowerGeometry will delete all its cell geometries at destruction time
      virtual ~EcalPreshowerHardcodedGeometry();

      void setzPlanes( float z1minus, 
		       float z2minus,
		       float z1plus, 
		       float z2plus ) ;

//Get the m_zplane content chosen by the user
double get_mzplane(int i){
if ((i >= 0)&&(i <= 3)){return m_zplane[i];}
   else {std::cout<<"Invalid m_zplane index"<<std::endl;
      return 0.;}
      }


      // Get closest cell
      virtual DetId getClosestCell( const GlobalPoint& r ) const ;

      // Get closest cell in arbitrary plane (1 or 2)
      DetId getClosestCellInPlane( const GlobalPoint& r     ,   int                plane   ) const ;

      virtual void initializeParms() ;
      virtual unsigned int numberOfTransformParms() const { return 3 ; }

      static std::string hitString() { return "EcalHitsES" ; }

      static std::string producerTag() { return "EcalPreshower" ; }

      static unsigned int numberOfAlignments() { return 8 ; }

      static unsigned int alignmentTransformIndexLocal( const DetId& id ) ;

      static unsigned int alignmentTransformIndexGlobal( const DetId& id ) ;

      static DetId detIdFromLocalAlignmentIndex( unsigned int iLoc ) ;

      static std::vector<HepGeom::Point3D<double> > localCorners( const double* pv,
						   unsigned int  i,
						   HepGeom::Point3D<double> &   ref ) ;

      static CaloCellGeometry* newCell( const GlobalPoint& f1 ,
					const GlobalPoint& f2 ,
					const GlobalPoint& f3 ,
					CaloCellGeometry::CornersMgr* mgr,
					const double*      parm ,
					const DetId&       detId   ) ;

   private:

      const double m_xWidWaf      ;
      const double m_xInterLadGap ;
      const double m_xIntraLadGap ;

      const double m_yWidAct      ;
      const double m_yCtrOff      ;

      double m_zplane[4];

      /// offline geometry
      ECALGeometry* offlineGeom;             // hardcoded geometry from root file

};


#endif

