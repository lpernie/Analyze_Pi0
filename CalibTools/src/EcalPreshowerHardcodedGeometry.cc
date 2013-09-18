#include <iostream>
#include "TFile.h"
#include "Geometry/CaloGeometry/interface/PreshowerStrip.h"
#include "Geometry/CaloGeometry/interface/CaloGenericDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "Analysis/CalibTools/interface/EcalPreshowerHardcodedGeometry.h"

EcalPreshowerHardcodedGeometry::EcalPreshowerHardcodedGeometry(ECALGeometry *geom) :
   m_xWidWaf      ( 6.3  ) ,
   m_xInterLadGap ( 0.05 ) , // additional gap between wafers in adj ladders
   m_xIntraLadGap ( 0.04 ) , // gap between wafers in same ladder
   m_yWidAct      ( 6.1  ) ,
   m_yCtrOff      ( 0.05 )   // gap at center
{
   offlineGeom = geom;
   if (offlineGeom == 0) {std::cout<<"EcalPreshowerHardcodedGeometry was not linked properly from the Geometry file" << std::endl;
       std::exit(-1);}
}


EcalPreshowerHardcodedGeometry::~EcalPreshowerHardcodedGeometry() {
}

unsigned int
EcalPreshowerHardcodedGeometry::alignmentTransformIndexLocal( const DetId& id )
{
   const CaloGenericDetId gid ( id ) ;

   assert( gid.isES() ) ;

// plane 2 is split into 2 dees along x=0 for both signs of z

// plane 1 at zsign=-1 is split into 2 dees between six=19 and six=20 for siy<=20,
//                                             and six=21 and 22 for siy>=21

// plane 1 at zsign=+1 is split into 2 dees between six=20 and six=21 for siy<=20,
//                                             and six=19 and 20 for siy>=21


// Desired numbering 
//                LEFT    RIGHT (as one faces the Dee from the IP)
//  ES-  pl=2     0       1
//       pl=1     2       3    the reversal of pl=2 and pl=1 is intentional here (CM Kuo)
//  ES+  pl=1     4       5
//       pl=2     6       7

   const ESDetId esid ( id ) ;
   const int jx ( esid.six() - 1 ) ;
   const int jy ( esid.siy() - 1 ) ;
   const int jz ( esid.zside() + 1 ) ;
   const int pl ( esid.plane() - 1 ) ;
   const bool second ( 1 == pl ) ;
   const bool top   ( 19 < jy ) ;
   const bool negz  ( 0 == jz ) ;
   const int lrl    ( 20>jx ? 0 : 1 ) ;
   const int lrr    ( 22>jx ? 0 : 1 ) ;

   return ( second ? jx/20 + 3*jz :  // 2nd plane split along middle
	    ( negz && !top ? lrl + 2 :  // 1st plane at neg z and bottom half split at six=19&20
	      ( negz && top ? lrr + 2 : // 1st plane at neg z and top half split at six=21&22
		( !negz && !top ? lrr + 4 : lrl + 4 ) ) ) ) ; // opposite at positive z
}

DetId 
EcalPreshowerHardcodedGeometry::detIdFromLocalAlignmentIndex( unsigned int iLoc )
{
   return ESDetId( 1, 10 + 20*( iLoc%2 ), 10, 2>iLoc || 5<iLoc ? 2 : 1, 2*( iLoc/4 ) - 1 ) ;
}

unsigned int
EcalPreshowerHardcodedGeometry::alignmentTransformIndexGlobal( const DetId& id )
{
   return (unsigned int)DetId::Ecal - 1 ;
}


void 
EcalPreshowerHardcodedGeometry::initializeParms() 
{
   typedef CaloSubdetectorGeometry::CellCont Cont ;
   unsigned int n1minus ( 0 ) ;
   unsigned int n2minus ( 0 ) ;
   unsigned int n1plus ( 0 ) ;
   unsigned int n2plus ( 0 ) ;
   double z1minus ( 0 ) ;
   double z2minus ( 0 ) ;
   double z1plus ( 0 ) ;
   double z2plus ( 0 ) ;
   //const Cont& con ( cellGeometries() ) ;
   //for( unsigned int i ( 0 ) ; i != offlineGeom->getPositionMap().size() ; ++i )

   std::map<DetId,GlobalPoint>::const_iterator it;
   for( it = (offlineGeom->getPositionMap()).begin(); it != (offlineGeom->getPositionMap()).end(); ++it)
   {
      
      /// selecting only preshower detIds
      if( it->first.subdetId() != EcalPreshower) continue; 

      const ESDetId esid ( it->first ) ;

      if( 1 == esid.plane() )
      {
	 if( 0 > esid.zside() )
	 {
	    z1minus += it->second.z() ;
	    ++n1minus ;
	 }
	 else
	 {
	    z1plus += it->second.z() ;
	    ++n1plus ;
	 }
      }
      if( 2 == esid.plane() )
      {
	 if( 0 > esid.zside() )
	 {
	    z2minus += it->second.z() ;
	    ++n2minus ;
	 }
	 else
	 {
	    z2plus += it->second.z() ;
	    ++n2plus ;
	 }
      }
//      if( 0 == z1 && 1 == esid.plane() ) z1 = fabs( i->second->getPosition().z() ) ;
//      if( 0 == z2 && 2 == esid.plane() ) z2 = fabs( i->second->getPosition().z() ) ;
//      if( 0 != z1 && 0 != z2 ) break ;
   }

   //std::cout << "initializeParms::total ES detid = " << n2plus+n2minus+n1plus+n1minus << std::endl;
   assert( 0 != n1minus &&
	   0 != n2minus &&
	   0 != n1plus  &&
	   0 != n2plus     ) ;
   z1minus /= (1.*n1minus) ;
   z2minus /= (1.*n2minus) ;
   z1plus /= (1.*n1plus) ;
   z2plus /= (1.*n2plus) ;
   assert( 0 != z1minus &&
	   0 != z2minus &&
	   0 != z1plus  &&
	   0 != z2plus     ) ;
   setzPlanes( z1minus, z2minus, z1plus, z2plus ) ;
}


void 
EcalPreshowerHardcodedGeometry::setzPlanes( float z1minus, 
				   float z2minus,
				   float z1plus, 
				   float z2plus ) 
{
   assert( 0 > z1minus &&
	   0 > z2minus &&
	   0 < z1plus  &&
	   0 < z2plus     ) ;

   m_zplane[0] = z1minus ;
   m_zplane[1] = z2minus ;
   m_zplane[2] = z1plus ;
   m_zplane[3] = z2plus ;
}

// Get closest cell, etc...
DetId 
EcalPreshowerHardcodedGeometry::getClosestCell( const GlobalPoint& point ) const
{
  return getClosestCellInPlane( point, 2 );
} 


DetId 
EcalPreshowerHardcodedGeometry::getClosestCellInPlane( const GlobalPoint& point,     int                plane          ) const
{
   const double x ( point.x() ) ;
   const double y ( point.y() ) ;
   const double z ( point.z() ) ;

   //std::cout << "EcalPreshowerGeometry::reference position:: "<< x << " " << y << " " << z << std::endl;

   if( 0 == z    ||
       1 > plane ||
       2 < plane    ) return DetId( 0 ) ;

   const unsigned int iz ( ( 0>z ? 0 : 2 ) + plane - 1 ) ;
   //std::cout<<"iz : "<< iz <<std::endl;

   const double ze ( m_zplane[iz]) ;
   const double xe ( x * ze/z ) ;
   const double ye ( y * ze/z ) ;

   //std::cout << "EcalPreshowerGeometry::m_zplane[" << iz << "] = " << m_zplane[iz] << "  (xe,ye,ze) = ("<< xe << "," << ye << "," << ze << ")" << std::endl;

   const double x0 ( 1 == plane ? xe : ye ) ;
   const double y0 ( 1 == plane ? ye : xe ) ;

   static const double xWid ( m_xWidWaf + m_xIntraLadGap + m_xInterLadGap ) ;

   const int row ( 1 + int( y0 + 20.*m_yWidAct - m_yCtrOff )/m_yWidAct ) ;
   const int col ( 1 + int( ( x0 + 20.*xWid )/xWid ) ) ;

   double closest ( 1e9 ) ; 

   DetId detId ( 0 ) ;

   const int jz ( 0 > ze ? -1 : 1 ) ;


   //std::cout<<"** p="<<point<<", ("<<xe<<", "<<ye<<", "<<ze<<"), row="<<row<<", col="<<col<<std::endl;

   for( int ix ( -1 ); ix != 2 ; ++ix ) // search within +-1 in row and col
   {
      for( int iy ( -1 ); iy != 2 ; ++iy )
      {
	 for( int jstrip ( ESDetId::ISTRIP_MIN ) ; jstrip <= ESDetId::ISTRIP_MAX ; ++jstrip )
	 {
	    const int jx ( 1 == plane ? col + ix : row + iy ) ;
	    const int jy ( 1 == plane ? row + iy : col + ix ) ;
	    if( ESDetId::validDetId( jstrip, jx, jy, plane, jz ) )
	    {
	       const ESDetId esId ( jstrip, jx, jy, plane, jz ) ;
	       //const unsigned int index ( esId.denseIndex() ) ;
	       //const GlobalPoint& p ( cellGeometries()[ index ]->getPosition() ) ;
	       //std::cout << "getClosestCell:: looking for position of detid " << esId.rawId() << std::endl;
	       const GlobalPoint& p ( offlineGeom->getPosition(esId) ) ; // hic
             //std::cout << "getClosestCell::point being considered (" << p.x() << "," <<p.y() << ","<<p.z() << ")" << std::endl; 
	       const double dist2 ( (p.x()-xe)*(p.x()-xe) + (p.y()-ye)*(p.y()-ye) ) ;
	       if( dist2 < closest )
	       {
		  closest = dist2 ;
		  detId   = esId  ;
	       }
	    }
         //else
         //{
         //   //std::cout << "-- EcalPreshowerGeometry::nearby detId is not valid" << std::endl;
         //}
	 }
      }
   }
   //std::cout << "EcalPreshowerGeometry::distance " << closest << std::endl;

   return detId ;
}

std::vector<HepGeom::Point3D<double> > 
EcalPreshowerHardcodedGeometry::localCorners( const double* pv,
				     unsigned int  i,
				     HepGeom::Point3D<double> &   ref )
{
   return ( PreshowerStrip::localCorners( pv, ref ) ) ;
}

CaloCellGeometry* 
EcalPreshowerHardcodedGeometry::newCell( const GlobalPoint& f1 ,
				const GlobalPoint& f2 ,
				const GlobalPoint& f3 ,
				CaloCellGeometry::CornersMgr* mgr,
				const double*      parm ,
				const DetId&       detId    ) 
{
   return ( new PreshowerStrip( f1, mgr , parm ) ) ;
}
