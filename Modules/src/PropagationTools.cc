#include "Analysis/Modules/interface/PropagationTools.h"



#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"


PropagationTools::PropagationTools(float radius) {

  const float epsilon = 0.001;
  Surface::RotationType rot; // unit rotation matrix
  const float barrelRadius = radius;
  const float barrelHalfLength = 270.9f;
   theBarrel_ = ReferenceCountingPointer<BoundCylinder>( 
           new BoundCylinder( Surface::PositionType(0,0,0), rot,
           SimpleCylinderBounds( barrelRadius-epsilon, barrelRadius+epsilon, 
           -barrelHalfLength, barrelHalfLength )  )  );


   const float endcapRadMin = 21.0f;
   const float endcapRadMax = 181.1f;
   const float endcapZ = 304.5f;

   EEpos_ = ReferenceCountingPointer<BoundDisk>( 
           new BoundDisk( Surface::PositionType(0,0,endcapZ), rot,
           SimpleDiskBounds( endcapRadMin, endcapRadMax, -epsilon , epsilon) ) );

   EEneg_ = ReferenceCountingPointer<BoundDisk>( 
           new BoundDisk( Surface::PositionType(0,0,-endcapZ), rot,
           SimpleDiskBounds( endcapRadMin, endcapRadMax, -epsilon , epsilon) ) );


}

bool PropagationTools::getTrackImpactPositionEB(const reco::Track* trk,
                                          const TrackerGeometry* trackerGeom, 
                                          const MagneticField* magField,
                                          math::XYZPoint& ew)
{
  PropagatorWithMaterial propag( alongMomentum, 0.105, magField );
  TrajectoryStateTransform transformer;

  //const TrajectoryStateOnSurface myTSOS = transformer.innerStateOnSurface(*(*ref), *trackerGeom, magField);
  const TrajectoryStateOnSurface myTSOS = transformer.outerStateOnSurface(*trk, *trackerGeom, magField);
  TrajectoryStateOnSurface  stateAtECAL;
  stateAtECAL = propag.propagate(myTSOS, *theBarrel_);
  if (stateAtECAL.isValid()){
    ew = stateAtECAL.globalPosition();
    return true;
  }
  else
    return false;
}

bool PropagationTools::getTrackImpactPositionEE(const reco::Track* trk,
                                          const TrackerGeometry* trackerGeom, 
                                          const MagneticField* magField,
                                          math::XYZPoint& ew)
{
  PropagatorWithMaterial propag( alongMomentum, 0.105, magField );
  TrajectoryStateTransform transformer;

  const TrajectoryStateOnSurface myTSOS = transformer.outerStateOnSurface(*trk, *trackerGeom, magField);
  TrajectoryStateOnSurface  stateAtECAL;
  if( trk->eta()>0.) stateAtECAL = propag.propagate(myTSOS, *EEpos_);
  else               stateAtECAL = propag.propagate(myTSOS, *EEneg_);
  if (stateAtECAL.isValid()){
    ew = stateAtECAL.globalPosition();
    return true;
  }
  else
    return false;
}

