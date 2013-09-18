#ifndef PropagationTools_HH
#define PropagationTools_HH

#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/Point3D.h"


class PropagationTools {

public:
  PropagationTools(float radius);
  bool getTrackImpactPositionEB(const reco::Track* tk_ref,
                              const TrackerGeometry* trackerGeom, const MagneticField* magField,
                              math::XYZPoint& ew);

  bool getTrackImpactPositionEE(const reco::Track* tk_ref,
                              const TrackerGeometry* trackerGeom, const MagneticField* magField,
                              math::XYZPoint& ew);


private:

   ReferenceCountingPointer<BoundCylinder>  theBarrel_;
   ReferenceCountingPointer<BoundDisk>  EEpos_;
   ReferenceCountingPointer<BoundDisk>  EEneg_;


};

#endif
