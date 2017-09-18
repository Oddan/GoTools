#ifndef _GOPARAMETRICTESSELABLEVOLUME_H
#define _GOPARAMETRICTESSELABLEVOLUME_H

#include <tuple>
#include "TesselableVolume.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/utils/Point.h"
#include "GoTools/trivariatemodel/ftVolume.h"

namespace TesselateUtils
{
  
struct GoParametricSpaceTraits {
  typedef std::shared_ptr<const Go::ParamCurve>   EdgeCurve;
  typedef std::shared_ptr<const Go::ParamSurface> FaceSurface;

  typedef Go::Point PointType;

  // geometry and start/end vertex indices
  typedef std::tuple<EdgeCurve, uint, uint> EdgeType; 

  // geometry, face edge index, and forward/reverse flag
  typedef struct {
    FaceSurface surf;
    vector<pair<uint, bool>> ix; // boolean indicate orientation of edge.  If
                                 // 'false', edge is oriented in the reverse
                                 // direction
  } FaceType;

  //typedef std::pair<FaceSurface, vector<std::pair<uint, bool>>> FaceType; 
  typedef std::shared_ptr<const Go::ParamVolume>  VolumeType;
};

using GoParametricTesselableVolume = TesselableVolume<GoParametricSpaceTraits>;  

// ----------------------------------------------------------------------------
// Additional constructor, making a TesselableVolume out of a Go::ftVolume
template<> template<>
GoParametricTesselableVolume::TesselableVolume(Go::ftVolume& fvol);
// ----------------------------------------------------------------------------

};

#endif
