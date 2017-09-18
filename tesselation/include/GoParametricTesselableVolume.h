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
  typedef std::tuple<EdgeCurve, uint, uint> EdgeType; // geometry and start/end vertex indices
  typedef std::pair<FaceSurface, vector<uint>> FaceType; // geometry and edge indices
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
