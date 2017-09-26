#ifndef _TESSELATE_PARAMETRIC_VOLUME_H
#define _TESSELATE_PARAMETRIC_VOLUME_H

#include <vector>
#include "common_defs.h"
#include "GoTools/geometry/ParamCurve.h"

namespace TesselateUtils {

  // tesselate parametric curve, returning a vector of parameter values
  std::vector<double>
  tesselateParametricCurve(const std::shared_ptr<const Go::ParamCurve> pc,
                           const double vdist);

  std::vector<Point2D>
  tesselateParametricSurface(const std::shared_ptr<const Go::ParamSurface> ps,
                             const double vdist);
  
                                               
};

#endif
