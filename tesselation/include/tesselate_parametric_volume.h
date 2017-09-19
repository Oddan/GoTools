#ifndef _TESSELATE_PARAMETRIC_VOLUME_H
#define _TESSELATE_PARAMETRIC_VOLUME_H

#include <vector>
#include "GoTools/geometry/ParamCurve.h"

namespace TesselateUtils {

  // Estimate the length of a parametric curve, based on a specified number of samples
  double estimateCurveLength(const Go::ParamCurve& pc, const unsigned int num_samples);
  
  // tesselate parametric curve, returning a vector of parameter values
  std::vector<double> tesselateParametricCurve(const Go::ParamCurve& pc,
                                               const double vdist);
                                               
};


#endif
