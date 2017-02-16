#ifndef _TESSELATE_CURVE_H
#define _TESSELATE_CURVE_H

#include <vector>

#include "GoTools/geometry/ParamCurve.h"

namespace Go {

  std::vector<double> tesselate_curve(const Go::ParamCurve& pc,
				      unsigned int num_internal_points);
  
};

#endif
