#ifndef _PARAMETRIC_OBJECT_ENERGIES_H
#define _PARAMETRIC_OBJECT_ENERGIES_H

#include <vector>
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "common_defs.h"

namespace TesselateUtils
{

// ============================================================================
ValAndDer<Point1D> parametric_curve_energy(const shared_ptr<const Go::ParamCurve> curve,
                                           const double startpar,
                                           const double endpar,
                                           const double* const ipar,
                                           const unsigned int num_ipar,
                                           const double vdist);
// ============================================================================

// ============================================================================

// ============================================================================  
  
}; // end namespace TesselateUtils

#endif
