#include <cmath>
#include "ParametricObjectEnergyFunctor.h"
#include "GoTools/utils/GeneralFunctionMinimizer.h"
#include "tesselate_parametric_volume.h"

using namespace TesselateUtils;
using namespace std;
using namespace Go;

namespace {

typedef shared_ptr<const ParamSurface> SurfPtr;

vector<Point3D> init_startpoints(const SurfPtr sp,
                                 const Point2D* const bpoints,
                                 const uint num_bpoints,
                                 const double vdist);
                                 
}; // end anonymous namespace

namespace TesselateUtils {

// ----------------------------------------------------------------------------
std::vector<Point2D>
tesselateParametricSurface(const SurfPtr surf,
                           const Point2D* const bpoints,
                           const uint num_bpoints,
                           const double vdist)
// ----------------------------------------------------------------------------
{
  vector<Point2D> ipoints = init_startpoints(ps, bpoints, num_bpoints, vdist);
}

  
}; // end namespace TesselateUtils

namespace {

// ----------------------------------------------------------------------------
vector<Point3D> init_startpoints(const SurfPtr sp,
                                 const Point2D* const bpoints,
                                 const uint num_bpoints,
                                 const double vdist)
// ----------------------------------------------------------------------------
{
  const double poly_area = polygon_area(polygon, num_corners);

  // determine area of _unbounded_ surface (we need the area corresponding to
  // the use of the full parameter domain)
  
}

  
};
