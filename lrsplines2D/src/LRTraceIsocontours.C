#include <functional>
#include "GoTools/lrsplines2D/LRTraceIsocontours.h"
#include "GoTools/lrsplines2D/SSurfTraceIsocontours.h"

using namespace std;
using namespace Go;

using LRSurfPtr = shared_ptr<LRSplineSurface>;

namespace {

// ----------------------------------------------------------------------------
// Wrapper to simplify syntax when the shared pointer itself is not required
// except to manage the allocated memory
inline shared_ptr<SplineSurface> as_spline_surf(const shared_ptr<LRSplineSurface> lrs)
// ----------------------------------------------------------------------------
{
  return shared_ptr<SplineSurface>(lrs->asSplineSurface());
}
  
}// end anonymous namespace


namespace Go
{

// ============================================================================
vector<CurveVec> LRTraceIsocontours(const LRSplineSurface& lrs,
				    const std::vector<double>& isovals,
				    const double tol,
				    const bool include_3D_curves,
				    const bool use_sisl_marching)
// ============================================================================
{
  //break LR surface up into individual patches
  const vector<LRSurfPtr> surf_fragments = lrs.subdivideIntoSimpler();

  // Compute isocontours for each surface fragment
  const function<vector<CurveVec>(LRSurfPtr)> compute_isovals {
    [&] (LRSurfPtr l){return SSurfTraceIsocontours(*as_spline_surf(l), isovals,
						   tol, include_3D_curves,
						   use_sisl_marching);}};

  // computing isocurves for each surface fragment (vector<vector<CurveVec>>)
  const auto curve_fragments = apply_transform(surf_fragments, compute_isovals);
  
  // merge isocontours across patches
  const vector<CurveVec> result ();// = merge_isocontours(curve_fragments, surf_fragments);

  return vector<CurveVec> ();// @@ DUMMY
}


  
}; // end namespace Go;


