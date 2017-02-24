#include <assert.h>
#include "tesselate_surface.h"
#include "tritools.h"

using namespace Go;
using namespace std;
using namespace TriTools;

namespace {
// ----------------------------------------------------------------------------
vector<double> define_parvec(double minparam, double maxparam, unsigned int num_intparams) {
// ----------------------------------------------------------------------------
  vector<double> result(num_intparams+2);
  result.front() = minparam;
  result.back()  = maxparam;
  for (size_t i = 1; i != num_intparams+1; ++i)
    result[i] = (maxparam-minparam)/(num_intparams + 1) * double(i);
  return result;
}

// ----------------------------------------------------------------------------
vector<Point> setup_parpoint_grid(const vector<double>& u, const vector<double>& v)
// ----------------------------------------------------------------------------
{
  vector<Point> result;
  for (vpar : v)
    for (upar : u)
      result.push_back({upar, vpar});
  return result;
}
  
}; 


namespace Go {

  // ----------------------------------------------------------------------------
SurfaceTriangulation tesselate_surface(const ParamSurface& ps,
				       const unsigned int num_internal_points_u,
				       const unsigned int num_internal_points_v,
				       const double margin_boundary,
				       const vector<vector<Point>>& bounding_curves)
// ----------------------------------------------------------------------------
{
  // Regularly sample surface
  const RectDomain dom = ps.containingDomain();
  const vector<double> u = define_parvec(dom.umin(), dom.umax(), num_internal_points_u);
  const vector<double> v = define_parvec(dom.vmin(), dom.vmax(), num_internal_points_v);

  vector<Point> parpoints = setup_parpoint_grid(u, v);
  
  // Clip surface against eventual bounding curves
  if (bounding_curves.size() > 0) {
    vector<int> keep = points_inside_loops(bounding_curves, parpoints, margin_boundary);
    vector<Point> parpoints_new;
    for (size_t i = 0; i != parpoints.size(); ++i)
      if (keep[i])
	parpoints_new.push_back(parpoints[i]);
    parpoints.swap(parpoints_new);
  }
  
  // Create triangulation
  return triangulate_with_boundaries(parpoints, bounding_curves);
}


  
}; // end namespace Go;
