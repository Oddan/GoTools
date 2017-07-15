#ifndef _TRIANGULATE_DOMAIN_H
#define _TRIANGULATE_DOMAIN_H

#include <array>
#include <vector>
#include "common_defs.h"

namespace TesselateUtils {

  using Triangle = std::array<uint, 3>;

  // Compute a generalized Delaunay triangulation of a not necessarily convex
  // domain with a set of given interior points.  The boundary points and
  // interior points are provided by 'points'.  The boundary points are listed
  // first, and there are 'num_bpoints' of them.  The total number of points
  // (boundary points and interior points) is given by 'tot_num_points'.  It is
  // assumed that all provided interior points lie within the boundary, and that
  // the spacing of the points is such that no triangle edge needs to be longer
  // than 'vdist'.  (If you are unsure, you can set 'vdist' as large as you
  // want, but runtime performance is improved by keeping 'vdist' low.  Boundary
  // points should be provided in clockwise order.  Each triangle is returned as
  // a set of three indices in to the list of points.  These indices are
  // provided in counterclockwise order for each triangle.
  std::vector<Triangle> triangulate_domain(const Point2D* const points,
					   const uint num_bpoints, 
					   const uint tot_num_points,
					   const double vdist);
};

#endif
