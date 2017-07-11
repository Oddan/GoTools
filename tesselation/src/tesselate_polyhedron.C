#include <iostream> // for debugging
#include <cmath>
#include <algorithm>
#include "tesselate_utils.h"
#include "tesselate_polyhedron.h"
#include "polyhedral_energies.h"

using namespace std;
using namespace TesselateUtils; 

namespace {

const double PI = 3.14159265358979323;
  
// ----------------------------------------------------------------------------  
vector<Point2D> init_startpoints(const Point2D* const polygon,
				 const unsigned int num_corners,
				 const double vdist)
// choose some initial startpoints as a basis for subsequent optimization  
// ----------------------------------------------------------------------------
{
  //const double poly_area = polygon_area(polygon, num_corners);
  const array<double, 4> bbox = bounding_box(polygon, num_corners);
  const double bbox_lx = bbox[1] - bbox[0];
  const double bbox_ly = bbox[3] - bbox[2];
  const double bbox_area =  bbox_lx * bbox_ly;

  // computing amount of points needed to approximately cover the bounding box, where points
  // are approximately equidistant by 'vdist'
  const int N_bbox = (int)ceil(bbox_area * 3 / (PI * vdist * vdist));
  const int nx = (int)ceil(sqrt(N_bbox * bbox_lx / bbox_ly));
  const int ny = (int)ceil(N_bbox/nx);

  // constructing regular grid of points within the bounding box
  const vector<Point2D> gridpoints = generate_grid_2D(Point2D {bbox[0], bbox[2]},
  						      Point2D {bbox[1], bbox[3]},
  						      nx, ny);
  
// keeping the points that fall within the original polygon
  const double TOL_FAC = 0.5 * vdist; // do not keep points too close to the boundary
  vector<Point2D> result = inpolygon(&gridpoints[0], (unsigned int)gridpoints.size(),
				     polygon,  num_corners, TOL_FAC);

  // add perturbation to avoid 'symmetry locking' during the optimization stage
  const double PM = 1.0e-1 * min(bbox_lx/nx, bbox_ly/ny); // perturbation magnitude
  transform(result.begin(), result.end(), result.begin(), 
	    [PM] (Point2D p) { return Point2D {p[0] + random_uniform(-PM, PM),
		                               p[1] + random_uniform(-PM, PM)};});

  return result;
}
  
}; // end anonymous namespace

namespace TesselateUtils {

// ============================================================================
vector<Point2D> tesselatePolygon2D(const Point2D* const polygon,
                                  const unsigned int num_corners,
                                  const double vdist)
// ============================================================================
{
  // computing the points defining the boundary polygon 

  vector<Point2D> bpoints; // "boundary points"
  for (unsigned int i = 0; i != num_corners; ++i) {
    const auto tseg = tesselateSegment2D(polygon[i], 
					 polygon[(i+1) % num_corners], 
					 vdist);
    bpoints.insert(bpoints.end(), tseg.begin(), tseg.end()-1);
  }

  // Choosing a set of interior points which we can later move around to optimal
  // locations 
  vector<Point2D> ipoints = init_startpoints(polygon, num_corners, vdist); 
  
  // Optimizing position of interior points

  // @@computing and reporting internal energy
  const double R = 2 * vdist;
  const auto e = polygon_energy(&bpoints[0], (unsigned int)bpoints.size(),
				&ipoints[0], (unsigned int)ipoints.size(), 
				R);

  cout << "Energy is:" << e.val << endl;

  vector<Point2D> points(bpoints);
  points.insert(points.end(), ipoints.begin(), ipoints.end());
  return points;
}

// ============================================================================
vector<Point2D> tesselateSegment2D(const Point2D& p1,
				   const Point2D& p2,
				   const double vdist)
// ============================================================================
{
  const double seg_len = dist(p1, p2);
  const unsigned int num_intervals = (unsigned int)ceil(seg_len/vdist);

  return interpolate(p1, p2, num_intervals - 1);
}
  
};
