#include <iostream> // for debugging
#include <cmath>
#include <algorithm>
#include "tesselate_utils.h"
#include "tesselate_polyhedron.h"
#include "polyhedral_energies.h"
#include "GoTools/utils/GeneralFunctionMinimizer.h"

using namespace std;
using namespace TesselateUtils; 

namespace {

const double PI = 3.14159265358979323;

class PolygonEnergyFunctor
{
public:
  PolygonEnergyFunctor(const Point2D* const polygon,
		       const unsigned int num_corners,
		       const unsigned int num_ipoints,
		       const double radius);

  double operator()(const double* const arg) const;
  void grad(const double* arg, double* grad) const;
  double minPar(int n) const;
  double maxPar(int n) const;
private:

  bool use_cached(const double* const arg) const;
  void update_cache(const double* const arg) const;
  const Point2D* const poly_;
  const unsigned int nc_; // num polygon corners
  const unsigned int ni_; // num interior points
  const double r_; // radius
  const std::array<double, 4> bbox_; // bounding box

  mutable ValAndDer cached_result_;
  mutable vector<double> cached_arg_;

};
  
// ----------------------------------------------------------------------------  
vector<Point2D> init_startpoints(const Point2D* const polygon,
				 const unsigned int num_corners,
				 const double vdist);
// choose some initial startpoints as a basis for subsequent optimization  
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void optimize_interior_points(const Point2D* const polygon,
			      const unsigned int num_corners,
			      Point2D* const ipoints,
			      const unsigned int num_ipoints,
			      const double vdist);
// ----------------------------------------------------------------------------  
  
}; // end anonymous namespace

namespace TesselateUtils {

// ============================================================================
Mesh2D tesselatePolygon2D(const Point2D* const polygon,
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
  optimize_interior_points(&bpoints[0], (uint)bpoints.size(),
			   &ipoints[0], (uint)ipoints.size(),
			   vdist*1.5);


  
  // @@computing and reporting internal energy
  const double R = 2 * vdist;
  const auto e = polygon_energy(&bpoints[0], (unsigned int)bpoints.size(),
				&ipoints[0], (unsigned int)ipoints.size(), 
				R);

  cout << "Energy is:" << e.val << endl;

  vector<Point2D> points(bpoints);
  points.insert(points.end(), ipoints.begin(), ipoints.end());

  const auto tris = triangulate_domain(&points[0],
				       (uint)bpoints.size(),
				       (uint)points.size(),
				       2*vdist);
  return {points, tris};
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

namespace {

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
  
// ----------------------------------------------------------------------------
void optimize_interior_points(const Point2D* const polygon,
			      const unsigned int num_corners,
			      Point2D* const ipoints,
			      const unsigned int num_ipoints,
			      const double vdist)
// ----------------------------------------------------------------------------  
{
  // Setting up function to minimize (energy function);
  auto efun = PolygonEnergyFunctor(polygon, num_corners, num_ipoints, vdist);
  
  // Setting up function minimizer
  Go::FunctionMinimizer<PolygonEnergyFunctor>
    funcmin(num_ipoints * 2, efun, (double* const)&ipoints[0], 1e-1); // @@ TOLERANCE?
  
  // do the minimization
  Go::minimise_conjugated_gradient(funcmin);

  // copying results back.  The cast is based on the knowledge that Points2D
  // consists of POD and that points are stored contiguously in memory
  // (as x1,y1,x2,y2, ...)
  double* const target = (double* const) ipoints;
  std::copy(funcmin.getPar(), funcmin.getPar() + funcmin.numPars(), target);
}

// ----------------------------------------------------------------------------    
PolygonEnergyFunctor::PolygonEnergyFunctor(const Point2D* const polygon,
					   const unsigned int num_corners,
					   const unsigned int num_ipoints,
					   const double radius)
// ----------------------------------------------------------------------------    
  : poly_(polygon), nc_(num_corners), ni_(num_ipoints), r_(radius),
    bbox_(bounding_box(polygon, num_corners))
{}

// ----------------------------------------------------------------------------    
void PolygonEnergyFunctor::update_cache(const double* const arg) const
// ----------------------------------------------------------------------------    
{
  if (cached_arg_.empty())
    cached_arg_.resize(2*ni_);

  copy(arg, arg + 2 * ni_, &cached_arg_[0]);
  cached_result_ = polygon_energy(poly_, nc_,
				  (const Point2D* const) &cached_arg_[0],
				  ni_, r_);
}

// ----------------------------------------------------------------------------    
double PolygonEnergyFunctor::operator()(const double* const arg) const
// ----------------------------------------------------------------------------    
{
  // Very inefficient to have separate methods for function value and function
  // gradient.  We have therefore implemented a caching system here.
  if (!use_cached(arg)) {
    update_cache(arg);
    //cout << "updating fval" << endl;
  } else {
    //cout << "using cached fval" << endl;
  }
  return cached_result_.val;
}

// ----------------------------------------------------------------------------    
void PolygonEnergyFunctor::grad(const double* arg, double* grad) const
// ----------------------------------------------------------------------------
{
  // Very inefficient to have separate methods for function value and function
  // gradient.  We have therefore implemented a caching system here.
  if (!use_cached(arg)) {
    update_cache(arg);
    //cout << "updating derivative" << endl;
  } else {
    //cout << "using cached derivative" << endl;
  }
  const double* const dp = (const double* const)&cached_result_.der[0];
  copy(dp, dp + 2 * ni_, grad);
}

// ----------------------------------------------------------------------------    
double PolygonEnergyFunctor::minPar(int n) const
// ----------------------------------------------------------------------------
{
  return bbox_[(n % 2) * 2];
}

// ----------------------------------------------------------------------------
double PolygonEnergyFunctor::maxPar(int n) const 
// ----------------------------------------------------------------------------  
{
  return bbox_[(n % 2) * 2 + 1];
}

// ----------------------------------------------------------------------------  
bool PolygonEnergyFunctor::use_cached(const double* const arg) const
// ----------------------------------------------------------------------------
{
  return (cached_arg_.size() > 0) && (std::equal(arg, arg + 2 * ni_, &cached_arg_[0]));
}

}; // end anonymous namespace
