#include <iostream> // for debugging
#include <fstream> // for debugging
#include <iterator> // for debugging
#include <assert.h>
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

// Class needed by the "optimize_interior_point" function for the call to
// Go::minimise_conjugated_gradient.
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

  mutable ValAndDer<Point2D> cached_result_;
  mutable vector<double> cached_arg_;

};

// Class needed by the "optimize_interior_point" function for the call to
// Go::minimise_conjugated_gradient.
class PolyhedronEnergyFunctor
{
public:
  PolyhedronEnergyFunctor(const Point3D* const bpoints,
                          const unsigned int num_bpoints,
                          const Triangle* const btris,
                          const unsigned int num_btris,
                          const Point3D* const ipoints,
                          const unsigned int num_ipoints,
                          const double radius);

  double operator()(const double* const arg) const;
  void grad(const double* arg, double* grad) const;
  double minPar(int n) const;
  double maxPar(int n) const;
private:

  bool use_cached(const double* const arg) const;
  void update_cache(const double* const arg) const;
  const Point3D* const bpoints_;
  const unsigned int nb_; // num boundary points
  const Triangle* const btris_;
  const unsigned int nt_; // num triangles
  const unsigned int ni_; // num interior points
  const double r_; // radius
  const std::array<double, 6> bbox_; // bounding box

  mutable ValAndDer<Point3D> cached_result_;
  mutable vector<double> cached_arg_;
};

// ----------------------------------------------------------------------------  
vector<Point2D> init_startpoints(const Point2D* const polygon,
				 const unsigned int num_corners,
				 const double vdist);
// choose some initial startpoints as a basis for subsequent optimization  
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
vector<Point3D> init_startpoints(const Point3D* const bpoints,
                                 const unsigned int num_bpoints,
                                 const Triangle* btris,
                                 const unsigned int num_btris,
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


// ----------------------------------------------------------------------------    
void optimize_interior_points(const Point3D* bpoints,
                              const unsigned int num_bpoints,
                              const Triangle* const btris,
                              const unsigned int num_btris,
                              Point3D* const ipoints,
                              const unsigned int num_ipoints,
                              const double vdist);
// ----------------------------------------------------------------------------    
  
}; // end anonymous namespace

namespace TesselateUtils {

// ============================================================================
Mesh2D tesselatePolygon2D(const Point2D* const polygon,
			  const unsigned int num_corners,
			  const double vdist,
			  const bool tesselate_boundary)
// ============================================================================
{
  // computing the points defining the boundary polygon 

  vector<Point2D> bpoints; // "boundary points"
  if (tesselate_boundary) {
    for (unsigned int i = 0; i != num_corners; ++i) {
      const auto tseg = tesselateSegment(polygon[i], 
					 polygon[(i+1) % num_corners], 
					 vdist);
      bpoints.insert(bpoints.end(), tseg.begin(), tseg.end()-1);
    }
  } else {
    bpoints = vector<Point2D>(polygon, polygon + num_corners);
  }

  // Choosing a set of interior points which we can later move around to optimal
  // locations 
  vector<Point2D> ipoints = init_startpoints(polygon, num_corners, vdist); 

  //ipoints = vector<Point2D>(ipoints.begin(), ipoints.begin()+2); // @@@@
  
  //Optimizing position of interior points
  optimize_interior_points(&bpoints[0], (uint)bpoints.size(),
  			   &ipoints[0], (uint)ipoints.size(),
  			   vdist * 1.5); //1); // * 1.5 @@
  
  // // @@@@@@computing and reporting internal energy
  // const double R = 1.5 * vdist;
  // const auto e = polygon_energy(&bpoints[0], (unsigned int)bpoints.size(),
  //       			&ipoints[0], (unsigned int)ipoints.size(), 
  //       			R);

  // // @@@@@ computing numerical derivatives
  // auto dpoints(ipoints);
  // double tiny = 1e-7;
  // dpoints[1][0] = dpoints[1][0] + tiny;
  // const auto dxe = polygon_energy(&bpoints[0], (unsigned int)bpoints.size(),
  //       			&dpoints[0], (unsigned int)dpoints.size(), 
  //       			R);
  // double numdiff_x = (e.val-dxe.val)/-tiny;
  // dpoints = ipoints;
  // dpoints[1][1] = dpoints[1][1] + tiny;
  // const auto dye = polygon_energy(&bpoints[0], (unsigned int)bpoints.size(),
  //                                 &dpoints[0], (unsigned int)dpoints.size(), 
  //                                 R);
  // double numdiff_y = (e.val-dye.val)/-tiny;

  
  
  // cout << "Energy is:" << e.val << endl;

  // Triangulating points
  vector<Point2D> points(bpoints);
  points.insert(points.end(), ipoints.begin(), ipoints.end());

  const auto tris = triangulate_domain(&points[0],
				       (uint)bpoints.size(),
				       (uint)points.size(),
				       3*vdist);
  return {points, tris};
}

// ============================================================================
Mesh3D tesselatePolyhedron3D(const Point3D* const bpoints,
                             const unsigned int num_bpoints,
                             const Triangle* const btris,
                             const unsigned int num_btris,
                             const double vdist)
// ============================================================================
{
  // choosing a set of interior points
  vector<Point3D> ipoints = init_startpoints(bpoints, num_bpoints, btris,
                                             num_btris, vdist);

  //ipoints = {{0.1, 0.1, 0.1}, {0.9, 0.1, 0.1}, {0.5, 0.5, 0.5}};
  //ipoints = {{0.1, 0.1, 0.1}}; //,{0.5, 0.5, 0.5}};//,
  //ipoints = {{0.0, 0.1, 0.1}}; //,{0.5, 0.5, 0.5}};//, 
  //ipoints = {{0.5, 0.5, 0.5}, {0.11, 0.11, 0.2}, {0.51, 0.5, 0.5}, {0.1, 0.2, 0.3}}; // @@@ hack  
  // for (int i =0; i != (int)ipoints.size(); ++i)
  //   cout << ipoints[i]; @@

  // ofstream os0("krull0.mat");
  // copy(ipoints.begin(), ipoints.end(), ostream_iterator<Point3D>(os0, " "));
  // os0.close(); // @@@
  
  
  // optimizing position of interior points.  We use a value of 'vdist' slightly
  // higher than what the interpoint distance goal is, to avoid points becoming
  // 'completely disconnected' from each other.
  optimize_interior_points(bpoints, num_bpoints, btris, num_btris,
                           &ipoints[0], (uint)ipoints.size(), vdist*2);

  cout << "Finished optimization of " << ipoints.size() << " internal points." << endl;

  // ofstream os("krull.mat");
  // copy(ipoints.begin(), ipoints.end(), ostream_iterator<Point3D>(os, " "));
  // os.close(); // @@@


  
  // // @@@@@@computing and reporting internal energy
  // const double R = 2 * vdist;
  // const auto e = polyhedron_energy(bpoints, num_bpoints, btris, num_btris, 
  //                                  &ipoints[0], (unsigned int)ipoints.size(), R);


  // optimize_interior_points(bpoints, num_bpoints, btris, num_btris,
  //                          &ipoints[0], (uint)ipoints.size(), vdist*2); // @@ vdist=1?

  // cout << "Finished optimization of " << ipoints.size() << " internal points." << endl;


  
  // // @@@@@@computing and reporting internal energy
  // const auto e2 = polyhedron_energy(bpoints, num_bpoints, btris, num_btris, 
  //                                  &ipoints[0], (unsigned int)ipoints.size(), R);

  // double tiny2 = 1e-5;
  // vector<Point3D> krull(ipoints);
  // for (int i = 1; i != 2; ++i)
  //   for (int d = 0; d != 3; ++d)
  //     krull[i][d] -= e2.der[i][d] * tiny2;

  // //const auto e3 = polyhedron_energy(bpoints, num_bpoints, btris, num_btris, &krull[0], krull.size(), R);
  
  
  // // optimize_interior_points(bpoints, num_bpoints, btris, num_btris,
  // //                          &ipoints[0], (uint)ipoints.size(), vdist*2); // @@ vdist=1?

  // // cout << "Finished optimization of " << ipoints.size() << " internal points." << endl;


  
  // // // @@@@@@computing and reporting internal energy
  // // const auto e3 = polyhedron_energy(bpoints, num_bpoints, btris, num_btris, 
  // //                                  &ipoints[0], (unsigned int)ipoints.size(), R);





  
  // // @@@@@ computing numerical derivatives
  // auto dpoints(ipoints);
  // int pt_ix = 1; // 0;
  // double tiny = 1e-7;
  // dpoints[pt_ix][0] = dpoints[pt_ix][0] + tiny;
  // const auto dxe = polyhedron_energy(bpoints, num_bpoints, btris, num_btris,
  //       			&dpoints[0], (unsigned int)dpoints.size(), R);
  // double numdiff_x = (e.val-dxe.val)/-tiny;
  // dpoints = ipoints;
  // dpoints[pt_ix][1] = dpoints[pt_ix][1] + tiny;
  // const auto dye = polyhedron_energy(bpoints, num_bpoints, btris, num_btris,
  //                                 &dpoints[0], (unsigned int)dpoints.size(), R);
  // double numdiff_y = (e.val-dye.val)/-tiny;

  // dpoints = ipoints;
  // dpoints[pt_ix][2] = dpoints[pt_ix][2] + tiny;
  // const auto dze = polyhedron_energy(bpoints, num_bpoints, btris, num_btris,
  //                                 &dpoints[0], (unsigned int)dpoints.size(), R);
  // double numdiff_z = (e.val-dze.val)/-tiny;

  
  // ofstream os("krull.mat");
  // copy(ipoints.begin(), ipoints.end(), ostream_iterator<Point3D>(os, " "));
  // os.close(); // @@@

  // const auto tull = polyhedron_energy(bpoints, btris, num_btris,
  //                                     &ipoints[0], (uint)ipoints.size(), vdist * 4);
  
  // In case boundary has much larger faces than 'vdist', we need to allow for
  // larger distances when constructing the tetrahedrons.  We therefore run a
  // check here, and adjust the 'vdist' parameter upwards if the boundary faces
  // require it.
  vector<double> lengths(num_btris, 0);
  transform(btris, btris + num_btris, lengths.begin(), [bpoints] (const Triangle& t) {
      return max({dist(bpoints[t[0]], bpoints[t[1]]),
                  dist(bpoints[t[1]], bpoints[t[2]]),
                  dist(bpoints[t[2]], bpoints[t[0]])});});
  const double L = *max_element(lengths.begin(), lengths.end());
  
  // constructing tets from points
  vector<Point3D> points(bpoints, bpoints + num_bpoints);
  points.insert(points.end(), ipoints.begin(), ipoints.end());
  const auto tets = construct_tets(&points[0], (uint)points.size(),
                                   btris, num_btris,
                                   10 * max(vdist,L)); 
  return{points, tets};
}
};

namespace {

// ----------------------------------------------------------------------------
vector<Point3D> init_startpoints(const Point3D* const bpoints,
                                 const unsigned int num_bpoints,
                                 const Triangle* btris,
                                 const unsigned int num_btris,
                                 const double vdist)
// choose some initial startpoints as a basis for subsequent optimization
// ----------------------------------------------------------------------------
{
  const array<double, 6> bbox = bounding_box_3D(bpoints, num_bpoints);
  const double bbox_lx = bbox[1] - bbox[0];
  const double bbox_ly = bbox[3] - bbox[2];
  const double bbox_lz = bbox[5] - bbox[4];

  // computing amount of points needed to approximately fill the bounding box,
  // where points are approximately equidistant by 'vdist'.
  // const double box_vol = bbox_lx * bbox_ly * bbox_lz;
  // const int N_bbox = (int)ceil(box_vol / (vdist * vdist * vdist));
  // const double Rxy = bbox_lx / bbox_ly;
  // const double Rxz = bbox_lx / bbox_lz;
  // const double Ryz = bbox_ly / bbox_lz;
  // const int nx = (int)ceil(pow(N_bbox * Rxy * Rxz, 1/3.));
  // const int ny = (int)ceil(pow(N_bbox / Rxy * Ryz, 1/3.));
  // const int nz = (int)ceil(pow(N_bbox / Rxz / Ryz, 1/3.));
  const int nx = (int)floor(bbox_lx/vdist);
  const int ny = (int)floor(bbox_ly/vdist);
  const int nz = (int)floor(bbox_lz/vdist);

  // constructing regular grid of points within the bounding box
  const vector<Point3D> gridpoints =
    generate_grid_3D(Point3D {bbox[0], bbox[2], bbox[4]},
                     Point3D {bbox[1], bbox[3], bbox[5]},
                     nx, ny, nz);

    // generate_grid_3D(Point3D {bbox[0] + vdist/2, bbox[2] + vdist/2, bbox[4] + vdist/2},
    //                  Point3D {bbox[1] - vdist/2, bbox[3] - vdist/2, bbox[5] - vdist/2},
    //                  nx, ny, nz);

  // keeping the points that fall within the shell
  const double TOL_FAC = 0.25 * vdist; // do not keep points too close to boundary
  vector<Point3D> result = inside_shell(&gridpoints[0],
                                        (unsigned int) gridpoints.size(),
                                        bpoints,
                                        btris,
                                        num_btris,
                                        TOL_FAC);

  // add perturbation to avoid 'symmetry locking' during the optimization stage
  const double PM = 1.0e-1 * min({bbox_lx/nx, bbox_ly/ny, bbox_lz/nz});
  transform(result.begin(), result.end(), result.begin(), [PM] (const Point3D& p) {
      return Point3D {p[0] + random_uniform(-PM, PM),
                      p[1] + random_uniform(-PM, PM),
                      p[2] + random_uniform(-PM, PM)};});
  return result;
}
   
// ----------------------------------------------------------------------------  
vector<Point2D> init_startpoints(const Point2D* const polygon,
				 const unsigned int num_corners,
				 const double vdist)
// choose some initial startpoints as a basis for subsequent optimization  
// ----------------------------------------------------------------------------
{
  //const double poly_area = polygon_area(polygon, num_corners);
  const array<double, 4> bbox = bounding_box_2D(polygon, num_corners);
  const double bbox_lx = bbox[1] - bbox[0];
  const double bbox_ly = bbox[3] - bbox[2];
  const double bbox_area =  bbox_lx * bbox_ly;

  // computing amount of points needed to approximately cover the bounding box, where points
  // are approximately equidistant by 'vdist'
  const int N_bbox = (int)ceil(bbox_area * 3 / (PI * vdist * vdist));
  const int nx = (int)floor(sqrt(N_bbox * bbox_lx / bbox_ly));
  const int ny = (int)floor(N_bbox/nx);

  // constructing regular grid of points within the bounding box
  const vector<Point2D> gridpoints =
    generate_grid_2D(Point2D {bbox[0] + vdist/2, bbox[2] + vdist/2},
		     Point2D {bbox[1] - vdist/2, bbox[3] - vdist/2},
		     nx, ny);
  
  // keeping the points that fall within the original polygon
  const double TOL_FAC = 0.25 * vdist; // do not keep points too close to the boundary
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
void optimize_interior_points(const Point3D* bpoints,
                              const unsigned int num_bpoints,
                              const Triangle* const btris,
                              const unsigned int num_btris,
                              Point3D* const ipoints,
                              const unsigned int num_ipoints,
                              const double vdist)
// ----------------------------------------------------------------------------
{
  // setting up function to minimize (energy function):
  auto efun = PolyhedronEnergyFunctor(bpoints,num_bpoints, btris, num_btris,
                                      ipoints, num_ipoints, vdist);

  // setting up function minimizer
  Go::FunctionMinimizer<PolyhedronEnergyFunctor>
    funcmin(num_ipoints * 3, efun, (double* const)&ipoints[0], 1e-1); // @@ tolerance?

  // do the minimization
  Go::minimise_conjugated_gradient(funcmin);

  // copying results back.  The cast is based on the knowledge that Points2D
  // consist of POD and that point coordinates are stored contiguously in memory
  // as (x1, y1, z1, x2, y2, z2, ...).
  double* const target = (double* const) ipoints;
  std::copy(funcmin.getPar(), funcmin.getPar() + funcmin.numPars(), target);
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
    funcmin(num_ipoints * 2, efun, (double* const)&ipoints[0], 1e-1);//1e-1); // @@ TOLERANCE?
  
  // do the minimization
  Go::minimise_conjugated_gradient(funcmin);

  // copying results back.  The cast is based on the knowledge that Points2D
  // consists of POD and that points are stored contiguously in memory
  // (as x1,y1,x2,y2, ...)
  double* const target = (double* const) ipoints;
  std::copy(funcmin.getPar(), funcmin.getPar() + funcmin.numPars(), target);
}

// ----------------------------------------------------------------------------
PolyhedronEnergyFunctor::PolyhedronEnergyFunctor(const Point3D* const bpoints,
                                                 const unsigned int num_bpoints,
                                                 const Triangle* const btris,
                                                 const unsigned int num_btris,
                                                 const Point3D* const ipoints,
                                                 const unsigned int num_ipoints,
                                                 const double radius)
// ----------------------------------------------------------------------------
  : bpoints_(bpoints), nb_(num_bpoints), btris_(btris), nt_(num_btris),
    ni_(num_ipoints), r_(radius), bbox_(bounding_box_3D(bpoints, num_bpoints))  
{ }

// ----------------------------------------------------------------------------
double PolyhedronEnergyFunctor::operator()(const double* const arg) const
// ----------------------------------------------------------------------------
{
  if (!use_cached(arg)) 
    update_cache(arg);
  return cached_result_.val;
}

// ----------------------------------------------------------------------------
void PolyhedronEnergyFunctor::grad(const double* arg, double* grad) const
// ----------------------------------------------------------------------------
{
  if (!use_cached(arg))
    update_cache(arg);
  const double* const dp = (const double* const)&cached_result_.der[0];
  copy(dp, dp + 3 * ni_, grad);
}

// ----------------------------------------------------------------------------
double PolyhedronEnergyFunctor::minPar(int n) const
// ----------------------------------------------------------------------------
{
  return bbox_[(n % 3) * 2];
}

// ----------------------------------------------------------------------------
double PolyhedronEnergyFunctor::maxPar(int n) const
// ----------------------------------------------------------------------------
{
  return bbox_[(n % 3) * 2 + 1];
}

// ----------------------------------------------------------------------------
bool PolyhedronEnergyFunctor::use_cached(const double* const arg) const
// ----------------------------------------------------------------------------
{
  return (cached_arg_.size() > 0) &&
         (std::equal(arg, arg + 3 * ni_, &cached_arg_[0]));
}

// ----------------------------------------------------------------------------  
void PolyhedronEnergyFunctor::update_cache(const double* const arg) const
// ----------------------------------------------------------------------------
{
  if (cached_arg_.empty())
    cached_arg_.resize(3 * ni_);
  copy(arg, arg + 3 * ni_, &cached_arg_[0]);
  cached_result_ = polyhedron_energy(bpoints_, nb_, btris_, nt_,
                                     (const Point3D* const) &cached_arg_[0],
                                     ni_, r_);
}

// ----------------------------------------------------------------------------    
PolygonEnergyFunctor::PolygonEnergyFunctor(const Point2D* const polygon,
					   const unsigned int num_corners,
					   const unsigned int num_ipoints,
					   const double radius)
// ----------------------------------------------------------------------------    
  : poly_(polygon), nc_(num_corners), ni_(num_ipoints), r_(radius),
    bbox_(bounding_box_2D(polygon, num_corners))
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
