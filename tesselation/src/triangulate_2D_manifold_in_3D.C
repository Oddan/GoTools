#include "triangulate_domain.h"
#include "tesselate_utils.h"

#include "GoTools/parametrization/PrFastUnorganized_OP.h"
#include "GoTools/parametrization/PrParametrizeInt.h"
#include "GoTools/parametrization/PrParametrizeBdy.h"
#include "GoTools/parametrization/PrPrmShpPres.h"
// #include "GoTools/parametrization/PrPrmUniform.h"
// #include "GoTools/parametrization/PrPrmLeastSquare.h"
// #include "GoTools/parametrization/PrPrmEDDHLS.h"
// #include "GoTools/parametrization/PrPrmMeanValue.h"

#include <memory>

using namespace Go;
using namespace std;
using namespace TesselateUtils;

namespace {
  vector<Point2D> compute_2D_paramerization(const Point3D* const points,
                                            const uint num_bpoints,
                                            const uint tot_num_points,
                                            const double vdist);

  vector<Point2D> compute_boundary_param(const vector<Point3D>& pts3D);
  vector<double> compute_approx_angles(const vector<Point3D>& pts3D);
  vector<Point2D> compute_2D_segments(const vector<Point3D>& pts3D,
                                      const vector<double>& angles);
  
  double compute_angle(const Point3D& p1, const Point3D& p2,
                       const Point3D& p3, const Point3D& midpt);
  
}; // end anonymous namespace 


namespace TesselateUtils {

// ============================================================================
std::vector<Triangle> triangulate_2D_manifold_in_3D(const Point3D* const points,
                                                    const uint num_bpoints,
                                                    const uint tot_num_points,
                                                    const double vdist)
// ============================================================================  
{
  // first, establish a reasonable 2D parameterization of the 3D points
  const auto uv = compute_2D_paramerization(points, num_bpoints, tot_num_points, vdist);
    
  // then triangulate them
  const double new_vdist = 4 / sqrt(max(tot_num_points - num_bpoints, (uint)1));
  return triangulate_domain(&uv[0], num_bpoints, tot_num_points, new_vdist);
}  

}; // end namespace TesselateUtils

namespace {

// ----------------------------------------------------------------------------
vector<Point2D> compute_2D_paramerization(const Point3D* const points,
                                          const uint num_bpoints,
                                          const uint tot_num_points,
                                          const double vdist)
// ----------------------------------------------------------------------------
{
  // the PrFastUnorganized_OP object below assumes boundary points are listed
  // last, so we need to reverse our usual order of boundary points and interior
  // points here
  vector<Point3D> tmp_points(points + num_bpoints, points + tot_num_points);
  tmp_points.insert(tmp_points.end(), points, points + num_bpoints);

  const uint num_ipoints = tot_num_points - num_bpoints;
  shared_ptr<PrFastUnorganized_OP> 
    u_points(new PrFastUnorganized_OP(tot_num_points,
                                      num_ipoints,
                                      (double*) &tmp_points[0]));
  u_points->setRadius(vdist);
  u_points->useRadius();
  u_points->initNeighbours();
  
  // Parametrize the boundary
  Point3D normal = compute_polygon_normal(points, num_bpoints);
  uint nix = uint(min_element(&normal[0], &normal[0]+3) - &normal[0]);
  Point3D sel_axis = {0, 0, 0};
  sel_axis[nix] = 1;

  // Point3D u = normal ^ sel_axis;  u = u / norm(u);
  // Point3D v = normal ^ u;;   v = v / norm(v);
  // vector<Point3D> tmp_bnd_pts {Point3D {0.0, 0.0, 0.0}, u, v};
  // tmp_bnd_pts.insert(tmp_bnd_pts.end(), points, points + num_bpoints);
  // vector<Point2D> bnd_par = transform_to_2D<Point2D, Point3D>(tmp_bnd_pts);
  // for (uint i = 0; i != num_bpoints; ++i) {
  //   u_points->setU(i + num_ipoints, bnd_par[i+3][0]);
  //   u_points->setV(i + num_ipoints, bnd_par[i+3][1]);
  // }


  const vector<Point3D> tmp_bnd_pts(points, points + num_bpoints);
  const vector<Point2D> bnd_par = compute_boundary_param(tmp_bnd_pts);
  
  for (uint i = 0; i != num_bpoints; ++i) {
    u_points->setU(i + num_ipoints, bnd_par[i][0]);
    u_points->setV(i + num_ipoints, bnd_par[i][1]);
  }

  // PrParametrizeBdy pb;
  // pb.setParamKind(PrCHORDLENGTHBDY);
  // pb.attach(u_points);
  // pb.parametrize();

  //Parametrize the interior
  //PrPrmMeanValue pi;
  PrPrmShpPres pi;
  pi.attach(u_points);
  pi.parametrize();

  // returning result
  vector<Point2D> result(tot_num_points);
  for (uint i = 0; i != num_bpoints; ++i)
    // paramter values of boundary points
    result[i] = {u_points->getU(i+num_ipoints), u_points->getV(i+num_ipoints)};

  for (uint i = 0; i != num_ipoints; ++i)
    // parameter values for interior points
    result [i + num_bpoints] = {u_points->getU(i), u_points->getV(i)};
                 
  return result;  
}

// ----------------------------------------------------------------------------
vector<Point2D> compute_boundary_param(const vector<Point3D>& pts3D)
// ----------------------------------------------------------------------------
{
  // Compute approximate 2D parameterization of boundary points, that tries as
  // much as possible to reflect the angles and distances between the 3D
  // boundary points.
  
  // computed midpoint (heuristically used to determine whether an angle bends
  // "inwards" or "outwards)
  const int N = (int)pts3D.size();
  
  // compute angles and ensure they sum up to 2 pi
  const vector<double> angles = compute_approx_angles(pts3D);

  // compute segment lengths and directions, and adjust to ensure endpoints meet
  const vector<Point2D> segs = compute_2D_segments(pts3D, angles);

  // compute 2D coordinates by adding up the segments
  vector<Point2D> result(N);
  result[0] = Point2D {0.0, 0.0};
  for (int i = 1; i != N; ++i) 
    result[i] = result[i-1] + segs[i];

  return result;
  
}

// ----------------------------------------------------------------------------  
vector<double> compute_approx_angles(const vector<Point3D>& pts3D)
// ----------------------------------------------------------------------------    
{
  const int N = (int)pts3D.size();
  const double PI = 3.1415926536;    
  const Point3D midpt = accumulate(pts3D.begin(),
                                   pts3D.end(),
                                   Point3D {0.0, 0.0, 0.0}) / (double) N;
  
  vector<double> angles(pts3D.size());
  for (int i = 0; i != N; ++i) {
    angles[i] = compute_angle(pts3D[(i-1+N)%N], pts3D[i], pts3D[(i+1)%N], midpt);
  }
  const double initial_angle_sum = accumulate(angles.begin(), angles.end(), 0.0);
  transform(angles.begin(), angles.end(), angles.begin(),
            [&](double x) {return x * 2 * PI / initial_angle_sum;});

  // the angles in 'angles' should now sum up to 2 PI
  return angles;
}

// ----------------------------------------------------------------------------  
double compute_angle(const Point3D& p1, const Point3D& p2, const Point3D& p3,
                     const Point3D& midpt)
// ----------------------------------------------------------------------------  
{
  const Point3D v1 = p2 - p1;
  const Point3D v2 = p3 - p2;
  const Point3D r  = p1 - midpt;

  const double v1_norm = norm(v1);
  const double v2_norm = norm(v2);

  const Point3D n = v1 ^ v2;
  const Point3D c =  r ^ v1; 

  const double sign = (n*c > 0) ? 1 : -1;

  const double result = sign * acos(v1 * v2 / (v1_norm * v2_norm));
  return result;
}

// ----------------------------------------------------------------------------
vector<Point2D> compute_2D_segments(const vector<Point3D>& pts3D,
                                    const vector<double>& angles)
// ----------------------------------------------------------------------------
{
  const int N = (int)pts3D.size();
  vector<Point2D> result(N);
  
  // First, compute segments that have the right orientation and length
  double cur_angle = 0;
  for (int i = 0; i != N; ++i) {
    const double len = sqrt(dist2(pts3D[i],
                                  pts3D[(i+1)%N]));
    cur_angle += angles[i];
    result[i] = {len * cos(cur_angle), len * sin(cur_angle)};
  }

  // Then, adjust segment lengths so that endpoints match in the end
  const double x_err = accumulate(result.begin(), result.end(), 0.0,
                                  [](double acc, const Point2D& p)
                                  {return acc + p[0];});
  const double y_err = accumulate(result.begin(), result.end(), 0.0,
                                  [](double acc, const Point2D& p)
                                  {return acc + p[1];});

  // distributing error evently across segments, to ensure a set of segments
  // that sum up to zero
  const double dx = x_err/N;
  const double dy = y_err/N;
  transform(result.begin(), result.end(), result.begin(),
            [dx, dy] (const Point2D& p)
            { return Point2D {p[0] - dx, p[1] - dy};});
  return result;
}

  
}; // end anonymous namespace 
