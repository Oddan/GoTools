#include "triangulate_domain.h"
#include "tesselate_utils.h"

#include "GoTools/parametrization/PrFastUnorganized_OP.h"
#include "GoTools/parametrization/PrParametrizeInt.h"
#include "GoTools/parametrization/PrParametrizeBdy.h"
// #include "GoTools/parametrization/PrPrmShpPres.h"
// #include "GoTools/parametrization/PrPrmUniform.h"
// #include "GoTools/parametrization/PrPrmLeastSquare.h"
// #include "GoTools/parametrization/PrPrmEDDHLS.h"
#include "GoTools/parametrization/PrPrmMeanValue.h"

#include <memory>

using namespace Go;
using namespace std;
using namespace TesselateUtils;

namespace {
  vector<Point2D> compute_2D_paramerization(const Point3D* const points,
                                            const uint num_bpoints,
                                            const uint tot_num_points,
                                            const double vdist);
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

  Point3D u = normal ^ sel_axis;  u = u / norm(u);
  Point3D v = normal ^ u;;   v = v / norm(v);
  vector<Point3D> tmp_bnd_pts {Point3D {0.0, 0.0, 0.0}, u, v};
  tmp_bnd_pts.insert(tmp_bnd_pts.end(), points, points + num_bpoints);
  vector<Point2D> bnd_par = transform_to_2D<Point2D, Point3D>(tmp_bnd_pts);
  for (uint i = 0; i != num_bpoints; ++i) {
    u_points->setU(i + num_ipoints, bnd_par[i+3][0]);
    u_points->setV(i + num_ipoints, bnd_par[i+3][1]);
  }

  // PrParametrizeBdy pb;
  // pb.setParamKind(PrCHORDLENGTHBDY);
  // pb.attach(u_points);
  // pb.parametrize();

  //Parametrize the interior
  PrPrmMeanValue pi;
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

}; // end anonymous namespace 
