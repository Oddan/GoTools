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
  u_points->initNeighbours();
  
  // Parametrize the boundary
  PrParametrizeBdy pb;
  pb.setParamKind(PrCHORDLENGTHBDY);
  pb.attach(u_points);
  pb.parametrize();

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
