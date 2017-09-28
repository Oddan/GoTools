#include "parametric_object_energies.h"
#include "tesselate_utils.h"
#include "interpoint_distances.h"
#include "distance_energy.h"

#include <array>
#include <tuple>
#include <vector>

using namespace Go;
using namespace std;
using namespace TesselateUtils;

namespace {
  
typedef vector<Point3D> PointVec3D;
typedef shared_ptr<const Go::ParamSurface> SurfPtr;

tuple<PointVec3D, PointVec3D, PointVec3D>
compute_3D_points_and_derivs(const Point2D* const pars,
                             const uint num_pars,
                             const SurfPtr surf);

PointVec3D compute_3D_points(const Point2D* const pars,
                           const uint num_pars,
                           const SurfPtr surf);
                                                                 
ValAndDer<Point2D> internal_energy(const PointVec3D& pts,
                                   const PointVec3D& uder,
                                   const PointVec3D& vder,
                                   const double vdist);
  
void accumulate_energy(const double dist,
                       const double vdist,
                       const Point3D& p1,
                       const Point3D& p2,
                       const Point3D& uder1,
                       const Point3D& uder2,
                       const Point3D& vder1,
                       const Point3D& vder2,
                       double& energy_acc,
                       Point2D& p1_der_acc,
                       Point2D& p2_der_acc);

ValAndDer<Point2D> boundary_energy(const PointVec3D& ipts3D,
                                   const PointVec3D& ipts_uder,
                                   const PointVec3D& ipts_vder,
                                   const Point2D* const ipts_par,
                                   const PointVec3D& bpts3D,
                                   const Point2D* const bpts_par,
                                   const double vdist,
                                   const ClippedGrid<2>* const cgrid);

ClippedDomainType point_domain_type(const Point2D& pt,
                                    const Point2D* const bpoints,
                                    const uint num_bpoints,
                                    const ClippedGrid<2>& cgrid);
};

namespace TesselateUtils {

// ----------------------------------------------------------------------------
ValAndDer<Point2D>
parametric_surf_energy(const shared_ptr<const Go::ParamSurface> surf,
                       const Point2D* const bpoints, // boundary polygon corners
                       const uint num_bpoints, // number of boundary corners
                       const Point2D* const ipoints, // internal points
                       const uint num_ipoints, // number of internal points
                       const double vdist,
                       const ClippedGrid<2>* const cgrid)
// ----------------------------------------------------------------------------
{
  // computing interior 3D points and their derivatives
  PointVec3D ip3D(num_ipoints), iu_der(num_ipoints), iv_der(num_ipoints);
  tie(ip3D, iu_der, iv_der) = compute_3D_points_and_derivs(ipoints, num_ipoints, surf);

  // compute boundary 3D points (derivatives not needed as boundary is fixed)
  const PointVec3D bp3D = compute_3D_points(bpoints, num_bpoints, surf);
  
  // compute potential energy between internal points
  const ValAndDer<Point2D> E_int = internal_energy(ip3D, iu_der, iv_der, vdist);

  // compute boundary energy (energy from interaction between internal points
  // and boundary points)
  const ValAndDer<Point2D> E_bnd = boundary_energy(ip3D, iu_der, iv_der, ipoints, 
                                                   bp3D, bpoints, 
                                                   vdist/1.5, cgrid);

  // adding up components and returning results
  ValAndDer<Point2D> E_tot = E_int;
  const double BND_FAC = 2; // increase penalty for approaching border
  E_tot.val += BND_FAC * E_bnd.val;
  for (uint i = 0; i != (uint)E_tot.der.size(); ++i)
    E_tot.der[i] += E_bnd.der[i] * BND_FAC;

  return E_tot;
}

  
}; //end namespace TesselateUtils


namespace {

// ----------------------------------------------------------------------------
ValAndDer<Point2D> boundary_energy(const PointVec3D& ipts3D,
                                   const PointVec3D& ipts_uder,
                                   const PointVec3D& ipts_vder,
                                   const Point2D* const ipts_par,
                                   const PointVec3D& bpts3D,
                                   const Point2D* const bpts_par,
                                   const double vdist,
                                   const ClippedGrid<2>* const cgrid)
// ----------------------------------------------------------------------------
{
  const auto num_ipts = ipts3D.size();
  const auto num_bpts = bpts3D.size();
  ValAndDer<Point2D> result {0, vector<Point2D>(num_ipts, {0.0, 0.0})};

  // detect points lying outside boundary
  vector<ClippedDomainType> point_status(num_ipts);
  if (cgrid == nullptr) 
    // simple, brute-force method that only distinguishes between points that
    // are clearly outside and other points
    transform(ipts_par, ipts_par + num_ipts, point_status.begin(),
              [&] (const Point2D& p) {
                bool on_bnd = false;
                bool inside = inpolygon(p, bpts_par, (uint)num_bpts, 0, on_bnd);
                return ((!inside) && (!on_bnd)) ? OUTSIDE : UNDETERMINED;
              });
  else
    // benefit from the precomputed informatoin in 'cgrid' to significantly
    // speed up computations
    transform(ipts_par, ipts_par + num_ipts, point_status.begin(),
              [&] (const Point2D& p) {
                return point_domain_type(p, bpts_par, (uint)num_bpts, *cgrid);
              });

  // looping across boundary segments and adding their energy contribution
  assert(false);
  
  // // looping across boundary segments and adding their energy contributions
  // for (uint i = 0; i != num_bpoints; ++i)
  //   add_boundary_contribution(bpoints[i], bpoints[(i+1)%num_bpoints],
  //       		      ipoints, num_ipoints, &point_status[0], vdist, result);

  // // adding energy for points having exited the domain
  // for (uint i = 0; i != num_ipoints; ++i) {
  //   if (point_status[i] == OUTSIDE)
  //     add_outside_penalty_energy(i, ipoints, bpoints, num_bpoints, vdist, result);
  // }
  
  return result;
}
  
// ----------------------------------------------------------------------------  
ValAndDer<Point2D> internal_energy(const PointVec3D& pts,
                                   const PointVec3D& uder,
                                   const PointVec3D& vder,
                                   const double vdist)
// ----------------------------------------------------------------------------  
{
  const auto dists = interpoint_distances(&pts[0], (uint)pts.size(), vdist);
  
  // accumulating energies and storing total value and partial derivatives in
  // 'result'
  ValAndDer<Point2D> result {0, vector<Point2D>((uint)pts.size(), {0.0, 0.0})};
  
  for (const auto& d : dists)
    accumulate_energy(d.dist, vdist,
                      pts[d.p1_ix], pts[d.p2_ix],
                      uder[d.p1_ix], uder[d.p2_ix],
                      vder[d.p1_ix], vder[d.p2_ix],
                      result.val,
                      result.der[d.p1_ix],
                      result.der[d.p2_ix]);
  return result;
}
  
// ----------------------------------------------------------------------------
void accumulate_energy(const double dist,
                       const double vdist,
                       const Point3D& p1,
                       const Point3D& p2,
                       const Point3D& uder1,
                       const Point3D& uder2,
                       const Point3D& vder1,
                       const Point3D& vder2,
                       double& energy_acc,
                       Point2D& p1_der_acc,
                       Point2D& p2_der_acc)
// ----------------------------------------------------------------------------
{
  const array<double, 2> e = distance_energy(dist, vdist);

  // accumulating total energy
  energy_acc += e[0];

  // Adding contribution to partial derivatives for the two points involved.
  // @@ Creation of temporary objects below.  Could be optimized!
  const Point3D dvec = (p2 - p1) * e[1] / dist; 

  // derivatives
  p1_der_acc[0] -= dvec * uder1;
  p1_der_acc[1] -= dvec * vder1;

  p2_der_acc[0] += dvec * uder2;
  p2_der_acc[1] += dvec * vder2;
    
}
  
// ----------------------------------------------------------------------------
tuple<PointVec3D, PointVec3D, PointVec3D>
compute_3D_points_and_derivs(const Point2D* const pars,
                             const uint num_pars,
                             const SurfPtr surf)
// ----------------------------------------------------------------------------
{
  tuple<PointVec3D, PointVec3D, PointVec3D> result;

  get<0>(result).resize(num_pars);
  get<1>(result).resize(num_pars);
  get<2>(result).resize(num_pars);

  vector<Go::Point> tmp_pts(3);
  for (uint i = 0; i != num_pars; ++i) {
    surf->point(tmp_pts, pars[i][0], pars[i][1], 1);

    copy(tmp_pts[0].begin(), tmp_pts[0].end(), &(get<0>(result)[i])[0]);
    copy(tmp_pts[1].begin(), tmp_pts[1].end(), &(get<1>(result)[i])[0]);
    copy(tmp_pts[2].begin(), tmp_pts[2].end(), &(get<2>(result)[i])[0]);
  }
  return result;
}

// ----------------------------------------------------------------------------  
PointVec3D compute_3D_points(const Point2D* const pars,
                           const uint num_pars,
                           const SurfPtr surf)
// ----------------------------------------------------------------------------
{
  PointVec3D result(num_pars);
  transform(pars, pars + num_pars, result.begin(), [surf] (const Point2D& p) {
      const auto tmp = surf->point(p[0], p[1]);
      return Point3D {tmp[0], tmp[1], tmp[2]};
    });
  return result;
}

// ----------------------------------------------------------------------------  
inline ClippedDomainType point_domain_type(const Point2D& pt,
                                           const Point2D* const bpoints,
                                           const uint num_bpoints,
                                           const ClippedGrid<2>& cgrid)
// ----------------------------------------------------------------------------  
{
  const double& dx = cgrid.cell_len[0];
  const double& dy = cgrid.cell_len[1];
  const uint ix = min((uint)max(floor((pt[0] - cgrid.bbox[0])/dx), 0.0), cgrid.res[0]-1);
  const uint iy = min((uint)max(floor((pt[1] - cgrid.bbox[2])/dy), 0.0), cgrid.res[1]-1);

  const auto type = (cgrid.type[iy * cgrid.res[0] + ix]);

  if (type == FAR_INSIDE || type == CLOSE_INSIDE || type == OUTSIDE)
    return type;

  // if we got here, type is either UNDETERMINED or INTERSECTED.  We must do an
  // explicit computation to determine its status vis-a-vis the boundary.
  bool on_bnd = false;
  const bool inside = inpolygon(pt, bpoints, num_bpoints, 0, on_bnd);

  return (inside || on_bnd) ? CLOSE_INSIDE : OUTSIDE;
}
  
  

  
}; // end anonymous namespace
