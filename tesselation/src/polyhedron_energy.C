#include "polyhedral_energies.h"
#include "interpoint_distances.h"
#include "tesselate_utils.h"
#include "common_defs.h"

using namespace std;
using namespace TesselateUtils;

namespace {

// // first entry in result array is the energy, the second the derivative
array<double, 2> energy(double dist, double R);

ValAndDer<Point3D> internal_energy(const Point3D* const ipoints,
                                   const unsigned int num_ipoints,
                                   const double vdist);
  
ValAndDer<Point3D> boundary_energy(const Point3D* const bpoints,
                                   const Triangle* const btris,
                                   const unsigned int num_btris,
                                   const Point3D* const ipoints,
                                   const unsigned int num_ipoints,
                                   const double vdist);

void add_boundary_contribution(const Point3D* const bpoints,
			       const Triangle& btri,
			       const Point3D* const ipoints,
			       const unsigned int num_ipoints,
			       const double vdist,
                               ValAndDer<Point3D>& result);

void accumulate_energy(const double dist,
		       const Point3D& p1,
		       const Point3D& p2,
		       const double vdist,
		       double& energy_acc,
		       Point3D& p1_der_acc,
		       Point3D& p2_der_acc,
                       const bool are_mirror_points = false);

}; // end anonymous namespace 

namespace TesselateUtils {

// ----------------------------------------------------------------------------
ValAndDer<Point3D> polyhedron_energy(const Point3D* const bpoints,
                                     const Triangle* const btris,
                                     const unsigned int num_btris,
                                     const Point3D* const ipoints,
                                     const unsigned int num_ipoints,
                                     const double vdist)
  // ----------------------------------------------------------------------------
{
  
  //compute internal energy (potential energy between internal points)
  const ValAndDer<Point3D> E_int = internal_energy(ipoints, num_ipoints, vdist);

  // compute boundary energy (energy from interaction between internal points
  // and boundary points, and from internal points and their mirror points)
  const ValAndDer<Point3D> E_bnd = boundary_energy(bpoints, btris, num_btris,
                                                   ipoints, num_ipoints,
                                                   vdist);

  // Adding up components and returning results
  ValAndDer<Point3D> E_tot = E_int;

  E_tot.val += E_bnd.val; 
  E_tot.der += E_bnd.der;

  return E_tot;
}

}; // end namespace Go


namespace {

// ----------------------------------------------------------------------------
ValAndDer<Point3D> internal_energy(const Point3D* const ipoints,
                                   const unsigned int num_ipoints,
                                   const double vdist)
// ----------------------------------------------------------------------------
{
  const auto dists = interpoint_distances(ipoints, num_ipoints, vdist, false);

  // accumulating energies and storing total value and partial derivatives in 'result'
  ValAndDer<Point3D> result {0, vector<Point3D>(num_ipoints, {0.0, 0.0, 0.0})};
  for (const auto& d : dists)
    accumulate_energy(d.dist, ipoints[d.p1_ix], ipoints[d.p2_ix], vdist, 
		      result.val, result.der[d.p1_ix], result.der[d.p2_ix]);

  return result;
}

// ----------------------------------------------------------------------------
void accumulate_energy(const double dist,
		       const Point3D& p1,
		       const Point3D& p2,
		       const double vdist,
		       double& energy_acc,
		       Point3D& p1_der_acc,
		       Point3D& p2_der_acc,
                       const bool are_mirror_points)
// ----------------------------------------------------------------------------  
{
  const array<double,2> e = energy(dist, vdist);

  // accumulating total energy
  energy_acc += e[0];
  
  // adding contribution to partial derivatives for the two points involved.
  // @@ Creation of temporary object on the line below.  If bottleneck, should
  // be rewritten.

  //If the points are mirror points of each other, all contributions to partial
  //derivatives should be doubled (and the accumulated values for the mirror
  //points are irrelevant)
  const double fac = are_mirror_points ? 2 : 1;
  
  const Point3D dvec = (p2 - p1) * (fac * e[1] / dist);
  p1_der_acc -= dvec;
  p2_der_acc += dvec;
  
}

// ----------------------------------------------------------------------------
ValAndDer<Point3D> boundary_energy(const Point3D* const bpoints,
                                   const Triangle* const btris,
                                   const unsigned int num_btris,
                                   const Point3D* const ipoints,
                                   const unsigned int num_ipoints,
                                   const double vdist)
// ----------------------------------------------------------------------------
{
  ValAndDer<Point3D> result {0, vector<Point3D>(num_ipoints, {0.0, 0.0, 0.0})};

  // looping across boundary faces and adding their energy contributions
  for (uint i = 0; i != num_btris; ++i) 
    add_boundary_contribution(bpoints, btris[i], ipoints, num_ipoints, vdist, result);
  
  return result;
}

// ----------------------------------------------------------------------------
void add_boundary_contribution(const Point3D* const bpts, // boundary points
			       const Triangle& btri,
			       const Point3D* const ipoints,
			       const unsigned int num_ipoints,
			       const double vdist,
                               ValAndDer<Point3D>& result)
{
  // identify the points within reach of the face
  const array<Point3D, 3> tripts { bpts[btri[0]], bpts[btri[1]], bpts[btri[2]]};
  const auto npoints = extract_from_range(ipoints, num_ipoints,
                                          [&tripts, &vdist] (const Point3D& p) {
                                            return point_on_triangle(p,
                                                                     tripts[0],
                                                                     tripts[1],
                                                                     tripts[2],
                                                                     vdist);});
  const auto& neigh_ixs = npoints.second;
  const auto& neigh_pts = npoints.first;

  // setting up a result structure only for the neigh points
  ValAndDer<Point3D> result_local {0, vector<Point3D>(neigh_pts.size(),
                                                      {0.0, 0.0, 0.0})};

  // computing the energy (and associated partial derivatives) between the
  // neighbor points and the two boundary point
  Point3D dummy; // used to store values we do not need
  for (uint i = 0; i != 3; ++i) {// loop over triangle corners
    const auto dvec = interpoint_distances(&neigh_pts[0], (uint)neigh_pts.size(),
                                           tripts[i], vdist);
    for (const DistanceEntry& d : dvec)
      accumulate_energy(d.dist, neigh_pts[d.p1_ix], tripts[i], vdist,
                        result_local.val, result_local.der[d.p1_ix], dummy);
  }
  
  // remapping results from result_local to result.  Multiplying by 1/6, since
  // the contribution from each corner will be counted as many times as the
  // valence of that corner. We do not keep track of valences; 6 is the typical
  // value for a regular mesh, which is why we divide it here.  @@ If this
  // simplification leads to bad results close to high valence (or low valence)
  // points, consider a more rigorous treatment.

  result.val += result_local.val * 0.1667; // 1/6
  for (uint i = 0; i != neigh_ixs.size(); ++i)
    result.der[neigh_ixs[i]] += result_local.der[i] * 0.1667;
  
  // Now, we compute the energy between neighbor points and their mirror images.
  // First, we remove the points that do not project perpendicularly to the
  // surface of the triangle itself.
  double dist_dummy;
  int sign_dummy;
  const double TOL = 1e-5 * vdist; // @@ should be enough?
  const auto per = extract_from_range(&neigh_pts[0], (uint)neigh_pts.size(),
                                      [&] (const Point3D& p) {
                                        return projects_to_triangle(p,
                                                                    &tripts[0],
                                                                    TOL,
                                                                    dist_dummy,
                                                                    sign_dummy);});
  const auto& per_pts_ixs = per.second;
  const auto& per_pts     = per.first;
  result_local.reset((uint)per_pts.size()); // reinitialize local result structure
 
  // computing mirror points and their energy contributions
  const auto mpoints = mirror_points_3D(&per_pts[0], (uint)per_pts.size(), &tripts[0]);
  const auto dvec = interpoint_distances(&per_pts[0], (uint)per_pts.size(),
                                         &mpoints[0], (uint)mpoints.size(), vdist);
  for (const DistanceEntry& d : dvec)  
    accumulate_energy(d.dist, per_pts[d.p1_ix], mpoints[d.p2_ix], vdist,
                      result_local.val, result_local.der[d.p1_ix], dummy, true);
  
  // Remapping results from 'results_local' to 'result'
  result.val += result_local.val;
  for (uint i = 0; i != per_pts.size(); ++i) 
    result.der[neigh_ixs[per_pts_ixs[i]]] += result_local.der[i];
}
  
// ----------------------------------------------------------------------------
// Energy as a function of distance.  Its support remains within the support of R
array<double, 2> energy(double dist, double R)
// ----------------------------------------------------------------------------
{
  const double tmp = max(R-dist, double(0));
  const double tmp2 = tmp*tmp;
  return {tmp2 * tmp2, -4 * tmp*tmp*tmp}; // energy and derivative
}

}; // end namespace Go
