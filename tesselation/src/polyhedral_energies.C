#include "polyhedral_energies.h"
#include "interpoint_distances.h"
#include "tesselate_utils.h"
#include "common_defs.h"

using namespace std;
using namespace TesselateUtils;

namespace {

// first entry in result array is the energy, the second the derivative
array<double, 2> energy(double dist, double R);

ValAndDer internal_energy(const Point2D* const ipoints,
			  const unsigned int num_ipoints,
			  const double vdist);
ValAndDer boundary_energy(const Point2D* const bpoints,
			  const unsigned int num_bpoints,
			  const Point2D* const ipoints,
			  const unsigned int num_ipoints,
			  const double vdist);

void add_boundary_contribution(const Point2D& bp1,
			       const Point2D& bp2,
			       const Point2D* const ipoints,
			       const unsigned int num_ipoints,
			       const double vdist,
			       ValAndDer& result);

void accumulate_energy(const double dist,
		       const Point2D& p1,
		       const Point2D& p2,
		       const double vdist,
		       double& energy_acc,
		       Point2D& p1_der_acc,
		       Point2D& p2_der_acc);

  
}; // end anonymous namespace 

namespace TesselateUtils {

// ----------------------------------------------------------------------------
ValAndDer polygon_energy(const Point2D* const bpoints,
			 const unsigned int num_bpoints,
			 const Point2D* const ipoints,
			 const unsigned int num_ipoints,
			 const double vdist)
// ----------------------------------------------------------------------------
{
  
  // compute internal energy (potential energy between internal points)
  const ValAndDer E_int = internal_energy(ipoints, num_ipoints, vdist);

  // compute boundary energy (energy from interaction between internal points
  // and boundary points, and from internal points and their mirror points)
  const ValAndDer E_bnd = boundary_energy(bpoints, num_bpoints,
  					  ipoints, num_ipoints, vdist);

  // Adding up components and returning results
  ValAndDer E_tot = E_int;

  E_tot.val += E_bnd.val; 
  E_tot.der += E_bnd.der;
  
  return E_tot;
}

}; // end namespace Go


namespace {

// ----------------------------------------------------------------------------
ValAndDer internal_energy(const Point2D* const ipoints,
			  const unsigned int num_ipoints,
			  const double vdist)
// ----------------------------------------------------------------------------
{
  const auto dists = interpoint_distances(ipoints, num_ipoints, vdist, false);

  // accumulating energies and storing total value and partial derivatives in 'result'
  ValAndDer result {0, vector<Point2D>(num_ipoints, {0.0, 0.0})};
  for (const auto& d : dists)
    accumulate_energy(d.dist, ipoints[d.p1_ix], ipoints[d.p2_ix], vdist, 
		      result.val, result.der[d.p1_ix], result.der[d.p2_ix]);

  return result;
}

// ----------------------------------------------------------------------------
void accumulate_energy(const double dist,
		       const Point2D& p1,
		       const Point2D& p2,
		       const double vdist,
		       double& energy_acc,
		       Point2D& p1_der_acc,
		       Point2D& p2_der_acc)
// ----------------------------------------------------------------------------  
{
  const array<double,2> e = energy(dist, vdist);

  // accumulating total energy
  energy_acc += e[0];
  
  // adding contribution to partial derivatives for the two points involved.
  // @@ Creation of temporary object on the line below.  If bottleneck, should
  // be rewritten.
  const Point2D dvec = (p2 - p1) * (e[1] / dist);
  p1_der_acc -= dvec;
  p2_der_acc += dvec;
  
}
						 
// ----------------------------------------------------------------------------
ValAndDer boundary_energy(const Point2D* const bpoints,
			  const unsigned int num_bpoints,
			  const Point2D* const ipoints,
			  const unsigned int num_ipoints,
			  const double vdist)
// ----------------------------------------------------------------------------
{
  ValAndDer result {0, vector<Point2D>(num_ipoints, {0.0, 0.0})};

  // looping across boundary segments and adding their energy contributions
  for (uint i = 0; i != num_bpoints; ++i)
    add_boundary_contribution(bpoints[i], bpoints[(i+1)%num_bpoints],
			      ipoints, num_ipoints, vdist, result);
  return result;
}

// ----------------------------------------------------------------------------  
void add_boundary_contribution(const Point2D& bp1,
			       const Point2D& bp2,
			       const Point2D* const ipoints,
			       const unsigned int num_ipoints,
			       const double vdist,
			       ValAndDer& result)
// ----------------------------------------------------------------------------
{
  // identify the points within reach of the boundary ("neighbor points"). 
  const auto npoints = extract_from_range(ipoints, num_ipoints,
  			     [&bp1, &bp2, vdist](const Point2D& p) {
        return point_on_line_segment(p, bp1, bp2, vdist, false);});
  const auto& neigh_ixs    = npoints.second;   
  const auto& neigh_points = npoints.first; 

  // Setting up a result structure only for the neigh poitns
  ValAndDer result_local {0, vector<Point2D>(neigh_points.size(), {0.0, 0.0})};
  Point2D dummy; // used to store values we don't want
  
  // compute the energy (and associated partial derivatives) between the
  // neighbor points and the two boundary points.  
  vector<DistanceEntry> dvec;
  dvec = interpoint_distances(&neigh_points[0], (uint)neigh_points.size(), bp1, vdist);
  for (const DistanceEntry& d : dvec)
    accumulate_energy(d.dist, neigh_points[d.p1_ix], bp1, vdist,
		      result_local.val, result_local.der[d.p1_ix], dummy);
  dvec = interpoint_distances(&neigh_points[0], (uint)neigh_points.size(), bp2, vdist);
  for (const auto& d : dvec)
    accumulate_energy(d.dist, neigh_points[d.p1_ix], bp2, vdist, 
		      result_local.val, result_local.der[d.p1_ix], dummy);

  // remapping results from result_local to result.  Multiplying by 0.5, since
  // segment endpoints will be counted again when neighbour segments are treated
  result.val += result_local.val * 0.5;
  for (uint i = 0; i != neigh_ixs.size(); ++i) 
    result.der[neigh_ixs[i]] += result_local.der[i] * 0.5;
  
  // Now, we compute the energy between the neighbor points and their mirror images.
  // first, remove points that do not project perpendicularly to within the segment.
  const auto per = extract_from_range(&neigh_points[0], (uint)neigh_points.size(),
				      [&bp1, &bp2] (const Point2D& p) {
					return projects_to_segment(p, bp1, bp2);});
  const auto& per_pts_ixs = per.second;
  const auto& per_pts     = per.first;
  result_local.reset((uint)per_pts.size()); // reinitialize local result structure

  const auto mpoints = mirror_points_2D(&per_pts[0], (uint)per_pts.size(), bp1, bp2);
  dvec = interpoint_distances(&per_pts[0], (uint)per_pts.size(),
			      &mpoints[0], (uint)per_pts.size(), vdist);

  for (const auto& d : dvec)
    accumulate_energy(d.dist, per_pts[d.p1_ix], mpoints[d.p2_ix], vdist,
		      result_local.val, result_local.der[d.p1_ix], dummy);

  // remapping results from result_local to result
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
