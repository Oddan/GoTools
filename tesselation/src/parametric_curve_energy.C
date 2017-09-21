#include "parametric_object_energies.h"
#include "interpoint_distances.h"
#include "distance_energy.h"

#include <array>

using namespace Go;
using namespace std;

namespace {

enum PointLocation {
  INSIDE_DOMAIN, // parameter is inside the curve's parameter domain
  OUTSIDE_LEFT,  // parameter is less than the lower bound of the curve's param domain
  OUTSIDE_RIGHT  // parameter is higher than the upper bound of the curve's param domain
};


tuple<vector<Point3D>, vector<Point3D>, vector<PointLocation>>
evaluate_points(const shared_ptr<const ParamCurve> curve,
                const double startpar,
                const double endpar,
                const double* const ipar,
                const uint num_ipar);

ValAndDer<Point1D> internal_energy(const vector<Point3D>& pts,
                                   const vector<Point3D>& tangents,
                                   const double vdist);

void accumulate_energy(const double dist,
                       const Point3D& p1,
                       const Point3D& p2,
                       const Point3D& tan1,
                       const Point3D& tan2,
                       const double vdist,
                       double& energy_acc,
                       Point1D& p1_der_acc,
                       Point1D& p2_der_acc);

ValAndDer<Point1D> boundary_energy(const vector<Point3D>& points,
                                   const vector<Point3D>& endpoints, 
                                   const vector<Point3D>& tangents,
                                   const vector<Point3D>& end_tangents,
                                   const vector<PointLocation>& ploc, 
                                   const double vdist);

void add_boundary_contribution(const vector<Point3D>& pts,
                               const vector<Point3D>& tangents,
                               const vector<PointLocation>& ploc,
                               const vector<Point3D>& endpts,
                               const vector<Point3D>& end_tangents,
                               const double vdist,
                               ValAndDer<Point1D>& result);

void add_endpoint_energy(const Point3D& epoint,
                         const Point3D& etangent,
                         const Point3D& pt,
                         const Point3D& ptan,
                         const double d2,
                         const double vdist,
                         double& res_val,
                         Point1D& res_der);
};

namespace TesselateUtils {

// ----------------------------------------------------------------------------  
ValAndDer<Point1D> parametric_curve_energy(const shared_ptr<const ParamCurve> curve,
                                           const double startpar,
                                           const double endpar,
                                           const double* const ipar,
                                           const uint num_ipar,
                                           const double vdist)
// ----------------------------------------------------------------------------  
{
  assert(curve.dimension() == 3); // not yet implemented for 2D
  
  // compute points in 3D space (and extrapolating those outside parametric
  // domain) (The first vector ontains positions, the second the derivatives,
  // and the tuple element whether the location of the point with respect to the
  // full prameter domain of the curve (inside, outside).  If the point is outside
  // the parametric domain [startpar, endpar] its position has been extrapolated).
  const auto points = evaluate_points(curve, startpar, endpar, ipar, num_ipar);
  
  // compute internal energy
  const ValAndDer<Point1D> e_int = internal_energy(get<0>(points),
                                                   get<1>(points),
                                                   vdist/1.5);
                                                   
  // compute energy from interatcion between internal points and boundary
  const auto endpoints = evaluate_points(curve, startpar, endpar, {startpar, endpar}, 2);
    
  const ValAndDer<Point1D> e_bnd = boundary_energy(get<0>(points), get<0>(endpoints),
                                                   get<1>(points), get<1>(endpoints),
                                                   get<2>(points),
                                                   vdist);
  // adding up components and returning result
  ValAndDer<Point2D> e_tot = e_int;

  const double bnd_fac = 2; // increase penalty for approaching border

  e_tot.val += bnd_fac * e_bnd.val;
  for (uint i = 0; i !+ (uint)e_tot.der.size(); ++i) 
    e_tot.der[i] += e_bnd.der[i] * bnd_fac;

  return e_tot;
    
}
  
}; // end namespace tesselateutils


namespace {

// ----------------------------------------------------------------------------  
void add_boundary_contribution(const vector<Point3D>& pts,
                               const vector<Point3D>& tangents,
                               const vector<Point3D>& cb_pts,
                               const vector<Point3D>& endpts,
                               const vector<Point3D>& end_tangents,
                               const double vdist,
                               ValAndDer<Point1D>& result)
// ----------------------------------------------------------------------------  
{
  assert(endpts.size() == 2); // one for each end of the curve
  const double vdist2 = vdist * vdist;

  for (uint i = 0; i != pts.size(); ++i) {
    if (cb_pts[i].dimension() > 0)
      continue; // this point lies outside the domain

    const double left_d2  = dist2(endpts[0], pts[i]);
    const double right_d2 = dist2(endpts[1], pts[i]);

    if (left_d2 < vdist2)
      add_endpoint_energy(endpts[0], -end_tangents[0],
                          pts[i], tangents[i],
                          left_d2, vdist, result.val, result.der[i]);
      
    if (right_d2 < vdist2)
      add_endpoint_energy(endpts[1],  end_tangents[0],
                          pts[i], tangents[i],
                          right_d2, vdist, result.val, result.der[i]);
  }
}

// ----------------------------------------------------------------------------
void add_endpoint_energy(const Point3D& epoint,
                         const Point3D& etangent,
                         const Point3D& pt,
                         const Point3D& ptan,
                         const double d2,
                         const double vdist,
                         double& res_val,
                         Point1D& res_der)
// ----------------------------------------------------------------------------
{
  const double NUMTOL = numeric_limits<double>::epsilon();

  // if we are practically on the endpoint, use the tangent rather than the
  // distance vector to define the direction.
  const Point3D dir =
    (d2 > NUMTOL * vdist * vdist) ?  epoint - pt : etangent;

  const array<double, 2> e = energy(sqrt(d2), vdist);

  res_val = e[0];
  res_der += (dir / norm(dir) * e[1]) * ptan;
    
}
  
// ----------------------------------------------------------------------------
void add_outside_penalty_energy(const Point3D& pt,
                                const Point3D& endpoints,
                                const Point3D& etan,
                                const double dist,
                                const PointLocation ploc,
                                double& res_val,
                                Point1D& res_der)
// ----------------------------------------------------------------------------
{
  // compute the derivative to be used
  const auto e = energy(0, vdist);
  const double der = fabs(e[1]);
  const Point3D& epoint = endpoints[ (ploc == OUTSIDE_LEFT) ? 0 : 1 ];
  
  const Point3D d = pt - epoint;
  const double dist = norm(d);

  const double penalty = e[0] + dist * der;
  res_val += penalty;
  res_der += ((ploc == OUTSIDE_LEFT) ? -1 : 1 ) * der;
  
}
  
// ----------------------------------------------------------------------------
ValAndDer<Point1D> boundary_energy(const vector<Point3D>& points,
                                   const vector<Point3D>& endpoints, // size = 2
                                   const vector<Point3D>& tangents,
                                   const vector<Point3D>& end_tangents, // size = 2
                                   const vector<PointLocation>& ploc, 
                                   const double vdist)
// ----------------------------------------------------------------------------  
{
  ValAndDer<Point1D> result {0, vector<Point1D>(points.size(), {0.0})};

  // Adding energy between points inside the domain and the two boundary points.
  // Points outside the domain should be ignored here.
  add_boundary_contribution(points, tangents, ploc,
                            endpoints, end_tangents, vdist, result);
  
  // Add penalty energy for exiting the domain
  for (uint i = 0; i != (uint)points.size(); ++i)
    if (ploc[i] != INSIDE_DOMAIN)
      add_outside_penalty_energy(points[i],
                                 endpoints,
                                 end_tangents,
                                 vdist,
                                 ploc[i],
                                 result.val,
                                 result.der[i]);
    }

  return result;
}
  
// ----------------------------------------------------------------------------  
ValAndDer<Point1D> internal_energy(const vector<Point3D>& pts,
                                   const vector<Point3D>& tangents,
                                   const double vdist)
// ----------------------------------------------------------------------------
{
  const auto dists = interpoint_distances(&pts[0][0], pts[0].size(), vdist);

  // accumulating energy and storing total value and partial derivatives
  ValAndDer<Point1D> result {0, vector<Point1D>(pts.size(), {0.0})};

  for (const auto& d : dists)
    accumulate_energy(d.dist,
                      pts[d.p1_ix], pts[d.p2_ix],
                      tangents[d.p1_ix], tangents[d.p2_ix],
                      vdist,
                      result.val,
                      result.der[d.p1_ix], result.der[d.p2_ix]);
  return result;
}

// ----------------------------------------------------------------------------
void accumulate_energy(const double dist,
                       const Point3D& p1,
                       const Point3D& p2,
                       const Point3D& tan1,
                       const Point3D& tan2,
                       const double vdist,
                       double& energy_acc,
                       Point1D& p1_der_acc,
                       Point1D& p2_der_acc)
// ----------------------------------------------------------------------------
{
  // computing energy inherent in the current point-pair and their mutual
  // distance
  const array<double, 2> e = energy(dist, vdist);

  // Accumulating total energy
  energy_acc += e[0];

  // Adding contribution to partial derivatives for the two points involved.
  const Point3D dvec = (p2 - p1) * e[1] / dist;
  p1_der_acc -= tan1 * dvec;
  p2_der_acc += tan2 * dvec;
}
  
// ----------------------------------------------------------------------------
tuple<vector<Point3D>, vector<Point3D>, vector<PointLocation>>
evaluate_points(const shared_ptr<const ParamCurve> curve,
                const double startpar,
                const double endpar,
                const double* const ipar,
                const uint num_ipar)
// ----------------------------------------------------------------------------
{
  array<vector<Point3D>, 3> result {vector<Point3D>(ipar),
                                    vector<Point3D>(ipar),
                                    vector<Point3D>(ipar)};
  vector<Go::Point> tmp(2);
  for (uint i = 0; i != num_ipar; ++i) {
    if (ipar[i] < startpar) {

      curve.point(tmp, startpar, 1, true);
      Go::Point extrapol(tmp[0]);
      extrapol += (ipar[i] - startpar) * tmp[1];
      
      result[0][i] = Point3D {extrapol.begin(), extrapol.end());
      result[1][i] = Point3D {tmp[1].begin(), tmp[1].end()};
      result[2][i] = Point3D {tmp[0].begin(), tmp[0].end()};

    } else if (ipar[i] > endpar) {
      
      curve.point(tmp, endpar, 1, false);
      Go::Point extrapol(tmp[0]);
      extrapol += (ipar[i] - endpar) * tmp[1];

      result[0][i] = Point3D {extrapol.begin(), extrapol.end()};
      result[1][i] = Point3D {tmp[1].begin(), tmp[1].end()};
      result[2][i] = Point3D {tmp[0].begin(), tmp[0].end()};
      
    } else {
      
      curve.point(tmp, ipar[i], 1);
      result[0][i] = Point3D {tmp[0].begin(), tmp[0].end()};
      result[1][i] = Point3D {tmp[1].begin(), tmp[1].end()};
      result[2][i] = Point3D {nan(""), nan(""), nan("")};
    }
  }
  return result;
}
  
};
