#include "parametric_object_energies.h"

#include <array>

using namespace Go;
using namespace std;

namespace {
  const vector<array<Point, 3>> evaluate_points(const shared_ptr<const ParamCurve> curve,
                                                const double startpar,
                                                const double endpar,
                                                const double* const ipar,
                                                const uint num_ipar);
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
  // compute points in 3D space (and extrapolating those outside parametric
  // domain) (The first point in the array is position, the second the
  // derivative, and the third the closest point on the boundary (only defined
  // if the point is outside the parametric domain [startpar, endpar] and has
  // thus been extrapolated).
  const vector<array<Point, 3>> evaluate_points(curve, startpar,
                                                endpar, ipar, num_ipar);
  
  // compute internal energy
  const ValAndDer<Point1D> e_int = internal_energy(curve, ipar, num_ipar, vdist);

  // compute energy from interatcion between internal points and boundary
  const ValAndDer<Point1D> e_bnd = boundary_energy(curve, ipar, num_ipar, vdist);

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
const vector<array<point, 3>> evaluate_points(const shared_ptr<const paramcurve> curve,
                                              const double startpar,
                                              const double endpar,
                                              const double* const ipar,
                                              const uint num_ipar)
// ----------------------------------------------------------------------------
{
  vector<array<point, 3>> result(num_ipar);
  vector<Point> tmp(2);
  for (uint i = 0; i != num_ipar; ++i) {
    if (ipar[i] < startpar) {

      curve.point(tmp, startpar, 1, true);
      Point extrapol(tmp[0]);
      extrapol += (ipar[i] - startpar) * tmp[1];
      
      result[i] = {extrapol, tmp[1], tmp[0]};
      
    } else if (ipar[i] > endpar) {
      
      curve.point(tmp, endpar, 1, false);
      Point extrapol(tmp[0]);
      extrapol += (ipar[i] - endpar) * tmp[1];
      
      result[i] = {extrapol, tmp[1], tmp[0]};
      
    } else {
      
      curve.point(tmp, ipar[i], 1);
      result[i] = {tmp[0], tmp[1], Point()};
      
    }
  }
  return result;
}
  
};
