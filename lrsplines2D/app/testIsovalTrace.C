#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"

#include "GoTools/lrsplines2D/LRSplineSurface.h"

using namespace Go;
using namespace std;

namespace {
  enum PointStatus {REGULAR = 0, BOUNDARY=1, CYCLIC_END = 3, SINGULAR = 4};
  using CurvePtr = shared_ptr<const SplineCurve>;
  using Array4   = array<double, 4>;

  CurvePtr traceIsoval(const SplineSurface& surf, double u, double v);

  Array4 get_derivs(const SplineSurface& surf, const Point& p, PointStatus& status);
  
}; // end anonymous namespace 


// ============================================================================
int main(int varnum, char* vararg[])
// ============================================================================
{
  const double kvec[] = {0, 0, 0, 0, 1, 1, 1, 1};
  const double coefs[] = {0, 0, 0, 0,
			  0, 0, 0, 0,
			  0, 0, 1, 0,
			  0, 0, 0, 0};

  const SplineSurface surf(4, 4, 4, 4, kvec, kvec, coefs, 1, false);

  
  const double u = 0.5;
  const double v = 0.5;

  Point p1;
  PointStatus status;
  const Array4 derivs = get_derivs(surf, {u, v}, status);
  const double  k = 3;
  const double tlim_1 = fabs(2 * derivs[0] / (k * derivs[2]));
  const double tlim_2 = fabs(2 * derivs[1] / (k * derivs[3]));
  //const double dt = 0.15;
  const double dt = min(tlim_1, tlim_2);
  cout << "Chosen steplength: " << dt << endl  << endl;

  const double K = 3;         // @@ this factor might be adapted
  const double d1 = max(fabs(derivs[0]), fabs(derivs[1]));
  const double d2 = max(fabs(derivs[2]), fabs(derivs[3]));

  const double dt_bis = 2 * d1 / (K * d2);

  cout << "Chosen steplength (bis) : " << dt_bis << endl;
  
  surf.point(p1, u, v);
  cout << "Point: " << p1 << endl;
  surf.point(p1, u+dt, v);
  cout << "Point (u+dt): " << p1 << endl;
  surf.point(p1, u+derivs[0] * dt, v + derivs[1] * dt);
  cout << "Point (tangent): " << p1 << endl;
  surf.point(p1,
	     u + derivs[0] * dt + 0.5 * derivs[2] * dt * dt,
	     v + derivs[1] * dt + 0.5 * derivs[3] * dt * dt);
  cout << "Point (curve): " << p1 << endl;
  
  //const CurvePtr result = traceIsoval(surf, u, v);


  // plot 3D version of surface (for use with goview only).  LR surface is here
  // used as goview supports plotting of 1D LR surfaces, but not regular
  // SplineSurfaces.
  const LRSplineSurface lrs(&surf, double(1e-6));
  ofstream os("dill.g2");
  lrs.writeStandardHeader(os);
  lrs.write(os);
  os.close();
  
  return 0;
};


namespace {

// ----------------------------------------------------------------------------
template<typename T>
vector<T> merge_vec(vector<T> v1, const vector<T>& v2)
// ----------------------------------------------------------------------------
{
  v1.insert(v1.end(), v2.begin(), v2.end()); return v1;
}

// ----------------------------------------------------------------------------
template<typename T>
vector<T> reverse_vec(vector<T> v)
// ----------------------------------------------------------------------------
{
  reverse(v.begin(), v.end()); return v;
}

// // ----------------------------------------------------------------------------  
// bool check_for_closed_cycle(const Point& p, const vector<Point>& vec)
// // ----------------------------------------------------------------------------  
// {
//   return true; // @@ dummy for now
// }

// ----------------------------------------------------------------------------  
Array4 get_derivs(const SplineSurface& surf, const Point& p, PointStatus& status)
// ----------------------------------------------------------------------------
{
  // If the parameter domain is described by (u, v) and the arc length
  // parameterization of the curve represented by 't', then the entries of the
  // returned array will be: [du/dt, dv/dt, d2u/dt2, d2v/dt2].
  static vector<Point> tmp(6, {0.0, 0.0});
  surf.point(tmp, p[0], p[1], 2);  // evaluate surface and its first and second
				   // derivatives
  const double ds_du    = tmp[1][0];
  const double ds_dv    = tmp[2][0];
  const double d2s_du2  = tmp[3][0];
  const double d2s_dudv = tmp[4][0];
  const double d2s_dv2  = tmp[5][0];

  const double dsdu2 = pow(ds_du, 2);
  const double dsdv2 = pow(ds_dv, 2);

  const double a = 1.0 / sqrt(dsdu2 + dsdv2);
  const double da_dt = pow(a, 4) *
    (d2s_dudv * (dsdu2 - dsdv2) + ds_du * ds_dv * (d2s_dv2 - d2s_du2));

  return Array4
    { a * ds_dv,  // u-component of tangent
     -a * ds_du, // v-component of tangent
      da_dt * ds_dv + pow(a,2) * (d2s_dudv * ds_dv - d2s_dv2 * ds_du), // double deriv.
     -da_dt * ds_du - pow(a,2) * (d2s_du2 * ds_dv - d2s_dudv * ds_du) // ditto
    };
}

// ----------------------------------------------------------------------------
double choose_steplength(const Array4& derivs)
// ----------------------------------------------------------------------------
{
  const double K = 3;         // @@ this factor might be adapted
  const double maxlen = 0.2;  // @@ this factor should be adapted to patch size

  const double d1 = max(fabs(derivs[0]), fabs(derivs[1]));
  const double d2 = max(fabs(derivs[2]), fabs(derivs[3]));

  const double dt = min(maxlen, 2 * d1 / (K * d2));

}
  
// ----------------------------------------------------------------------------
Point find_next_point(const SplineSurface& surf, const vector<Point>& prev_points,
		      bool forward, double isoval, Array4& derivs, PointStatus& status)
// ----------------------------------------------------------------------------
{
  // upon entry, 'derivs' should contain the derivatives corresponding to the
  // last point in the 'prev_points' vector.  At exit, 'derivs' will contain the
  // derivatives corresponding to the last point added.  The use of this
  // variable is for efficiency only.
  const double dt = choose_steplength(derivs);
  return Point(); //@@ Dummy
}

// ----------------------------------------------------------------------------
double value_at(const SplineSurface& s, const Point& par)
// ----------------------------------------------------------------------------
{
  static Point tmp;
  s.point(tmp, par[0], par[1]);
  return tmp[0];
}
  
// ----------------------------------------------------------------------------
vector<Point> trace_unidir(const SplineSurface& surf, const Point& startpoint,
			   const bool forward, PointStatus& last_point_status)
// ----------------------------------------------------------------------------
{
  vector<Point> result(1, startpoint);
  const double isoval = value_at(surf, startpoint); // the isovalue at which the
						    // curve should lie
  
  // get curve tangent, second derivatives, and check the status of the
  // startpoint (REGULAR, BOUNDARY, etc.)
  auto derivs = get_derivs(surf, startpoint, last_point_status); 
  
  // we are at a regular point inside the domain, keep traceing
  while (last_point_status == REGULAR)
    result.emplace_back(find_next_point(surf, result, forward, isoval,
					derivs, last_point_status));

  // only include the startpoint itself if we have been tracing forward
  if (!forward)
    result.erase(result.begin()); 
				  
  return result;
}

// ----------------------------------------------------------------------------
CurvePtr curve_from_points(const SplineSurface& surf, const vector<Point>& pvec)
// ----------------------------------------------------------------------------  
{
  // make a hermite spline based on these points
  const CurvePtr dummy; //@@
  return dummy;
}

// ============================================================================
CurvePtr traceIsoval(const SplineSurface& surf, double u, double v)
// ============================================================================
{
  assert(surf.dimension()==1);
  PointStatus ps;

  // trace in the first direction
  const auto res1 = trace_unidir(surf, {u, v}, true, ps);

  // if we did not get back to where we began, trace in second direction also
  const auto res =
    (ps == CYCLIC_END) ?
    res1 : merge_vec(reverse_vec(trace_unidir(surf, {u, v}, false, ps)), res1);

  return curve_from_points(surf, res);
  
}

} // end anonymous namespace

