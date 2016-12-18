#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <list>

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"

#include "GoTools/lrsplines2D/LRSplineSurface.h"

using namespace Go;
using namespace std;

namespace {
  enum PointStatus {REGULAR = 0, BOUNDARY=1, CYCLIC_END = 3, SINGULAR = 4};
  enum PointIterationOutcome {OK = 0, ITERATION_EXCEED = 1, SINGULARITY_ENCOUNTERED = 2};
  using CurvePtr = shared_ptr<const SplineCurve>;
  using Array4   = array<double, 4>;

  CurvePtr traceIsoval(const SplineSurface& surf, double u, double v, double tol);
  Array4 get_derivs(const SplineSurface& surf, const Point& p, PointStatus& status);

  PointIterationOutcome move_point_to_isocontour(const SplineSurface& surf,
						 const double isoval, const double tol,
						 Point& uv);
  void save_sampled_surface(const SplineSurface& surf, const int samples, string filename);
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

  // // sample and save surface, for debug purposes
  // save_sampled_surface(surf, 200, "sampled.mat");
  // return 0;
  
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
  const double isoval = p1[0];
  cout << "Point: " << p1 << endl;
  Point ipoint = {u+dt, v};
  surf.point(p1, ipoint[0], ipoint[1]);
  cout << "Point (u+dt): " << p1 << endl;
  bool success = (move_point_to_isocontour(surf, isoval, 1e-12, ipoint) == OK);
  cout << "Convergence: " << (success ? " Yes!" : " No!") << endl;
  surf.point(p1, ipoint[0], ipoint[1]);
  cout << "Point (parameter line): " << p1 << " " <<  ipoint << endl;
  
  Point tpoint = {u+derivs[0] * dt, v + derivs[1] * dt};
  surf.point(p1, tpoint[0], tpoint[1] );
  cout << "Point (tangent): " << p1 << endl;
  success = (move_point_to_isocontour(surf, isoval, 1e-12, tpoint) == OK);
  cout << "Convergence: " << (success ? " Yes!" : " No!") << endl;
  surf.point(p1, tpoint[0], tpoint[1]);
  cout << "Point (tangent): " << p1 << " " << tpoint << endl;

  Point cpoint = {u + derivs[0] * dt + 0.5 * derivs[2] * dt * dt,
		  v + derivs[1] * dt + 0.5 * derivs[3] * dt * dt};
  surf.point(p1, cpoint[0], cpoint[1]);
  cout << "Point (curve): " << p1 << endl;

  success = (move_point_to_isocontour(surf, isoval, 1e-12, cpoint) == OK);
  cout << "Convergence: " << (success ? " Yes!" : " No!") << endl;
  surf.point(p1, cpoint[0], cpoint[1]);
  cout << "Point (iterated): " << p1 << " " << cpoint << endl;
  
  const CurvePtr result = traceIsoval(surf, u, v, 1e-10);

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
double value_at(const SplineSurface& s, const Point& par)
// ----------------------------------------------------------------------------
{
  static Point tmp;
  s.point(tmp, par[0], par[1]);
  return tmp[0];
}
  
// ----------------------------------------------------------------------------
inline Point tangent_at(const SplineSurface& surf, const Point& p)
// ----------------------------------------------------------------------------
{
  static vector<Point> tmp(3);
  surf.point(tmp, p[0], p[1], 1);
  Point res {tmp[2][0], -tmp[1][0]};
  res = res / res.length(); // @@ What to do in case of singularity?
  return res;
}

  
// ----------------------------------------------------------------------------
void save_sampled_surface(const SplineSurface& surf, const int samples, string filename)
// ----------------------------------------------------------------------------
{
  ofstream os(filename.c_str());

  const double umin = surf.startparam_u();
  const double vmin = surf.startparam_v();
  const double du   = (surf.endparam_u() - umin) / (samples-1);
  const double dv   = (surf.endparam_v() - vmin) / (samples-1);
  
  for (int i = 0; i < samples; ++i)
    for (int j = 0; j < samples; ++j)
      os << value_at(surf, Point { umin + i * du, vmin + j * dv}) << " ";
		       
  os.close();
}
  
  
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

// ----------------------------------------------------------------------------    
inline void truncate_derivs_to_domain(const SplineSurface& surf, const Point& uv,
				      double& du, double& dv)
// ----------------------------------------------------------------------------
{
  du = ((uv[0] > surf.startparam_u()) & (uv[0] < surf.endparam_u())) ? du : 0;
  dv = ((uv[1] > surf.startparam_v()) & (uv[1] < surf.endparam_v())) ? dv : 0;
}

// ----------------------------------------------------------------------------  
PointIterationOutcome move_point_to_isocontour(const SplineSurface& surf,
					       const double isoval,
					       const double tol, Point& uv)
// ----------------------------------------------------------------------------  
{
  const double SING_TOL = 1e-9; // @@ passed as parameter?
  const int MAX_ITER = 10;
  static vector<Point> cur_val(3);
  surf.point(cur_val, uv[0], uv[1], 1);

  truncate_derivs_to_domain(surf, uv, cur_val[1][0], cur_val[2][0]);
  
  const double& s     = cur_val[0][0];
  const double& s_du  = cur_val[1][0];
  const double& s_dv  = cur_val[2][0];
  
  double grad2 = s_du * s_du + s_dv * s_dv;
  if (grad2 < SING_TOL)
    // vanishing gradient.  Singularity encountered.
    return SINGULARITY_ENCOUNTERED;

  double err = s - isoval;
  int iter = 0;

  while (pow(err, 2) > (tol * tol * grad2)) {

    const double grad2_inv = double(1)/grad2;
    
    uv[0] = uv[0] - err * s_du * grad2_inv;
    uv[1] = uv[1] - err * s_dv * grad2_inv;

    surf.point(cur_val, uv[0], uv[1], 1);
    truncate_derivs_to_domain(surf, uv, cur_val[1][0], cur_val[2][0]);
    
    err = s - isoval;
    grad2 = s_du * s_du + s_dv * s_dv;

    if (++iter > MAX_ITER)
      // did not converge within the allowed number of iterations
      return ITERATION_EXCEED;
  }
  // converged
  cout << "Iterations: " << iter << endl;
  return OK;
}

// ----------------------------------------------------------------------------  
Array4 get_derivs(const SplineSurface& surf, const Point& p, PointStatus& status)
// ----------------------------------------------------------------------------
{
  // If the parameter domain is described by (u, v) and the arc length
  // parameterization of the curve represented by 't', then the entries of the
  // returned array will be: [du/dt, dv/dt, d2u/dt2, d2v/dt2].
  static vector<Point> tmp(6, {0.0, 0.0});
  const double SING_TOL = 1e-7; // @@ is this value reasonable?  Should it be a
				// parameter to be passed
    
  surf.point(tmp, p[0], p[1], 2);  // evaluate surface and its first and second
				   // derivatives
  const double ds_du    = tmp[1][0];
  const double ds_dv    = tmp[2][0];
  const double d2s_du2  = tmp[3][0];
  const double d2s_dudv = tmp[4][0];
  const double d2s_dv2  = tmp[5][0];

  const double dsdu2 = pow(ds_du, 2);
  const double dsdv2 = pow(ds_dv, 2);

  const double n = sqrt(dsdu2 + dsdv2);
  const double a = 1.0 / n;
  const double da_dt = pow(a, 4) *
    (d2s_dudv * (dsdu2 - dsdv2) + ds_du * ds_dv * (d2s_dv2 - d2s_du2));

  // determininng point status
  status =
    (n < SING_TOL)                                                    ? SINGULAR :
    (((p[0] <= surf.startparam_u()) | (p[0] >= surf.endparam_u())) ||
    ((p[1] <= surf.startparam_v()) | (p[1] >= surf.endparam_v())))    ? BOUNDARY :
                                                                        REGULAR;

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
  return dt;
}

// ----------------------------------------------------------------------------
Point extrapolate_point(const SplineSurface& surf, double dt, const Point& cur_point,
			const Array4& derivs, const double isoval, double tol, PointStatus& status)
// ----------------------------------------------------------------------------
{
  // determine a new test point based on previous point and derivatives.
  // Converge the new point so that it lies on the isocurve.  Check status of
  // the new point (could be REGULAR, BOUNDARY or SINGULAR).

  // Loop until a suitable point is found
  bool found = false;
  Point new_point;
  while (!found) {
    // First, use derivatives to estimate position of new point
    new_point = cur_point;
    for (int i = 0; i != 2; ++i)
      new_point[i] += dt * derivs[i] + 0.5 * dt * dt * derivs[i+2];
    
    // move estimated point to an actual isocontour position
    switch (move_point_to_isocontour(surf, isoval, tol, new_point)) {
    case OK:
      status = (new_point[0] <= surf.startparam_u() ||
		new_point[0] >= surf.endparam_u()   ||
		new_point[1] <= surf.startparam_v() ||
		new_point[1] >= surf.endparam_v()) ? BOUNDARY : REGULAR;
      found = true;
      break;
    case ITERATION_EXCEED:
      // Could not converge to new point.  Shorten steplength and try again
      dt = dt/2;
      break;
    case SINGULARITY_ENCOUNTERED:
      // check if singularity is on the curve.  If not, shorten timestep and try again
      const double err = fabs(isoval - value_at(surf, new_point));
      if (err < tol) {
	found = true;
	status = SINGULAR;
      } else {
	dt = dt/2;
      }
      break;
    }
  }
  // If we have not been able to converge to a new point, flag it as singular.
  // @@ Validity of this should be tested
  if (cur_point.dist2(new_point) < tol * tol)
    status = SINGULAR;
  return new_point;
}

// ----------------------------------------------------------------------------
double estimate_arclength(const Point& p1, const Point& p2,
			  const Point& t1, const Point& t2)
// ----------------------------------------------------------------------------
{
  const double EPS_ANGLE = 1e-5;
  const double theta = t1.angle(t2);
  const double d = p1.dist(p2);

  if (theta < EPS_ANGLE)
    // assume linear
    return d;

  // Radius of circle segment
  const double R = 0.5 * d / sin(0.5 * theta);

  return R * theta;
}  

// ----------------------------------------------------------------------------
bool check_midpoint(const SplineSurface& surf, const Point& pstart, const Point& pend,
		    const bool forward, const double isoval, const double tol,
		    Point& midpoint, PointIterationOutcome& outcome)
// ----------------------------------------------------------------------------
{
  // Function returns 'true' if midpoint is already OK
  outcome = OK; // this is the default case
  
  // Computing tangents @@ could be passed along, for optimization purposes...
  const Point der_start = tangent_at(surf, pstart) * (forward ? 1 : -1);
  const Point der_end   = tangent_at(surf, pend) * (forward ? 1 : -1);

  const double arclen = estimate_arclength(pstart, pend, der_start, der_end);
  
  // Estimating position of new point (assuming a cubic hermite spline with arc
  // length parameterization, tangent scaled by approximate arc length)
  midpoint = 0.5 * (pstart + pend) + 0.125 * (der_start - der_end) * arclen;

  // Check whether the new point is close enough to the isoval
  vector<Point> tmp(3); surf.point(tmp, midpoint[0], midpoint[1], 1);
  const double err = tmp[0][0] - isoval;
  const double grad2 = (tmp[1][0] * tmp[1][0]) + (tmp[2][0] * tmp[2][0]);
  if (pow(err, 2) <= (tol * tol * grad2))
    return true;

  outcome = move_point_to_isocontour(surf, isoval, tol, midpoint);
  return false;
}

// ----------------------------------------------------------------------------
vector<Point> intermediate_points(const SplineSurface& surf, const Point& old_point,
				  const Point& new_point, const bool forward,
				  const double isoval, const double tol, bool& ok)
// ----------------------------------------------------------------------------
{
  // Find enough intersection points between 'new_point', and 'old_point', to
  // ensure that the hermite interpolation will remain sufficiently close to the
  // real intersection curve at any time.

  // First, check if a point needs to be inserted at all.  If true, we will
  // establish a list of new points to insert

  ok = true;
  Point midpoint;
  PointIterationOutcome outcome;
  if (!check_midpoint(surf, old_point, new_point, forward, isoval, tol, midpoint, outcome) &
      (outcome == OK)) {
    
    // this curve needs to be refined.  But more than one point might be required.
    list<Point> tmp {old_point, midpoint, new_point};
    for (auto it = tmp.begin(), it_prev = it++; it != tmp.end(); it_prev = it++) {

      if (!check_midpoint(surf, *it_prev, *it, forward, isoval, tol, midpoint, outcome) &
	  (outcome == OK)) {
	tmp.insert(it, midpoint);
	--it; // the new interval has to be checked, so we decrement the iterator again
      }
      if (outcome != OK)
	break;
    }

    ok = (outcome == OK);
    return vector<Point>(++tmp.begin(), --tmp.end()); // Skip first and last point
  } 
  
  // if we got here, no midpoint was needed, and we can return an empty vector
  ok = (outcome == OK);
  return vector<Point>();

}

  

// ----------------------------------------------------------------------------
inline double signed_angle(const Point& dir1, const Point& dir2)
// ----------------------------------------------------------------------------
{
  double res = dir1.angle2(dir2);
  return (res < M_PI) ? res : 2 * M_PI - res;
}

// ----------------------------------------------------------------------------
inline bool in_between(const SplineSurface& surf, const Point& testp,
		       const Point& startp, const Point& endp)
// ----------------------------------------------------------------------------
{
  const double start_end_dist2 = startp.dist2(endp);
  const double dist2_1 = testp.dist2(startp);
  const double dist2_2 = testp.dist2(endp);

  // First, check whether the testpoint is "in-between" the start and end point,
  // distance-wise.
  if (dist2_1 > start_end_dist2 || dist2_2 > start_end_dist2)
    return false;

  // compare directions of tangents to intersection curve
  const Point start_dir = tangent_at(surf, startp);
  const Point end_dir   = tangent_at(surf, endp);
  const Point test_dir  = tangent_at(surf, testp);
  const double angle_diff_1 = signed_angle(start_dir, test_dir);
  const double angle_diff_2 = signed_angle(test_dir, end_dir);

  // check the continuity of directional change
  if (angle_diff_1 * angle_diff_2 < 0)
    return false;

  // angle and distance are consistent with testpoint being on the intersetion
  // curve between the start and end points.  @@ This test may be too lax - check!
  return true;
}

// ----------------------------------------------------------------------------
int closed_cycle_check(const SplineSurface& surf, const vector<Point>& points,
		       const size_t old_size)
// ----------------------------------------------------------------------------
{
  // check if the first point of the vector (start point) is found in-between
  // the newly inserted points, in which case we assume the cycle to have been
  // closed.
  const Point& startpoint = points.front();

  // Initial check to see if we have moved far enough along the curve to be in
  // the vicinity of the curve's start point.
  if (startpoint.dist2(points.back()) > points[old_size-1].dist2(points.back()))
    return 0;
  
  for (size_t it = old_size; it != points.size(); ++it) 
    if (in_between(surf, startpoint, points[it-1], points[it])) 
      return int(it-1);

  // no end-of-cycle found
  return 0; 
}
				  
// ----------------------------------------------------------------------------
PointStatus find_next_point(const SplineSurface& surf, vector<Point>& prev_points,
			    bool forward, double isoval, Array4& derivs, double tol)
// ----------------------------------------------------------------------------
{
  PointStatus status = REGULAR; // Will be changed in functions call below
  // upon entry, 'derivs' should contain the derivatives corresponding to the
  // last point in the 'prev_points' vector.  At exit, 'derivs' will contain the
  // derivatives corresponding to the last point added.  The use of this
  // variable is for efficiency only.
  const double dt = choose_steplength(derivs);

  // returned status here can be REGULAR, BOUNDARY or SINGULAR (CYCLIC_END will
  // be checked for later).  Since the current point is REGULAR, we should be
  // guaranteed that another point can be produced.
  const Point new_pt = extrapolate_point(surf, (forward ? dt : -dt),
					 prev_points.back(), derivs, isoval, tol, status);

  // ensure that the curve between previous and new point will stick closely
  // enough to the actual intersection
  bool ok; // if a topology problem is detected (points lying on different
	   // curves), ok will be set to 'false' upon return of the function
	   // call below. In that case, the returned vector will contain the
	   // points assumed to be lying on the present curve (but the new point
	   // will not, and should thus be discarded).
  const size_t prev_size = prev_points.size();
  const vector<Point> ipoints = intermediate_points(surf, prev_points.back(),
						    new_pt, forward, isoval, tol, ok);

  prev_points.insert(prev_points.end(), ipoints.begin(), ipoints.end());
  if (ok)
    prev_points.emplace_back(new_pt);
		     
  // check for possibly closed cycle
  const int ix = closed_cycle_check(surf, prev_points, prev_size);
  if (ix > 0) {
    status = CYCLIC_END;
    assert(ix >= int(prev_size)-1);
    prev_points.erase(prev_points.begin() + ix, prev_points.end());
  }

  PointStatus ps_tmp;
  derivs = get_derivs(surf, prev_points.back(), ps_tmp);
  if (status != CYCLIC_END)
    status = ps_tmp;
  
  return status;
}

// ----------------------------------------------------------------------------
vector<Point> trace_unidir(const SplineSurface& surf, const Point& startpoint,
			   const bool forward, double tol, PointStatus& last_point_status)
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
    last_point_status = find_next_point(surf, result, forward, isoval, derivs, tol);

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
CurvePtr traceIsoval(const SplineSurface& surf, double u, double v, double tol)
// ============================================================================
{
  assert(surf.dimension()==1);
  PointStatus ps;

  // trace in the first direction
  const auto res1 = trace_unidir(surf, {u, v}, true, tol, ps);

  // if we did not get back to where we began, trace in second direction also
  const auto res =
    (ps == CYCLIC_END) ?
    res1 : merge_vec(reverse_vec(trace_unidir(surf, {u, v}, false, tol, ps)), res1);


  // debug
  cout << "Number of points: " << res1.size() << endl;
  ofstream os("krull.mat");
  for (auto p : res1)
    os << p[0] <<  " " << p[1] << '\n';

  os.close();

  return curve_from_points(surf, res);
  
}

} // end anonymous namespace

