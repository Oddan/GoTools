#include <cmath>
#include <list>

#include "sislP.h"
#include "GoTools/geometry/SISLconversion.h"

#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/SSurfTraceIsocontours.h"

using namespace std;
using namespace Go;

namespace {
enum PointStatus {REGULAR = 0, BOUNDARY=1, CYCLIC_END = 3, SINGULAR = 4};

// OK: point could be iterated to convergence, ITERATION_EXCEED: could not
// converge, FLAT_REGION_ENCOUNTERED: Either we are at a flat internal point, or
// at the boundary with a "flat" tangent vector along the boundary.
enum PointIterationOutcome {OK = 0, ITERATION_EXCEED = 1, FLAT_REGION_ENCOUNTERED = 2};
using Array4   = array<double, 4>;
using PandDer  = pair<Point, Array4>; // point with its derivatives

CurveVec compute_levelset(const SplineSurface& ss,
			  SISLSurf* ss_sisl, // SISL surface corresponding to ss
			  SISLSurf* ss_sisl_3D, // only needed if use_sisl_marching = true
			  const double isoval,
			  const double tol, 
			  const bool include_3D_curves,
			  const bool use_sisl_marching);

// ----------------------------------------------------------------------------
inline double value_at(const SplineSurface& s, const Point& par)
// ----------------------------------------------------------------------------
{
  static Point tmp;
  s.point(tmp, par[0], par[1]);
  return tmp[0];
}

// ----------------------------------------------------------------------------
SISLSurf* make_sisl_3D(const SplineSurface& ss)
// ----------------------------------------------------------------------------
{
  // convert ss to 3D and produce a SISLSurf object from it.  Since the
  // conversion to 3D is already defined for LRSplineSurface, we re-use that
  // functionality here, by temporarily converting to LRSplineSurface.  @@ This
  // could easily be changed later by reimplementing the 3D conversion directly
  // here.
  
  assert(ss.dimension() == 1);
  const double knot_tol = 1e-8;
  LRSplineSurface tmp(&ss, knot_tol);
  tmp.to3D(); // 3D conversion takes place here

  return GoSurf2SISL(*shared_ptr<SplineSurface>(tmp.asSplineSurface()), true);
}  

  
}; // end anonymous namespace 

namespace Go {


// ============================================================================
vector<CurveVec> SSurfTraceIsocontours(const SplineSurface& ss,
				       const vector<double>& isovals,
				       const double tol, 
				       bool include_3D_curves,
				       bool use_sisl_marching)
// ============================================================================
{
  assert(ss.dimension() == 1); // only intended to work for spline functions
  
  // Compute topology for each requested level-set.  We use SISL for this
  SISLSurf* sislsurf1D = GoSurf2SISL(ss, false);
  SISLSurf* sislsurf3D = use_sisl_marching ? make_sisl_3D(ss) : nullptr;
  
  // Defining function tracing out the level set for a specified isovalue
  const function<CurveVec(double)> comp_lset = [&] (double ival)
    {return compute_levelset(ss, sislsurf1D, sislsurf3D, ival, tol,
			     include_3D_curves, use_sisl_marching);};

  // Computing all level-set curves for all isovalues ("transforming" each
  // isovalue into its corresponding level-set)
  const vector<CurveVec> result = apply_transform(isovals, comp_lset);

  // Cleaning up after use of SISL objects
  freeSurf(sislsurf1D);
  if (sislsurf3D)
    freeSurf(sislsurf3D);

  // Returning result
  return result;
}

} // end namespace Go;


namespace {

// ----------------------------------------------------------------------------  
pair<SISLIntcurve**, int> compute_topology(SISLSurf* ss_sisl, double isoval)
// ----------------------------------------------------------------------------
{
  // This function does the equivalent of SISL sh1851, but assumes that the
  //provided surface is 1D (bivariate spline function), noncyclic and k-regular.

  
  // Finding intersections, using SISL function sh1761.
  int jpt = 0; // number of single intersection point
  int jcrv = 0; // number of intersection curves
  int jsurf = 0; 
  double* gpar = SISL_NULL; // parameter values of the single intersection points
  double* spar = SISL_NULL; // dummy array
  int *pretop=SISL_NULL;
  SISLIntcurve** wcurve; // array containing descriptions of the intersection curves
  SISLIntsurf** wsurf; 

  auto qp  = newPoint(&isoval, 1, 1);
  auto qo1 = newObject(SISLSURFACE);
  qo1->s1 = ss_sisl;
  qo1->o1 = qo1;
  auto qo2 = newObject(SISLPOINT);
  qo2->p1 = qp;
  
  SISLIntdat* qintdat = SISL_NULL; // intersection result
  auto freeall = [&] (bool skip_wcurves = false) {
    if (gpar)                         free(gpar);
    if (spar)                         free(spar);
    if (pretop)                       free(pretop);
    if (qintdat)                      freeIntdat(qintdat);
    if ((bool)wcurve && (jcrv > 0) && !skip_wcurves) freeIntcrvlist(wcurve, jcrv);
    for (int i = 0; i < jsurf; ++i)   freeIntsurf(wsurf[i]);
    if ((bool)wsurf && (jsurf > 0))          free(wsurf);
  };
  
  auto cleanup_and_throw = [&freeall] (string str) {
    freeall();
    throw runtime_error(str.c_str());
  };
  
  const double epsge = 1e-6;
  int kstat = 0;
  
  // find intersections
  sh1761(qo1, qo2, epsge, &qintdat, &kstat);
  if (kstat < 0) cleanup_and_throw("SISL error in sh1761.");

  const int kdeg = 1;
    
  // express intersections on output format (SISL)
  if (qintdat)
    hp_s1880 (qo1, qo2, kdeg, 2, 0, qintdat, &jpt, &gpar, &spar, &pretop,
	      &jcrv, &wcurve, &jsurf, &wsurf, &kstat);
  
  freeall(true);
  
  return pair<SISLIntcurve**, int>(wcurve, jcrv);
}

// ----------------------------------------------------------------------------    
inline void truncate_derivs_to_domain(const SplineSurface& surf, const Point& uv,
				      double& du, double& dv)
// ----------------------------------------------------------------------------
{
  du = ((uv[0] > surf.startparam_u()) && (uv[0] < surf.endparam_u())) ? du : 0;
  dv = ((uv[1] > surf.startparam_v()) && (uv[1] < surf.endparam_v())) ? dv : 0;
}

// ----------------------------------------------------------------------------
void truncate_point_to_domain(const SplineSurface& surf, Point& uv)
// ----------------------------------------------------------------------------
{
  uv[0] = uv[0] < surf.startparam_u() ? surf.startparam_u() :
	  uv[0] > surf.endparam_u()   ? surf.endparam_u() : uv[0];
  uv[1] = uv[1] < surf.startparam_v() ? surf.startparam_v() :
	  uv[1] > surf.endparam_v()   ? surf.endparam_v() : uv[1];
}


// ----------------------------------------------------------------------------  
PointIterationOutcome move_point_to_isocontour(const SplineSurface& surf,
					       const double isoval,
					       const double tol, Point& uv)
// ----------------------------------------------------------------------------  
{
  const double SING_TOL = 1e-7; // @@ passed as parameter?
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
    return FLAT_REGION_ENCOUNTERED;

  double err = s - isoval;
  int iter = 0;

  while (pow(err, 2) > (tol * tol * grad2)) {

    const double grad2_inv = double(1)/grad2;
    
    uv[0] = uv[0] - err * s_du * grad2_inv;
    uv[1] = uv[1] - err * s_dv * grad2_inv;
    truncate_point_to_domain(surf, uv);

    surf.point(cur_val, uv[0], uv[1], 1);
    truncate_derivs_to_domain(surf, uv, cur_val[1][0], cur_val[2][0]);
    
    err = s - isoval;
    grad2 = s_du * s_du + s_dv * s_dv;

    if (grad2 < SING_TOL)
      // vanishing gradient.  Singularity encountered.
      return FLAT_REGION_ENCOUNTERED;

    if (++iter > MAX_ITER)
      // did not converge within the allowed number of iterations
      return ITERATION_EXCEED;
  }
  // converged
  //cout << "Iterations: " << iter << endl;
  return OK;
}

// ----------------------------------------------------------------------------
inline bool is_point_on_boundary(const SplineSurface& surf, const Point& uv)
// ----------------------------------------------------------------------------
{
  return (uv[0] <= surf.startparam_u() || uv[0] >= surf.endparam_u()   ||
	  uv[1] <= surf.startparam_v() || uv[1] >= surf.endparam_v());
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
  status = (n*n < SING_TOL)              ? SINGULAR :
           is_point_on_boundary(surf, p) ? BOUNDARY :
                                           REGULAR;
  return Array4
    { a * ds_dv,  // u-component of tangent
     -a * ds_du, // v-component of tangent
      da_dt * ds_dv + pow(a,2) * (d2s_dudv * ds_dv - d2s_dv2 * ds_du), // double deriv.
     -da_dt * ds_du - pow(a,2) * (d2s_du2 * ds_dv - d2s_dudv * ds_du) // ditto
    };
}

// ----------------------------------------------------------------------------
inline double signed_angle(const Point& dir1, const Point& dir2)
// ----------------------------------------------------------------------------
{
  double res = dir1.angle2(dir2);
  return (res < M_PI) ? res : 2 * M_PI - res;
}

// ----------------------------------------------------------------------------
inline bool in_between(const SplineSurface& surf, const PandDer& testpad,
		       const PandDer& startpad, const PandDer& endpad)
// ----------------------------------------------------------------------------
{
  const Point& testp  = testpad.first;
  const Point& startp = startpad.first;
  const Point& endp   = endpad.first;
  
  const double start_end_dist2 = startp.dist2(endp);
  const double dist2_1 = testp.dist2(startp);
  const double dist2_2 = testp.dist2(endp);

  // First, check whether the testpoint is "in-between" the start and end point,
  // distance-wise.
  if (dist2_1 > start_end_dist2 || dist2_2 > start_end_dist2)
    return false;

  // compare directions of tangents to intersection curve
  const Point start_dir {startpad.second[0], startpad.second[1]};
  const Point end_dir   {endpad.second[0],   endpad.second[1]};
  const Point test_dir  {testpad.second[0],  testpad.second[1]};

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
int closed_cycle_check(const SplineSurface& surf, const vector<PandDer>& points,
		       const size_t old_size)
// ----------------------------------------------------------------------------
{
  // check if the first point of the vector (start point) is found in-between
  // the newly inserted points, in which case we assume the cycle to have been
  // closed.
  const PandDer& startpoint = points.front();

  // if (points.size() > 20000) {
  //   cout << "exiting using hack " << endl;
  //   return 20000;
  // }
  
  // Initial check to see if we have moved far enough along the curve to be in
  // the vicinity of the curve's start point.
  if (startpoint.first.dist2(points.back().first) >
      points[old_size-1].first.dist2(points.back().first))
    return 0;
  
  for (size_t it = old_size; it != points.size(); ++it) 
    if (in_between(surf, startpoint, points[it-1], points[it])) 
      return int(it-1);

  // no end-of-cycle found
  return 0; 
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
Point mirror_tangent(const Point& p1, const Point& p2, const Point& tan2)
// ----------------------------------------------------------------------------
{
  // it is assumed that tan2 is of unit length
  Point d = p2-p1; d.normalize();
  const double sprod = tan2 * d;
  const double beta = sqrt(1-sprod*sprod);
  return Point {sprod * d[0] - beta * d[1],
                sprod * d[1] + beta * d[0]};
}

// ----------------------------------------------------------------------------
pair<Point, Point> identify_tangents(const PandDer& pstart, const PandDer& pend,
				     const bool forward)
// ----------------------------------------------------------------------------
{
  const int dir = forward ? 1 : -1;
  Point der_start {pstart.second[0] * dir, pstart.second[1] * dir};
  Point der_end   {pend.second[0]   * dir, pend.second[1]   * dir};

  // If a tangent is missing (singular point), mirror it with the other
  if (der_start.length2() == 0) 
    der_start = mirror_tangent(pstart.first, pend.first, der_end);
  else if (der_end.length2() == 0) 
    der_end   = mirror_tangent(pend.first, pstart.first, der_start);

  return {der_start, der_end};
}

// ----------------------------------------------------------------------------
bool check_midpoint(const SplineSurface& surf, const PandDer& pstart, const PandDer& pend,
		    const bool forward, const double isoval, const double tol,
		    PandDer& midpoint, PointIterationOutcome& outcome)
// ----------------------------------------------------------------------------
{
  // Function returns 'true' if midpoint is already OK
  outcome = OK; // this is the default case
  
  // Identifying tangents, and constructing one in case of degeneracy
  const pair<Point, Point> tangents = identify_tangents(pstart, pend, forward);
  const Point& der_start  = tangents.first;
  const Point& der_end    = tangents. second;
  
  const double arclen = estimate_arclength(pstart.first, pend.first, der_start, der_end);
  
  // Estimating position of new point (assuming a cubic hermite spline with arc
  // length parameterization, tangent scaled by approximate arc length)
  midpoint.first = 0.5 * (pstart.first + pend.first) + 0.125 * (der_start - der_end) * arclen;
  truncate_point_to_domain(surf, midpoint.first);
  
  // Check whether the new point is close enough to the isoval
  vector<Point> tmp(3); surf.point(tmp, midpoint.first[0], midpoint.first[1], 1);
  const double err = tmp[0][0] - isoval;
  const double grad2 = (tmp[1][0] * tmp[1][0]) + (tmp[2][0] * tmp[2][0]);
  if (pow(err, 2) <= (tol * tol * grad2))
    return true;

  // the midpoint is not already accurate enough.  Move it to isocontour and
  // compute associated derivative information
  
  outcome = move_point_to_isocontour(surf, isoval, tol, midpoint.first);
  PointStatus status;  // @@ Should we keep/return this value?
  midpoint.second = get_derivs(surf, midpoint.first, status); 
  return false;
}


// ----------------------------------------------------------------------------
vector<PandDer> intermediate_points(const SplineSurface& surf, const PandDer& old_point,
				    const PandDer& new_point, const bool forward,
				    const double isoval, const double tol, bool& ok)
// ----------------------------------------------------------------------------
{
  // Find enough intersection points between 'new_point', and 'old_point', to
  // ensure that the hermite interpolation will remain sufficiently close to the
  // real intersection curve at any time.

  // First, check if a point needs to be inserted at all.  If true, we will
  // establish a list of new points to insert

  ok = true;
  
  PandDer midpoint;
  PointIterationOutcome outcome;
  if (!check_midpoint(surf, old_point, new_point, forward, isoval, tol, midpoint, outcome) &
      (outcome == OK)) {
    
    // this curve needs to be refined.  But more than one point might be required.
    list<PandDer> tmp {old_point, midpoint, new_point};
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
    return vector<PandDer>(++tmp.begin(), --tmp.end()); // Skip first and last point
  } 
  
  // if we got here, no midpoint was needed, and we can return an empty vector
  ok = (outcome == OK);
  return vector<PandDer>();

}


// ----------------------------------------------------------------------------
double patchsize(const SplineSurface& surf, const Point& uv)
// ----------------------------------------------------------------------------
{
  const double TOL = 1e-8;// @@ knot tolerance.  Reasonable value?

  const double umin = surf.nextSegmentVal(0, uv[0], false, TOL);
  const double umax = surf.nextSegmentVal(0, uv[0], true, TOL);
  const double vmin = surf.nextSegmentVal(1, uv[1], false, TOL);
  const double vmax = surf.nextSegmentVal(1, uv[1], true, TOL);
  return min((umax-umin), (vmax-vmin));
}

// ----------------------------------------------------------------------------
inline double boundary_distance(const SplineSurface& surf, const PandDer& p,
				const bool v_dir, const bool forward)
// ----------------------------------------------------------------------------
{
  const double HUGEVAL = 1e100; // @@ something more elegant here?
  const int ix = v_dir ? 1 : 0;
  const double dir_component = p.second[ix] * (forward ? 1 : -1);
  const double boundary_comp_val =
    v_dir ? ((dir_component < 0) ? surf.startparam_v() : surf.endparam_v()) :
            ((dir_component < 0) ? surf.startparam_u() : surf.endparam_u());

  const double min_dist = fabs(p.first[ix] - boundary_comp_val);
  const double travel_dist = (dir_component == 0) ? HUGEVAL : min_dist/fabs(dir_component);

  return travel_dist;
}

// ----------------------------------------------------------------------------
double choose_steplength(const SplineSurface& surf, const PandDer& p,
			 const bool forward, const double tol)
// ----------------------------------------------------------------------------
{
  const double psize = patchsize(surf, p.first);
  const double K = 3;         // @@ this factor might be adapted
  const double maxlen = 0.2 * psize; 
  const Array4& derivs = p.second;
  const double bd_u = boundary_distance(surf, p, false, forward);
  const double bd_v = boundary_distance(surf, p, true,  forward);
    
  const double d1 = max(fabs(derivs[0]), fabs(derivs[1]));
  const double d2 = max(fabs(derivs[2]), fabs(derivs[3]));

  const double dt = min({bd_u+tol, bd_v+tol, maxlen, 2 * d1 / (K * d2)});

  return dt;
}

// ----------------------------------------------------------------------------
PandDer extrapolate_point(const SplineSurface& surf, double dt, const PandDer& cur_pander,
			  const double isoval, double tol, PointStatus& status)
// ----------------------------------------------------------------------------
{
  // determine a new test point based on previous point and derivatives.
  // Converge the new point so that it lies on the isocurve.  Check status of
  // the new point (could be REGULAR, BOUNDARY or SINGULAR).
  const Point& cur_point = cur_pander.first;
  const Array4& derivs   = cur_pander.second;
  
  // Loop until a suitable point is found
  bool found = false;
  Point new_point;
  while (!found) {
    // First, use derivatives to estimate position of new point
    new_point = cur_point;
    for (int i = 0; i != 2; ++i)
      new_point[i] += dt * derivs[i] + 0.5 * dt * dt * derivs[i+2];
    truncate_point_to_domain(surf, new_point);
    
    // move estimated point to an actual isocontour position
    switch (move_point_to_isocontour(surf, isoval, tol, new_point)) {
    case OK:
      // check if tangent is consistent with previous one (otherwise we have
      // likely jumped to a neighbor curve
      {
	const auto derivs_new = get_derivs(surf, new_point, status);
	const double sprod = derivs_new[0] * derivs[0] + derivs_new[1] * derivs[1];
	if (sprod < 0) {
	  // the new point is most likely not on the correct curve
	  dt = dt/2;
	  break;
	} else {
	  status =  is_point_on_boundary(surf, new_point) ? BOUNDARY : REGULAR;
	  found = true;
	}
      }
      break;
    case ITERATION_EXCEED:
      // Could not converge to new point.  Shorten steplength and try again
      dt = dt/2;
      break;
    case FLAT_REGION_ENCOUNTERED:
      // check if singularity is on the curve.  If not, shorten timestep and try again
      if (fabs(dt) < tol) {  
	found = true;
	status = is_point_on_boundary(surf, new_point) ? BOUNDARY : SINGULAR;
      } else {
	dt = dt/2;
      }
      break;
    }
  }
  // If we have not been able to converge to a new point, flag it as singular.
  // @@ Validity of this should be tested
  // if (cur_point.dist2(new_point) < tol * tol)
  //   status = SINGULAR;

  PointStatus tmp_status; // just a throwaway value - we already have a more
			  // precise assessment of the status of this point
  return PandDer
    { new_point,
      get_derivs(surf, new_point, tmp_status) // in case the last point is not
	                                      // our original extrapolated
					      // point, we have to do the
					      // gradient calculation again
    };
}

// ----------------------------------------------------------------------------
PointStatus find_next_point(const SplineSurface& surf, vector<PandDer>& prev_points,
			    const bool forward, const double isoval,
			    const double tol, const double fac)
// ----------------------------------------------------------------------------
{
  PointStatus status = REGULAR; // Will be changed in functions call below
  // upon entry, 'derivs' should contain the derivatives corresponding to the
  // last point in the 'prev_points' vector.  At exit, 'derivs' will contain the
  // derivatives corresponding to the last point added.  The use of this
  // variable is for efficiency only.
  const double dt =
    choose_steplength(surf, prev_points.back(), forward, tol) * fac;

  // returned status here can be REGULAR, BOUNDARY or SINGULAR (CYCLIC_END will
  // be checked for later).  Since the current point is REGULAR, we should be
  // guaranteed that another point can be produced.
  const PandDer new_pt = extrapolate_point(surf, (forward ? dt : -dt),
					   prev_points.back(), isoval, tol, status);

  // ensure that the curve between previous and new point will stick closely
  // enough to the actual intersection
  bool ok; // if a topology problem is detected (points lying on different
	   // curves), ok will be set to 'false' upon return of the function
	   // call below. In that case, the returned vector will contain the
	   // points assumed to be lying on the present curve (but the new point
	   // will not, and should thus be discarded).
  const size_t prev_size = prev_points.size();
  const vector<PandDer> ipoints =
    intermediate_points(surf, prev_points.back(), new_pt, forward, isoval, tol, ok);

  prev_points.insert(prev_points.end(), ipoints.begin(), ipoints.end());

  // Add the new point if no topology issue arose, and ascertain that at least
  // one point has been added.
  if (ok) 
    prev_points.emplace_back(new_pt);
  else if (ipoints.empty() && (fabs(dt) > 2 * tol)) 
    // no next point was ultimately added.  Unless step size limit is reached,
    // call routine again, with smaller step size.
    return find_next_point(surf, prev_points, forward, isoval, tol, fac/2);
  
  // check for possibly closed cycle
  const int ix = closed_cycle_check(surf, prev_points, prev_size);
  if (ix > 0) {
    status = CYCLIC_END;
    assert(ix >= int(prev_size)-1);
    prev_points.erase(prev_points.begin() + ix, prev_points.end());
  }

  // Ensuring point status is correct
  // @@ could this information be preserved from earlier computations?
  PointStatus ps_tmp;
  get_derivs(surf, prev_points.back().first, ps_tmp);
  if (status != CYCLIC_END)
    status = ps_tmp;
  
  return status;
}

// ----------------------------------------------------------------------------
inline bool moving_inwards(const SplineSurface& surf,
			   const bool forwards,
			   const PandDer& pt)
// ----------------------------------------------------------------------------
{
  const int sign = forwards ? 1 : -1;
  if ((pt.first[0] <= surf.startparam_u()) && (pt.second[0] * sign < 0))
    return false; // exiting u_min boundary
  if ((pt.first[0] >= surf.endparam_u()) && (pt.second[0] * sign > 0))
    return false; // exiting u_max boundary
  if ((pt.first[1] <= surf.startparam_v()) && (pt.second[1] * sign < 0))
    return false; // exiting v_min boundary
  if ((pt.first[1] >= surf.endparam_v()) && (pt.second[1] * sign > 0))
    return false; // exiting v_max boundary

  return true;
}

// ----------------------------------------------------------------------------
vector<PandDer> trace_unidir(const SplineSurface& surf, const Point& startpoint,
			   const bool forward, double tol, PointStatus& last_point_status)
// ----------------------------------------------------------------------------
{
  const double isoval = value_at(surf, startpoint); // the isovalue at which the
						    // curve should lie
  
  // get curve tangent, second derivatives, and check the status of the
  // startpoint (REGULAR, BOUNDARY, etc.)
  const Array4 derivs = get_derivs(surf, startpoint, last_point_status); 

  // Establishing the result vector
  vector<PandDer> result(1, PandDer {startpoint, derivs});
  
  // we are at a regular point inside the domain, keep tracing
  while ((last_point_status == REGULAR) |
	 ((last_point_status == BOUNDARY) && (moving_inwards(surf, forward, result.back()))))
    last_point_status = find_next_point(surf, result, forward, isoval, tol, 1);

  // only include the startpoint itself if we have been tracing forward
  if (!forward)
    result.erase(result.begin()); 
				  
  return result;
}

// ----------------------------------------------------------------------------
pair<CurvePtr, CurvePtr> curve_from_points(const SplineSurface& surf,
					   const vector<PandDer>& pvec,
					   bool include_3D)
// ----------------------------------------------------------------------------  
{
  if (pvec.size() < 2) {
    cout << "Warning: degenerate curve (single point).  Returning empty curve." << endl;
    return pair<CurvePtr, CurvePtr>();
  }
    
  // make a cubic hermite spline based on these points

  // compute the arc length for each interval
  const int num_intervals = (int)pvec.size() - 1;
  vector<double> alengths(num_intervals);
  transform(pvec.begin(), pvec.end()-1, pvec.begin()+1, alengths.begin(),
	    [] (const PandDer& p1, const PandDer& p2) {
	      return estimate_arclength(p1.first, p2.first,
					Point {p1.second[0], p1.second[1]},
					Point {p2.second[0], p2.second[1]});});

  // compute individual knot values
  vector<double> kvals(pvec.size(), 0); // first knotval will be zero.  Compute the rest
  partial_sum(alengths.begin(), alengths.end(), kvals.begin()+1);

  // Constructing knot vector (triple interior knots, quadruple knots at endpoints)
  vector<double> kvec(kvals.size() * 3 + 2, 0);
  for (auto vals = kvals.begin(), target = kvec.begin() + 1; vals != kvals.end(); ++vals) {
    for (int i = 0; i != 3; ++i)
      *target++ = *vals;
  }
  kvec.back() = kvec[kvec.size() - 2]; // quadruple last know
  
  // computing control point values
  const int dim = 2;
  vector<double> coefs(dim * (3 * num_intervals + 1));

  // filling in the control points lying on the curve
  for (size_t i = 0; i != pvec.size(); ++i) {
    const size_t ix = dim * 3 * i;
    coefs[ix]   = pvec[i].first[0];
    coefs[ix+1] = pvec[i].first[1];
  }

  // filling in the control points that specify tangent directions

  // left tangent (nb: note the end criteron in the for loop)
  for (size_t i = 0; i != pvec.size()-1; ++i) {
    const size_t left_ix = dim * (3 * i + 1);
    const double len = alengths[i]/3;
    coefs[left_ix    ] = coefs[dim * (3 * i)    ] + len * pvec[i].second[0];
    coefs[left_ix + 1] = coefs[dim * (3 * i) + 1] + len * pvec[i].second[1];
  }

  // right tangent (nb: note different range in the for loop)
  for (size_t i = 1; i != pvec.size(); ++i) {
    const size_t right_ix = dim * (3 * i - 1);
    const double len = alengths[i-1]/3;
    coefs[right_ix    ] = coefs[dim * (3 * i)    ] - len * pvec[i].second[0];
    coefs[right_ix + 1] = coefs[dim * (3 * i) + 1] - len * pvec[i].second[1];
  }

  // Constructing the final spline curve in parameter space
  const int num_pts = (int)coefs.size() / dim;    
  const CurvePtr paramcurve {new SplineCurve(num_pts, 4, kvec.begin(), coefs.begin(), dim)};
  
  // constructing 3D curve if requested
  CurvePtr spacecurve;
  if (include_3D) {
    const double isoval = value_at(surf, pvec[0].first);
    vector<double> coefs3D(num_pts * 3);
    for (int i = 0; i != num_pts; ++i) {
      coefs3D[3 * i    ] = coefs[2 * i];
      coefs3D[3 * i + 1] = coefs[2 * i + 1];
      coefs3D[3 * i + 2] = isoval;
    }
    spacecurve = CurvePtr{new SplineCurve(num_pts, 4, kvec.begin(), coefs3D.begin(), 3)};
  }
  
  return pair<CurvePtr, CurvePtr> {paramcurve, spacecurve};
}


// ----------------------------------------------------------------------------  
pair<CurvePtr, CurvePtr> trace_isoval(double u, double v, const SplineSurface& surf,
				      const double tol, const bool include_3D)
// ----------------------------------------------------------------------------
{
  PointStatus ps;
    
  // trace in the first direction
  const auto res1 = trace_unidir(surf, {u, v}, true, tol, ps);

  // if we did not get back to where we began, trace in second direction also,
  // otherwise, repeat initial point to create a closed loop
  const auto res =
    (ps == CYCLIC_END) ?
    insert_back(res1, res1.front()) :
    merge_vec(reverse_vec(trace_unidir(surf, {u, v}, false, tol, ps)), res1);

  return curve_from_points(surf, res, include_3D);
}

// ----------------------------------------------------------------------------
pair<CurvePtr, CurvePtr> trace_isoval_sisl(SISLIntcurve* ic, SISLSurf* s,
					   const double isoval,
					   const double tol, const bool include_3D)
// ----------------------------------------------------------------------------
{
  double pnt[] = {0, 0, isoval};
  double nrm[] = {0, 0, 1};
  const double epsge = tol; //1e-6;
  const double epsco = 1e-15; // Not used
  const double maxstep = 0.0;
  const int dim = 3;
  const int makecurv = 2; // make both geometric and parametric curves
  int stat;
  
  s1314(s, pnt, nrm, dim, epsco, epsge, maxstep, ic, makecurv, 0, &stat);

  SISLCurve* sc = ic->pgeom;
  SISLCurve* sp = ic->ppar1;
  if (sc == 0) {
    MESSAGE("s1314 returned code: " << stat << ", returning.");
    return pair<CurvePtr, CurvePtr>();
  }
  
  assert(sc->rcoef == 0); // otherwise, slight modification of the below is necessary
  const bool rational = false;

  // 3D curve will be produced whether requested or not, since s1314 cannot generate parameter
  // curve without also generating 3D curve.
  return pair<CurvePtr, CurvePtr>
  {
    CurvePtr(new SplineCurve(sp->in, sp->ik, sp->et, sp->ecoef, 2, rational)), // pcurve
    CurvePtr(new SplineCurve(sc->in, sc->ik, sc->et, sc->ecoef, 3, rational)) // scurve
  };
}

// ----------------------------------------------------------------------------
CurveVec compute_levelset(const SplineSurface& ss,
			  SISLSurf* ss_sisl, // SISL surface corresponding to ss
			  SISLSurf* ss_sisl_3D, // only nonzero if use_sisl_marching = true
			  const double isoval,
			  const double tol, 
			  const bool include_3D_curves,
			  const bool use_sisl_marching)
// ----------------------------------------------------------------------------
{
  //cout << "Now computing level set for " << isoval << endl;
  //compute topology
  const pair<SISLIntcurve**, int> topo_pts = compute_topology(ss_sisl, isoval);

  // marching out curves, using either SISL routine s1314 or (quicker) native routine
  using MarchFun = function<pair<CurvePtr, CurvePtr>(SISLIntcurve*)>;
  const MarchFun sisl_mfun {[&] (SISLIntcurve* ic)
      {return trace_isoval_sisl(ic, ss_sisl_3D, isoval, tol, include_3D_curves);}};
  const MarchFun ntve_mfun {[&] (SISLIntcurve* ic)
      {return trace_isoval(ic->epar1[0], ic->epar1[1], ss, tol, include_3D_curves);}};
    
  const CurveVec result = use_sisl_marching ?
    apply_transform(topo_pts.first, topo_pts.first + topo_pts.second, sisl_mfun) : 
    apply_transform(topo_pts.first, topo_pts.first + topo_pts.second, ntve_mfun);
    
  // cleaning up
  if (topo_pts.second > 0)
    freeIntcrvlist(topo_pts.first, topo_pts.second);

  return result;
}

}; // end anonymous namespace 
