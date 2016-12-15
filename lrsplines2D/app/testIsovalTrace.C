#include <iostream>
#include <fstream>
#include <vector>

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"

#include "GoTools/lrsplines2D/LRSplineSurface.h"

using namespace Go;
using namespace std;

namespace {
  enum PointStatus {REGULAR = 0, BOUNDARY=1, CYCLIC_END = 3, SINGULAR = 4};
  using CurvePtr = shared_ptr<const SplineCurve>;

  CurvePtr traceIsoval(const SplineSurface& surf, double u, double v);
  
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
  const CurvePtr result = traceIsoval(surf, u, v);


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
array<double, 4> get_derivs(const SplineSurface& surf, const Point& p,
			    PointStatus& status)
// ----------------------------------------------------------------------------
{
  return array<double, 4> {0, 0, 0, 0}; // @@ dummy
}

// ----------------------------------------------------------------------------
Point find_next_point(const SplineSurface& surf, const vector<Point>& prev_points,
		      bool forward, PointStatus& status)
// ----------------------------------------------------------------------------
{
  return Point(); //@@ Dummy
}
// ----------------------------------------------------------------------------
vector<Point> trace_unidir(const SplineSurface& surf, const Point& startpoint,
			   const bool forward, PointStatus& last_point_status)
// ----------------------------------------------------------------------------
{
  vector<Point> result(1, startpoint);

  // check the status of the startpoint (REGULAR, BOUNDARY, etc.)
  get_derivs(surf, startpoint, last_point_status); 

  while (last_point_status == REGULAR) {
    // we are at a regular point inside the domain, keep tracing
    result.emplace_back(find_next_point(surf, result, forward, last_point_status));
  }

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
