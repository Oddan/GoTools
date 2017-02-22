#include <algorithm>
#include "tritools.h"

#include "ttl/ttl.h"
#include <ttl/halfedge/HeDart.h>
#include <ttl/halfedge/HeTraits.h>


using namespace Go;
using namespace std;
using namespace ttl;

namespace {

inline double cross_value (const Point& p1, const Point& p2)
{
  return p1[0]*p2[1] - p2[0] * p1[1];
}  

// this function assumes that the line equation is on normalized form, i.e. a^2 + b^2 = 1,
// where a = line_equation[0] and b = line_equation[1].
inline double signed_distance(const array<double, 3>& line_equation, const Point& p)
{
  return line_equation[0] * p[0] + line_equation[1] * p[1] + line_equation[2];
}

inline array<double, 4> get_tribox(const Point& p1, const Point& p2, const Point& p3)
{
  return
    { min(min(p1[0], p2[0]), p3[0]),  // xmin
      max(max(p1[0], p2[0]), p3[0]),  // xmax
      min(min(p1[1], p2[1]), p3[1]),  // ymin 
      max(max(p1[1], p2[1]), p3[1])   // ymax
    }; 
}
  
}; // end anonymous namespace 

namespace TriTools {

// ----------------------------------------------------------------------------
array<double, 3> normalized_line_equation(const Point& pt_first,
					  const Point& pt_second)
// ----------------------------------------------------------------------------
{
  const Point dp = pt_second - pt_first;
  const double inv_length = 1.0/dp.length();

  return {
     dp[1] * inv_length,
    -dp[0] * inv_length,
    -cross_value(pt_first, pt_second) * inv_length
  };
}

// ----------------------------------------------------------------------------
vector<int> points_inside_triangle(const Point& p1,
				   const Point& p2,
				   const Point& p3,
				   const vector<Point>& points,
				   const array<double, 3>& eps)
// ----------------------------------------------------------------------------
{
  vector<int> result(points.size(), 0);
  transform(points.begin(), points.end(), result.begin(), [&] (const Point& p) {
      return (int)(
	(signed_distance(normalized_line_equation(p1, p2), p) <= -eps[0]) &
	(signed_distance(normalized_line_equation(p2, p3), p) <= -eps[1]) &
	(signed_distance(normalized_line_equation(p3, p1), p) <= -eps[2]));});
  return result;
}

// ----------------------------------------------------------------------------
vector<int> points_inside_loops(const vector<vector<Point>>& loops,
				const vector<Point>& points,
				const double margin)
// ----------------------------------------------------------------------------
{
  const vector<pair<array<Point, 3>, array<bool, 3>> triangles =
	       triangulate_boundary(loops);

  vector<int> result(points.size(), 0); // we initially consider all points outside
  for (const auto tri : triangles) {
    auto cur_inside = points_inside_triangle(tri.first[0],
					     tri.first[1],
					     tri.first[2],
					     points,
					     {(tri.second[0] ? margin : 0),
					      (tri.second[1] ? margin : 0),
					      (tri.second[2] ? margin : 0)});
    transform(result.begin(), result.end(), cur_inside.begin(), result.begin(),
	      [](int a, int b) {return int(a || b);});
  }
  return result;
}

// ----------------------------------------------------------------------------
std::vector<std::pair<std::array<Go::Point, 3>, std::array<bool, 3>>>
triangulate_boundary(const std::vector<std::vector<Go::Point>>& loops)
// ----------------------------------------------------------------------------
{
  // Create delauney triangulation containing all boundary points
  Triangulation triang = create_delaunay_from_loops(loops);

  // constrain boundary edges
  for (l : loops)
    impose_boundary(triang, l.begin(), l.end());

  // generate triangles of the result
  return export_triangles(triang);
  
}

  
}; // end namespace TriTools
