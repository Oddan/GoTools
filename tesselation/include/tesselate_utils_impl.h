#ifndef _TESSELATE_UTILS_IMPL_H
#define _TESSELATE_UTILS_IMPL_H


#include <assert.h>
#include <cmath>
#include <random>
namespace TesselateUtils
{

// ----------------------------------------------------------------------------
template<typename T> inline
void operator += (std::vector<T>& lhs, const std::vector<T>& rhs)
// ----------------------------------------------------------------------------
{
  assert(lhs.size() == rhs.size());
  auto rhs_it = rhs.cbegin();
  for(auto lhs_it = lhs.begin(); lhs_it != lhs.end(); ++lhs_it, ++rhs_it)
    *lhs_it += *rhs_it;
}

// ----------------------------------------------------------------------------
template<typename P>
std::vector<P> interpolate(const P& p1, const P& p2, unsigned int num)
// ----------------------------------------------------------------------------  
{
  std::vector<double> param(num, double(1)/(num+1));
  std::partial_sum(param.begin(), param.end(), param.begin());

  std::vector<P> result(1, p1);
  std::transform(param.begin(), param.end(), back_inserter(result),
		 [&p1, &p2](double t) {return p1 * (1-t) + p2 * t;});
  
  result.push_back(p2);
  return result;
};

// ----------------------------------------------------------------------------
template<typename P>
double polygon_area(const P* const poly, const unsigned int num_corners)
// ----------------------------------------------------------------------------  
{
  double area = 0;

  for (int i = 0; i != (int)num_corners; ++i) {
    const P& p1 = poly[i];
    const P& p2 = poly[(i+1)%num_corners];

    area += 0.5 * (-p1[1] * p2[0] + p1[0] * p2[1]);
  }

  return area;
};

// ----------------------------------------------------------------------------    
template<typename P>  
std::array<double, 4> bounding_box(const P* const points,
				   const unsigned int num_points)
// ----------------------------------------------------------------------------    
{
  assert(num_points > 0);

  const auto minmax_x = std::minmax_element(points, points + num_points,
					    [](P p1, P p2) {return p1[0] < p2[0];});
  const auto minmax_y = std::minmax_element(points, points + num_points,
					    [](P p1, P p2) {return p1[1] < p2[1];});

  return std::array<double, 4> { (*(minmax_x.first))[0], (*(minmax_x.second))[0],
                                 (*(minmax_y.first))[1], (*(minmax_y.second))[1]};
}

// ----------------------------------------------------------------------------
template<typename P>
std::vector<P> generate_grid_2D(const P& c1, const P& c2, int nx, int ny)
// ----------------------------------------------------------------------------  
{
  // Generate vectors of points demarcating the leftmost and rightmost grid columns.
  const auto yminvec = interpolate(c1, {c1[0], c2[1]}, ny);
  const auto ymaxvec = interpolate({c2[0], c1[1]}, c2, ny);

  // Generate grid rows
  std::vector<std::vector<P>> gridrows;
  std::transform(yminvec.begin(), yminvec.end(), ymaxvec.begin(),
		 std::back_inserter(gridrows),
		 [nx](P p1, P p2) {return interpolate(p1, p2, nx);});

  return flatten(gridrows);
}

// ----------------------------------------------------------------------------  
template<typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& arg)
// ----------------------------------------------------------------------------    
{
  std::vector<T> result;
  for (auto s = arg.begin(); s != arg.end(); ++s)
    result.insert(result.end(), s->begin(), s->end());
  return result;
}

// ----------------------------------------------------------------------------
// Check if a point is inside a polygon
template<typename P>
bool inpolygon(const P& pt, const P* const poly, 
	       const unsigned int num_corners, const double tol)
// ----------------------------------------------------------------------------      
{
  // count the number of intersections of polygonal line segments and the ray
  // from pt "eastward" (i.e. [pt, {inf, pt[1]}]).  An odd number of
  // intersections signifies that the point is inside the polygon.  

  int isects = 0;
  for (unsigned int i = 0; i != num_corners; ++i) {
    const P& p1 = poly[i];
    const P& p2 = poly[(i+1) % num_corners];

    // First, check if point is _on_ the edge, in which case we exclude it (we
    // only seek those that are true interior points)
    if (point_on_line_segment(pt, p1, p2, tol))
      return false;
    
    // case where startpoint of segment is exactly on the ray
    if (p1[1] == pt[1]) {
      isects += int(p1[0] > pt[0]); // must be to the right to be on the ray
      continue;
    } else {
      // we know that the startpoint is not _exactly_ on the ray
      if ((p1[1] - pt[1]) * (p2[1] - pt[1]) > 0)
	// Both segment points are on the same vertical side of ray, no
	// intersection possible
	continue;
      else if ((p1[0] <= pt[0]) && (p2[0] <= pt[0])) 
	// both segments are to the left of ray origin - no intersection
	// possible
	continue;
      else {
	//there is one point below and one above the ray, and at least one of
	//them are at the right of the ray origin
	if ((p1[0] > pt[0]) && (p2[0] > pt[0])) {
	  // both points are to the right of ray origin.  Intersection is
	  // ensured
	  isects += 1;
	  continue;
	} else {
	  // One point is to the left and one to the right of ray origin.
	  // Moreover, one point is above and one below the ray origin.
	  
	  // compute x-coordinate where segment crosses y=pt[1] axis
	  const double t = (pt[1] - p1[1]) / (p2[1] - p1[1]);
	  const double x = t * p2[0] + (1-t) * p1[0];
	  isects += int(x >= pt[0]);
	}
      }
    }
  }
  
  return (isects % 2 == 1);
  
}
  
// ----------------------------------------------------------------------------
// Return the points that are inside the polygon
template<typename P>
std::vector<P> inpolygon(const P* const pts, const unsigned int num_pts,
			 const P* const poly, const unsigned int num_corners,
			 const double tol)
// ----------------------------------------------------------------------------    
{
  std::vector<P> result;
  std::copy_if(pts, pts + num_pts, std::back_inserter(result),
	       [&](const P& p) {return inpolygon(p, poly, num_corners, tol);});
  return result;
}

// ----------------------------------------------------------------------------    
template<typename P> inline
double projected_distance_to_line(P p, P a, P b)
// ----------------------------------------------------------------------------    
{
  const double dab = dist(a,b);
  p -= a; // now represents ap
  b -= a; // now represents ab

  return (p[0]*b[1] - p[1]*b[0]) / dab;

  // const P ap = p-a;
  // const P ab = b-a;
  
  //  return std::abs((ap[0]*ab[1] - ap[1]*ab[0]) / dist(a,b));
}

// ----------------------------------------------------------------------------    
// Check if the point p is within distance 'tol' of the line segment defined by
// points a and b
template<typename P> inline
bool point_on_line_segment(const P& p, const P& a, const P& b, double tol)
// ----------------------------------------------------------------------------    
{
  const double tol2 = tol * tol;
  const double dpa2 = dist2(p,a);
  if (dpa2 < tol2) return true;

  const double dpb2 = dist2(p,b);
  if (dpb2 < tol2) return true;

  // we need to establish a maximum distance from endpoints such that no point
  // outside this radius have a chance of being close enough to the line
  // segment.  We ensure this by multiplying the distance by sqrt(5/4), which
  // becomes 5/4 since we work with squared distances.  The factor 5/4 is easily
  // derived using Pythagoras and considering the distance to the line segment
  // midpoint
  const double dab2 = dist2(a,b);
  const double max_rad2 = std::max(dab2, tol2) * 1.25; // 1.25 = 5/4.
  if ((dpa2 > max_rad2) || (dpb2 > max_rad2)) return false;

  // we now know the point is _not_ coincident with the segment endpoints, but
  // still close enough that it remains a candidate.

  // does the point project outside the segment, if so, it is not on the segment
  if ((dab2 < dpa2) || (dab2 < dpb2))
    return false;

  // The only thing that remains to be checked now is whether the orthogonally
  // projected distances is within the tolerance
  const double d = std::abs(projected_distance_to_line(p, a, b));

  return (d < tol);
  
  // if (d > tol) return false;

  // // the point is close enough to the infinite line, but is it close enough to
  // // the _segment_?
  // return (dab2 > dpa2) && (dab2 > dpb2);
  // // const double ap = dist2(a,p);  // @ we can use dist2 instead of dist, since 
  // // const double ab = dist2(a,b);  //   the square function is convex (here, it
  // // const double bp = dist2(b,p);  //   prevents us from doing the costly square
  // // return (ab > ap) && (ab > bp); //   root operation)
}

// ----------------------------------------------------------------------------
template<typename P> inline
bool acute_angle(const P& a, const P& b, const P& c)
// ----------------------------------------------------------------------------
{
  // angle is acute if scalar product of vectors ba and bc is positive
  return (c[0] - b[0]) * (a[0] - b[0]) + (c[1] - b[1]) * (a[1] - b[1]) > 0;
}

// ----------------------------------------------------------------------------
template<typename P> inline
bool projects_to_segment(const P& p, const P& a, const P&b)
// ----------------------------------------------------------------------------    
{
  return acute_angle(p, a, b) && acute_angle(p, b, a);
}

// ----------------------------------------------------------------------------    
inline double random_uniform(double minval, double maxval)
// ----------------------------------------------------------------------------    
{
  assert(maxval > minval);
  static std::uniform_real_distribution<double> unif(minval, maxval);
  static std::default_random_engine re;
  
  return unif(re);
}

// ----------------------------------------------------------------------------    
template<typename T> inline
std::vector<unsigned int> locate_nonzeros(const std::vector<T> vec)
// ----------------------------------------------------------------------------    
{
  std::vector<unsigned int> result;
  locate_nonzeros(vec, result);
  return result;
}

// ----------------------------------------------------------------------------    
template<typename T> inline
void locate_nonzeros(const std::vector<T> vec, std::vector<unsigned int>& target)
// ----------------------------------------------------------------------------
{
  target.clear();
  const unsigned int N = (unsigned int) vec.size();
  for (unsigned int i = 0; i != N; ++i) 
    if (vec[i] != 0)
      target.push_back(i);
}

// ----------------------------------------------------------------------------
template<typename T, typename Pred> inline
std::pair<std::vector<T>, std::vector<unsigned int>>
extract_from_range(const T* const start, unsigned int len, const Pred& fun)
// ----------------------------------------------------------------------------  
{
  std::pair<std::vector<T>, std::vector<unsigned int>> result;

  std::vector<int> flag(len);
  transform(start, start+len, flag.begin(), fun);
  result.second = locate_nonzeros(flag);
  result.first.reserve(result.second.size());
  transform(result.second.begin(), result.second.end(),
	    back_inserter(result.first),
	    [start](unsigned int ix) {return start[ix];});
  return result;
}

// ----------------------------------------------------------------------------
template<typename P> inline
P mirror_point(P p, const P& a, const P& b)
// ----------------------------------------------------------------------------
{
  // compute normalized tangent vector and normal vector
  const P t = (a - b) / dist(a, b); // @@ unnecessary intermediary.  May be optimized
  const P n {-t[1], t[0]};

  // translate so that segment line passes through origo at 'a'
  p -= a;

  // compute mirroring matrix (from Householder).  Column-based storage.
  const double nxny = n[0] * n[1];
  const std::array<double, 4> H {1 - 2*n[0]*n[0], -2*nxny, -2*nxny, 1-2*n[1]*n[1]};

  // compute mirrored point by multplying with matrix and translate back by adding a
  return P { H[0] * p[0] + H[2] * p[1] + a[0],   
             H[1] * p[0] + H[3] * p[1] + a[1]};
}

// ----------------------------------------------------------------------------        
template<typename P> inline
std::vector<P> mirror_points(const P* const pts,
			     const unsigned int num_pts,
			     const P& a,
			     const P& b)
// ----------------------------------------------------------------------------
{
  std::vector<P> result(num_pts);
  transform(pts, pts + num_pts, result.begin(),
	   [&a, &b] (const P& p) {return mirror_point(p, a, b);});
  return result;
}

// ----------------------------------------------------------------------------
template<typename P> inline
P solve_2D_matrix(const P& Mcol1, const P& Mcol2, const P& rhs)
// ----------------------------------------------------------------------------
{
  const double det = Mcol1[0] * Mcol2[1] - Mcol1[1] * Mcol2[0];
  if (det==0)
    return P {std::nan(""), std::nan("")};

  const double Dx = rhs[0] * Mcol2[1] - rhs[1] * Mcol2[0];
  const double Dy = Mcol1[0] * rhs[1] - Mcol1[1] * rhs[0];

  return P {Dx/det, Dy/det};
}

// ----------------------------------------------------------------------------
template<typename P> inline
bool circumscribe_triangle(const P& p1, const P& p2, const P& p3,
			   P& circ_center, double& radius2)
// ----------------------------------------------------------------------------
{
  const double p1_2 = norm2(p1);
  const double p2_2 = norm2(p2);
  const double p3_2 = norm2(p3);
  circ_center = solve_2D_matrix(P {p1[0] - p2[0], p1[0] - p3[0]},
				P {p1[1] - p2[1], p1[1] - p3[1]},
				P {0.5 * (p1_2 - p2_2), 0.5 * (p1_2 - p3_2)});
  radius2 = dist2(p1, circ_center);
  return !(std::isnan(circ_center[0]));
}

// ----------------------------------------------------------------------------
// check whether two segments intersect.  A positive value for the tolerance
// means that a vincinity within the tolerance counts as an intersection.  A
// negative value for the tolerance, on the other hand, means that each line
// segment must cross the other by a length of at least |tol| times its length.
// In this case, a mere touch is not sufficient to count as a intersection.
template<typename P> inline
bool segments_intersect(const P& seg1_a, const P& seg1_b,
			const P& seg2_a, const P& seg2_b, const double tol)
// ----------------------------------------------------------------------------
{
  // first, compare bounding boxes, to eliminate obvious cases
  if (std::min(seg1_a[0], seg1_b[0]) > std::max(seg2_a[0], seg2_b[0]) + tol)
    return false;
  if (std::min(seg1_a[1], seg1_b[1]) > std::max(seg2_a[1], seg2_b[1]) + tol)
    return false;

  // OK, bounding boxes overlap.  Check for intersection.
  // We determine u and v such that: 
  // (1-u) * seg1_a + u * seg1_b = (1-v) * seg2_a + v * seg2_b

  const P uv = solve_2D_matrix( P {seg2_b[0] - seg2_a[0], seg2_b[1] - seg2_a[1]},
                                P {seg1_a[0] - seg1_b[0], seg1_a[1] - seg1_b[1]},
				P {seg1_a[0] - seg2_a[0], seg1_a[1] - seg2_a[1]});
  if (std::isnan(uv[0])) {
    // degenerate case - segments are collinear.  If tolerance is > 0, we must
    // check whether they are indeed on the same line.  (If tolerance is < 0, we
    // should not consider this and intersection).
    if (tol < 0)
      return false;
      
    const double l = std::min(dist(seg1_a, seg1_b), dist(seg2_a, seg2_b));
    return ( std::abs(projected_distance_to_line(seg1_a, seg2_a, seg2_b)) < (tol * l) );
  }

  // OK.  Our system is not degenerate.  Let us check whether an actual
  // intersection occurs.
  return ((uv[0] > -tol) && (uv[0] < 1 + tol) && (uv[1] > -tol) && (uv[1] < 1 + tol));
}
  

}; // end namespace


#endif
