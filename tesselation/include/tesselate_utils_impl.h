#ifndef _TESSELATE_UTILS_IMPL_H
#define _TESSELATE_UTILS_IMPL_H


#include <assert.h>
#include <random>
namespace TesselateUtils
{

// ----------------------------------------------------------------------------
template<typename T> inline
void operator += (std::vector<T>& lhs, const std::vector<T>& rhs)
// ----------------------------------------------------------------------------
{
  auto rhs_it = rhs.cbegin();
  for(auto lhs_it = lhs.begin(); lhs_it != lhs.end(); ++lhs_it, ++rhs_it)
    *lhs_it += *rhs_it;
}

// ----------------------------------------------------------------------------
template<typename P> inline
void operator += (P& p1, const P& p2)
// ----------------------------------------------------------------------------
{ 
  p1[0] += p2[0]; 
  p1[1] += p2[1];
}

// ----------------------------------------------------------------------------
template<typename P> inline
P operator - (const P& p1, const P& p2)
// ----------------------------------------------------------------------------
{
  return P {p1[0] - p2[0], p1[1] - p2[1]};
}

// ----------------------------------------------------------------------------
template<typename P> inline 
P operator + (const P& p1, const P& p2)
// ----------------------------------------------------------------------------
{
  return P {p1[0] + p2[0], p1[1] + p2[1]};
}

// ----------------------------------------------------------------------------
template<typename P> inline
P operator*(const P& p, double t)
// ----------------------------------------------------------------------------
{
  return P {t*p[0], t*p[1]};
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
template<typename P>
double projected_distance_to_line(const P& p, const P& a, const P& b)
// ----------------------------------------------------------------------------    
{
  const P ap = p-a;
  const P ab = b-a;
  
  return std::abs((ap[0]*ab[1] - ap[1]*ab[0]) / dist(a,b));
}

// ----------------------------------------------------------------------------    
// Check if the point p is within distance 'tol' of the line segment defined by
// points a and b
template<typename P>
bool point_on_line_segment(const P& p, const P& a, const P& b, double tol)
// ----------------------------------------------------------------------------    
{
  if (dist(p, a) < tol) return true;
  if (dist(p, b) < tol) return true;

  // we now know the point is _not_ coincident with the segment endpoints
  const double d = projected_distance_to_line(p, a, b);
  if (d > tol) return false;

  // the point is close enough to the infinite line, but is it close enough to
  // the _segment_?
  const double ap = dist(a,p);
  const double ab = dist(a,b);
  const double bp = dist(b,p);
  return (ab > ap) && (ab > bp);
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


    
}; // end namespace


#endif
