#ifndef _TESSELATE_UTILS_H
#define _TESSELATE_UTILS_H

#include <array>
#include <vector>
#include <algorithm>

namespace TesselateUtils {

  const double PI = 3.14159265358979323;

template<typename P> inline double dist2(const P& p1, const P& p2)
// ----------------------------------------------------------------------------
{
  const double dx = p1[0] - p2[0];
  const double dy = p1[1] - p2[1];
  return dx*dx + dy*dy;
}

// ----------------------------------------------------------------------------
template<typename P> inline double dist(const P& p1, const P& p2) 
// ----------------------------------------------------------------------------
{
  return sqrt(dist2(p1, p2));
}

// ----------------------------------------------------------------------------
template<typename P>
std::vector<P> interpolate(const P& p1, const P& p2, unsigned int num);
// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------
// Returns signed area of polygon.  A counterclockwise-ordered polygon has
// positive area.
template<typename P>
double polygon_area(const P* const poly, const unsigned int num_corners);
// ----------------------------------------------------------------------------  


// ----------------------------------------------------------------------------
// Returns a bounding box on the form [xmin, xmax, ymin, ymax]
template<typename P>  
std::array<double, 4> bounding_box(const P* const points,
				   const unsigned int num_points);
// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------
// Generate grid of points with resolution (nx+2) x (ny+2) in the rectangle
// defined by the two opposite corners c1 and c2.  nx and ny denote the number
// of _internal_ grid points, hence the (nx+2) x (ny+2) resolution
template<typename P>
std::vector<P> generate_grid_2D(const P& c1, const P& c2, int nx, int ny);
// ----------------------------------------------------------------------------  

// ----------------------------------------------------------------------------  
template<typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& arg);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------
// Check if a point is inside a polygon
template<typename P> bool inpolygon(const P& pt, const P* const poly, 
				    const unsigned int num_corners, const double tol);
// ----------------------------------------------------------------------------      
  
// ----------------------------------------------------------------------------
// Return the points that are inside the polygon
template<typename P>
std::vector<P> inpolygon(const P* const pts, const unsigned int num_pts,
			 const P* const poly, const unsigned int num_corners, 
			 const double tol);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------    
// 'p' is the point we want to compute distance from.  'a' and 'b' are two
// distinct points on the infinite line (thus defining the line)
template<typename P>
double projected_distance_to_line(P p, P a, P b);
// ----------------------------------------------------------------------------    


// ----------------------------------------------------------------------------    
// Check if the point p is within distance 'tol' of the line segment defined by
// points a and b
template<typename P> inline
bool point_on_line_segment(const P& p, const P& a, const P& b, double tol);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------
// Returns 'true' if the line through point 'p' that perpendicularly intersects
// the line passing through 'a' and 'b' has its itnersection point within the
// segment ab.
template<typename P> inline
bool projects_to_segment(const P& p, const P& a, const P&b);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------
// return 'true' if the angle abc is acute
template<typename P> inline
bool acute_angle(const P& a, const P& b, const P& c);
// ----------------------------------------------------------------------------      

// ----------------------------------------------------------------------------
// Mirror point 'p' accross the line passing through 'a' and 'b', and return the
// result.
template<typename P> inline
P mirror_point(P p, const P& a, const P& b);
// ----------------------------------------------------------------------------
  
  
// ----------------------------------------------------------------------------
// mirror the provided set of points across the line passing through a and b,
// and return the result (the mirror points) as a vector
template<typename P> inline
std::vector<P> mirror_points(const P* const pts,
			     const unsigned int num_pts,
			     const P& a,
			     const P& b);
// ----------------------------------------------------------------------------

  
// ----------------------------------------------------------------------------    
double random_uniform(double minval, double maxval);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------
// similar to MATLAB's "find".  There are two versions, depending on whether
// clarity of code (result in return value) or speed (result as output argument)
// is favored.
template<typename T> inline
std::vector<unsigned int> locate_nonzeros(const std::vector<T> vec);

template<typename T> inline
void locate_nonzeros(const std::vector<T> vec, std::vector<unsigned int>& target);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------
// Extract all elements from a range that satisfy a predicate. The predicate
// should take a reference to T and return a bool (or int).  The first vector in
// the result contains the extracted values, the second vector contain the
// respective indices.
template<typename T, typename Pred> inline
std::pair<std::vector<T>, std::vector<unsigned int>>
extract_from_range(const T* const range_start,
		   unsigned int range_length, const Pred& fun);
// ----------------------------------------------------------------------------    

  
}; // end namespace TesselateUtils

#include "tesselate_utils_impl.h"

#endif
