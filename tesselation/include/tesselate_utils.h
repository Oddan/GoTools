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
double projected_distance_to_line(const P& p, const P& a, const P& b);
// ----------------------------------------------------------------------------    


// ----------------------------------------------------------------------------    
// Check if the point p is within distance 'tol' of the line segment defined by
// points a and b
template<typename P>
bool point_on_line_segment(const P& p, const P& a, const P& b, double tol);
// ----------------------------------------------------------------------------    

// ----------------------------------------------------------------------------    
double random_uniform(double minval, double maxval);
// ----------------------------------------------------------------------------    

}; // end namespace TesselateUtils

#include "tesselate_utils_impl.h"

#endif
