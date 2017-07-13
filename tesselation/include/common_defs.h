#ifndef _COMMON_DEFS_H
#define _COMMON_DEFS_H

#include <array>
#include <vector>

namespace TesselateUtils {

using uint = unsigned int;
using Point2D = std::array<double,2>;
struct ValAndDer {
  double val;
  std::vector<Point2D> der;

  void reset(uint num_der) {
    val = 0;
    der = std::vector<Point2D>(num_der, Point2D {0, 0});
  }
};



// ----------------------------------------------------------------------------
inline void operator += (Point2D& p1, const Point2D& p2)
// ----------------------------------------------------------------------------
{ 
  p1[0] += p2[0]; 
  p1[1] += p2[1];
}

// ----------------------------------------------------------------------------
inline void operator -= (Point2D& p1, const Point2D& p2)
// ----------------------------------------------------------------------------
{ 
  p1[0] -= p2[0]; 
  p1[1] -= p2[1];
}

// ----------------------------------------------------------------------------
inline void operator *= (Point2D& p1, const double t)
// ----------------------------------------------------------------------------
{ 
  p1[0] *= t; 
  p1[1] *= t;
}

// ----------------------------------------------------------------------------
inline void operator /= (Point2D& p1, const double t)
// ----------------------------------------------------------------------------
{ 
  p1[0] /= t; 
  p1[1] /= t;
}

// ----------------------------------------------------------------------------
inline Point2D operator + (const Point2D& p1, const Point2D& p2)
// ----------------------------------------------------------------------------
{
  Point2D tmp(p1);
  tmp += p2;
  return tmp;
}

// ----------------------------------------------------------------------------
inline Point2D operator - (const Point2D& p1, const Point2D& p2)
// ----------------------------------------------------------------------------
{
  Point2D tmp(p1);
  tmp -= p2;
  return tmp;
}

// ----------------------------------------------------------------------------
inline Point2D operator*(const Point2D& p, double t)
// ----------------------------------------------------------------------------
{
  Point2D tmp(p);
  tmp *= t;
  return tmp;
}

// ----------------------------------------------------------------------------
inline Point2D operator/(const Point2D& p, double t)
// ----------------------------------------------------------------------------
{
  Point2D tmp(p);
  tmp /= t;
  return tmp;
}



  
}

#endif
