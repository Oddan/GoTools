#ifndef _TESSELATE_POLYHEDRON_H
#define _TESSELATE_POLYHEDRON_H

#include <vector>
#include <array>
#include "common_defs.h" 
#include "triangulate_domain.h"

namespace TesselateUtils {

  //typedef std::array<double,2> Point2D;
  
  //struct Point2D {double x; double y;};
  
  // struct Polyhedron {

  //   typedef std::pair<int, int> IntPair;
  //   typedef std::vector<IntPair> Loop; // <edge index, orientation flag>
    
  //   std::vector<double> points;
  //   std::vector<IntPair> edges; // indices into 'points'
  //   std::std::vector<Loop> faces: 

  // };

  struct Mesh2D {
    std::vector<Point2D> points;
    std::vector<Triangle> tris;
  };
    
  Mesh2D tesselatePolygon2D(const Point2D* const polygon,
			    const unsigned int num_corners,
			    const double vdist);

  
  std::vector<Point2D> tesselateSegment2D(const Point2D& p1,
					  const Point2D& p2,
					  const double vdist);
				       

  
};

#endif
