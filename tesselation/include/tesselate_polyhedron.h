#ifndef _TESSELATE_POLYHEDRON_H
#define _TESSELATE_POLYHEDRON_H

#include <vector>
#include <array>
#include <cmath>
#include "common_defs.h"
#include "tesselate_utils.h"
#include "triangulate_domain.h"

namespace TesselateUtils {

struct Mesh2D {
  std::vector<Point2D> points;
  std::vector<Triangle> tris;
};

struct Mesh3D {
  std::vector<Point3D> points;
  std::vector<Tet> tets;
};

template<typename PointXD> 
std::vector<PointXD> tesselateSegment(const PointXD& p1,
					const PointXD& p2,
					const double vdist);

// Boundary polygon should be given counterclockwise
Mesh2D tesselatePolygon2D(const Point2D* const polygon,
			  const unsigned int num_corners,
			  const double vdist,
			  const bool tesselate_boundary=true);

// All  face triangle should be presented with face normal pointing outwards.
Mesh3D tesselatePolyhedron3D(const Point3D* const bpoints,
                             const unsigned int num_bpoints,
                             const Triangle* const btris,
                             const unsigned int num_tris,
                             const double vdist);


// ========================= TEMPLATE IMPLEMENTATIONS =========================

template<typename PointXD> inline
std::vector<PointXD> tesselateSegment(const PointXD& p1,
				      const PointXD& p2,
				      const double vdist)
{
  const double seg_len = dist(p1, p2);
  const unsigned int num_intervals = (unsigned int)std::ceil(seg_len/vdist);

  return interpolate(p1, p2, num_intervals - 1);
}
};






#endif
