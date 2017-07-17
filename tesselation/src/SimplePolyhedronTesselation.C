#include <iostream> // for debugging
#include <array>
#include "tesselate_polyhedron.h"
#include "SimplePolyhedronTesselation.h"
#include "fit_points_to_plane.h"

using namespace std;
using namespace TesselateUtils;

namespace {
  // return all the corners, in order
  vector<uint> assemble_segments(const vector<Segment>& segs);
}; // end anonymous namespace

namespace TesselateUtils
{

// ----------------------------------------------------------------------------  
template<>  
array<uint, 2> SimplePolyhedron::edgeCornerIndices(uint edge_ix) const
// ----------------------------------------------------------------------------  
{
  return edges_[edge_ix];
}

// ----------------------------------------------------------------------------  
template<>  
vector<uint> SimplePolyhedron::faceBoundaryPointIndices(uint face_ix) const
// ----------------------------------------------------------------------------
{
  // first, we add the relevant corner points
  const auto face_edges = faces_[face_ix].edge_ixs;
  vector<uint> result = assemble_segments(extract_from_range(edges_, face_edges));

  // then, we add the indices of the points internal to each edge
  for(const auto& e_ix : face_edges) {
    vector<uint> ixs(edge_ipoints_[e_ix].size(), 0);
    iota(ixs.begin(), ixs.end(), edge_ipoints_start_ixs_[e_ix]);
    result.insert(result.end(), ixs.begin(), ixs.end());
  }
  return result;
}

// ----------------------------------------------------------------------------
template<> 
void SimplePolyhedron::compute_tesselation(const array<Point3D, 2>& boundary,
                                           const Segment& edge,
                                           const double vdist,
                                           vector<Point3D>& ipoints)
// ----------------------------------------------------------------------------
{
  vector<Point3D> points = tesselateSegment(boundary[0], boundary[1], vdist);
  ipoints.resize(points.size()-2);
  copy(points.begin() + 1, points.end()-1, ipoints.begin());
}

// ----------------------------------------------------------------------------
template<> 
void SimplePolyhedron::compute_tesselation(const vector<Point3D>& boundary,
                                           const FaceLoop& face,
                                           const double vdist,
                                           vector<Point3D>& ipoints,
                                           vector<Triangle>& triangles)
// ----------------------------------------------------------------------------
{
  // fit a plane to the points on the boundary
  const array<double, 4> pl = fit_points_to_plane(&boundary[0],
                                                  (uint)boundary.size());

  // check if any point is too far off the fitted plane
  for (const auto& p : boundary) {
    const double dist = p[0] * pl[0] + p[1] * pl[1] + p[2] * pl[2] + pl[3];
    if (dist > vdist/2)
      throw runtime_error("The boundary points of the face are not sufficiently"
                          " close to a common plane.");
  }
  // if we got here, we know that the plane is OK

  // tesselate face in plane

  // lift face back up into 3D
}

// ----------------------------------------------------------------------------  
template<> 
void SimplePolyhedron::compute_tesselation(const vector<Point3D>& bpoints,
                                           const vector<Triangle>& btris,
                                           const VolumeType& volume,
                                           const double vdist,
                                           vector<Point3D>& ipoints,
                                           vector<Tet>& tets)
// ----------------------------------------------------------------------------
{

}
  
}; // end namespace TesselateUtils

namespace {

// ----------------------------------------------------------------------------
// return all the corners, in order
vector<uint> assemble_segments(const vector<Segment>& segs)
// ----------------------------------------------------------------------------
{
  vector<uint> result;
  if (!segs.empty()) {
    result.push_back(segs[0][0]);
    for (const auto s : segs)
      result.push_back(s[0] == result.back() ? s[1] : s[0]);
  }
  return result;
}
    
}; // end anonymous namespace

