#include <assert.h>
#include <algorithm>
#include "GoParametricTesselableVolume.h"

using namespace std;
using namespace Go;

namespace TesselateUtils
{

// ----------------------------------------------------------------------------
// Additional constructor
template<> template<>
GoParametricTesselableVolume::TesselableVolume(ftVolume& fvol)
// ----------------------------------------------------------------------------
{
  const auto shell = fvol.getOuterShell();
  
  // Setting the parametric volume
  volume_ = fvol.getVolume();

  // setting the corners (vertices)
  const auto vertices = fvol.vertices();
  transform(vertices.begin(), vertices.end(), back_inserter(corners_),
            [] (const shared_ptr<Vertex> v) { return v->getVertexPoint();});
  
  // setting the edges and faces
  auto unique_edges = shell->getUniqueInnerEdges();

  shared_ptr<ParamCurve> p = unique_edges[0]->geomCurve();
  
  // setting the faces
}

// ----------------------------------------------------------------------------

template <> array<uint, 2>
GoParametricTesselableVolume::edgeCornerIndices(uint edge_ix) const
// ----------------------------------------------------------------------------  
{
  return {1, 1}; //@@ dummy
}

// ----------------------------------------------------------------------------  
template <> vector<uint>
GoParametricTesselableVolume::faceBoundaryPointIndices(uint face_ix) const
// ----------------------------------------------------------------------------
{
  return {1}; // @@ dummy
}

// ----------------------------------------------------------------------------
template<> 
void GoParametricTesselableVolume::
compute_tesselation(const array<Go::Point, 2>& boundary,
                    const GoParametricTesselableVolume::EdgeType& edge,
                    const double vdist,
                    vector<Point>& ipoints)
// ----------------------------------------------------------------------------
{

}

// ----------------------------------------------------------------------------  
template<>
void GoParametricTesselableVolume::
compute_tesselation(const vector<Point>& boundary,
                    const GoParametricTesselableVolume::FaceType& face,
                    const double vdist,
                    vector<Point>& ipoints,
                    vector<Triangle>& triangles)
// ----------------------------------------------------------------------------  
{

}

// ----------------------------------------------------------------------------  
template<>
void GoParametricTesselableVolume::
compute_tesselation(const vector<Point>& bpoints,
                    const vector<Triangle>& btris,
                    const GoParametricTesselableVolume::VolumeType& volume,
                    const double vdist,
                    vector<Point>& ipoints,
                    vector<Tet>& tets)
// ----------------------------------------------------------------------------  
{

}

  
};
