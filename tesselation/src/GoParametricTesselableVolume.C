#include <assert.h>
#include <algorithm>
#include "GoParametricTesselableVolume.h"
#include "tesselate_parametric_volume.h"

using namespace std;
using namespace Go;

namespace {

// ----------------------------------------------------------------------------
uint find_vertex_index(const vector<shared_ptr<Vertex>>& vertices,
                       const shared_ptr<Vertex> vertex)
// ----------------------------------------------------------------------------
{
  const auto ptr = find(vertices.begin(), vertices.end(), vertex);
  assert(ptr != vertices.end()); // it should be in there somewhere
  return uint(ptr - vertices.begin());
}

// ----------------------------------------------------------------------------  
vector<pair<uint, bool>>
construct_edge_index_vector(const vector<shared_ptr<ftEdge>>& edges,
                            const shared_ptr<Loop> loop)
// ----------------------------------------------------------------------------  
{
  vector<pair<uint, bool>> result;
  const auto loop_edges = loop->getEdges(); // vector of shared_ptr<ftEdgeBase>
  for (auto e : loop_edges) {
    pair<uint, bool> cur_elem;

    // find the current edge in the edge list
    bool found = false;
    bool found_was_twin = false;
    for (uint i = 0; (i != edges.size()) && !found; ++i) {
      if (edges[i] == e) {
        cur_elem.first = i;
        found = true;
        found_was_twin = false;
      } else if (edges[i]->twin() == e.get()) {
        cur_elem.first = i;
        found = true;
        found_was_twin = true;
      }
    }
    assert(found);

    // determine orientation
    const auto ft_e = e->geomEdge();
    bool is_reversed = found_was_twin ^ ft_e->isReversed();
    cur_elem.second = !is_reversed;

    result.push_back(cur_elem);
  }
    
  return result;
}
                            
  
}; // end anonymous namespace

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
  const auto vertices = fvol.vertices(); // shared_ptrs to vertices
  transform(vertices.begin(), vertices.end(), back_inserter(corners_),
            [] (const shared_ptr<Vertex> v) { return v->getVertexPoint();});
  
  // setting the edges 
  auto unique_edges = shell->getUniqueInnerEdges();
  transform(unique_edges.begin(), unique_edges.end(), back_inserter(edges_),
            [&vertices] (const shared_ptr<ftEdge> e) {
              // @@ for the moment, we require that the ParamCurve is not clipped
              assert(e->tMin() == e->geomCurve()->startparam());
              assert(e->tMax() == e->geomCurve()->endparam());
                
              return EdgeType { e->geomCurve(),
                                find_vertex_index(vertices, e->getVertex(true)), 
                                find_vertex_index(vertices, e->getVertex(false))};});
      
  // setting the faces.  We assume here that there are no inner loops
  auto ft_faces = shell->allFaces();
  transform(ft_faces.begin(), ft_faces.end(), back_inserter(faces_),
            [&unique_edges] (const shared_ptr<ftSurface> fs) {
              //assert(fs->nmbBoundaryLoops() == 1); // not yet support for inner loops
              return FaceType {
                fs->surface(),
                construct_edge_index_vector(unique_edges,
                                            fs->getBoundaryLoop(0))};});
  
}

// ----------------------------------------------------------------------------

template <> array<uint, 2>
GoParametricTesselableVolume::edgeCornerIndices(uint edge_ix) const
// ----------------------------------------------------------------------------  
{
  const auto e = edges_[edge_ix];
  return {get<1>(e), get<2>(e)};
}

// ----------------------------------------------------------------------------  
template <> vector<uint>
GoParametricTesselableVolume::faceBoundaryPointIndices(uint face_ix) const
// ----------------------------------------------------------------------------
{
  const auto edge_ixs = faces_[face_ix].ix;
  vector<uint> result;

  for (auto e : edge_ixs) {
    const auto edge = edges_[e.first];

    // internal points
    vector<uint> ixs(edge_ipoints_[e.first].size(), 0);
    iota(ixs.begin(), ixs.end(), edge_ipoints_start_ixs_[e.first]);

    
    if (e.second) { // oriented in the correct way

      result.push_back(get<1>(edge));

    } else { // reversely oriented

      result.push_back(get<2>(edge));
      reverse(ixs.begin(), ixs.end());
    }
    
    result.insert(result.end(), ixs.begin(), ixs.end());
  }

  return result;  
}

// ----------------------------------------------------------------------------
template<> 
void GoParametricTesselableVolume::
compute_tesselation(const array<Go::Point, 2>& boundary,
                    const EdgeType& edge,
                    const double vdist,
                    vector<Point>& ipoints)
// ----------------------------------------------------------------------------
{
  const auto param = tesselateParametricCurve(get<0>(edge), vdist);
  ipoints.resize(param.size() - 2);
  transform(param.begin()+1, param.end()-1, ipoints.begin(),
            [&edge] (const double d) { return get<0>(edge)->point(d); });
}

// ----------------------------------------------------------------------------  
template<>
void GoParametricTesselableVolume::
compute_tesselation(const vector<Point>& boundary,
                    const FaceType& face,
                    const double vdist,
                    vector<Point>& ipoints,
                    vector<Triangle>& triangles)
// ----------------------------------------------------------------------------  
{
  // @@ dummy
  ipoints = boundary;
}

// ----------------------------------------------------------------------------  
template<>
void GoParametricTesselableVolume::
compute_tesselation(const vector<Point>& bpoints,
                    const vector<Triangle>& btris,
                    const VolumeType& volume,
                    const double vdist,
                    vector<Point>& ipoints,
                    vector<Tet>& tets)
// ----------------------------------------------------------------------------  
{
  // @@ dummy
  ipoints = bpoints;
}

  
};
