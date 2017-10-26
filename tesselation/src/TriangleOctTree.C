#include "TriangleOctTree.h"
#include <assert.h>

using namespace std;
using namespace TesselateUtils;

namespace {

// ----------------------------------------------------------------------------
inline uint max_ix(const Triangle* const tris, uint num_tris)
// ----------------------------------------------------------------------------
{

  const auto me = max_element(tris, tris + num_tris,
                        [](const Triangle& t1, const Triangle& t2) {
                           return max_element(t1.begin(), t1.end()) <
                                  max_element(t2.begin(), t2.end());});
  return *max_element(me->begin(), me->end());
}

// ----------------------------------------------------------------------------
inline array<double, 3> compute_midvals(const array<double, 6>& bbox)
// ----------------------------------------------------------------------------
{
  array<double, 3> result;
  for (uint i = 0; i != 3; ++i)
    result[i] = 0.5 * (bbox[2*i] + bbox[2*i+1]);
  return result;
}

// ----------------------------------------------------------------------------  
inline array<double, 6> sub_box(const array<double, 6>& bbox, uint oct_ix)
// ----------------------------------------------------------------------------  
{ 
  array<double, 6> result;

  // xmin or xmax?
  if (oct_ix%2 == 0) {
    result[0] = bbox[0];
    result[1] = 0.5 * (bbox[0] + bbox[1]);
  } else { 
    result[0] = 0.5 * (bbox[0] + bbox[1]);
    result[1] = bbox[1];
  }

  // ymin or ymax?
  if ( (oct_ix/2) % 2 == 0) {
    result[2] = bbox[2];
    result[3] = 0.5 * (bbox[2] + bbox[3]);
  } else {
    result[2] = 0.5 * (bbox[2] + bbox[3]);
    result[3] = bbox[3];
  }
  
  // zmin or zmax?
  if ( (oct_ix % 4) == 0) {
    result[4] = bbox[4];
    result[5] = 0.5 * (bbox[4] + bbox[5]);
  } else {
    result[4] = 0.5 * (bbox[4] + bbox[5]);
    result[5] = bbox[5];
  }
  return result;
}
  
}; // end anonymous namespace

namespace TesselateUtils
{

// ============================================================================
TriangleOctTree::TriangleOctTree(const Point3D* const pts,
                                 const Triangle* const tris,
                                 const uint num_tris,
                                 const uint max_num_ixs)
// ============================================================================
  : max_num_ixs_(max_num_ixs), // how many indices the volume can contain before it subdivides
    points_(pts),
    tris_(shared_ptr<vector<Triangle>>(new vector<Triangle>(tris, tris + num_tris))),
    bbox_(bounding_box_3D(pts, max_ix(tris, num_tris))),
    midvals_(compute_midvals(bbox_)), 
    indices_(num_tris)
{
  for (uint i = 0; i != num_tris; ++i)
    indices_[i] = i; // all triangles belong to this oct-tree

  reorganize_if_necessary();
}

// ============================================================================
void TriangleOctTree::addTriangle(const Triangle& t)
// ============================================================================
{
  tris_->push_back(t); // add the triangle itself
  include_last_triangle();  // update the indices_ (of itself or of its children)
  reorganize_if_necessary();
}

// ----------------------------------------------------------------------------
void TriangleOctTree::include_last_triangle()
// ----------------------------------------------------------------------------
{
  if (is_subdivided()) {
    const auto included_octs = determine_octs(tris_->back());
    for (uint i = 0; i != 8; ++i)
      if (included_octs[i]) 
        children_[i]->include_last_triangle();
  } else {
    indices_.push_back(uint(tris_->size()));
    reorganize_if_necessary();
  }
}

// ----------------------------------------------------------------------------  
bool TriangleOctTree::is_subdivided() const
// ----------------------------------------------------------------------------
{
  return (bool) children_[0];
}
  
// ============================================================================
void TriangleOctTree::getIntersectionCandidates(const Triangle& t,
                                                vector<uint>& candidates,
                                                bool clear) const
// ============================================================================
{
  candidates.resize(tris_->size());
  if (clear)
    fill(candidates.begin(), candidates.end(), 0);
  
  if (is_subdivided()) {
    // collect relevant indices from each of the children
    const auto octs = determine_octs(t);
    for (uint i = 0; i != 8; ++i) 
      if (octs[i])
        children_[i]->getIntersectionCandidates(t, candidates, false);
  } else {
    for (auto ix : indices_)
      candidates[ix] = 1;
  }
}

// ----------------------------------------------------------------------------
array<bool, 8> TriangleOctTree::determine_octs(const Triangle& t) const
// ----------------------------------------------------------------------------
{
  array<bool, 8> result;
  determine_octs(t, result);
  return result;
}
  
// ----------------------------------------------------------------------------  
void TriangleOctTree::determine_octs(const Triangle& t,
                                     array<bool, 8>& result) const
// ----------------------------------------------------------------------------
{
  fill(result.begin(), result.end(), true);

  // eliminate octs that are guaranteed not to contain any part of the triangle
  const Point3D& p1 = points_[t[0]];
  const Point3D& p2 = points_[t[1]];
  const Point3D& p3 = points_[t[2]];
  uint sum[3];
  for (uint dim = 0; dim != 3; ++dim) 
    sum[dim] = (p1[dim] < midvals_[dim]) &&
               (p2[dim] < midvals_[dim]) &&
               (p3[dim] < midvals_[dim]);

  for (uint i = 0; i != 8; ++i) {
    const uint xi = i%2;     // zero or one
    const uint yi = (i/2)%2; // zero or one
    const uint zi = (i/4);   // zero or one

    if ((sum[0] == (xi ? 3 : 0)) ||
        (sum[1] == (yi ? 3 : 0)) ||
        (sum[2] == (zi ? 3 : 0)))
      result[i] == false;
  }
}
  
// ----------------------------------------------------------------------------
void TriangleOctTree::reorganize_if_necessary()
// ----------------------------------------------------------------------------
{
  if ((uint)indices_.size() > max_num_ixs_) {
    assert(!is_subdivided()); // indices_ should be empty otherwise

    // construct children
    for (uint i = 0; i != 8; ++i) 
      children_[i] = shared_ptr<TriangleOctTree>
        (new TriangleOctTree(points_, tris_, sub_box(bbox_, i), max_num_ixs_));
    
    // distribute triangles
    array<bool, 8> octs;
    for (const auto ix : indices_) {
      determine_octs((*tris_)[ix], octs);
      for (uint i = 0; i != 8; ++i) 
        if (octs[i])
          (children_[i]->indices_).push_back(ix);
    }
  }
}

// ----------------------------------------------------------------------------
TriangleOctTree::TriangleOctTree(const Point3D* const pts,
                                 shared_ptr<vector<Triangle>> tris,
                                 const std::array<double, 6>& bbox,
                                 const uint max_num_ixs)
// ----------------------------------------------------------------------------
  : max_num_ixs_(max_num_ixs),
    points_(pts),
    tris_(tris),
    bbox_(bbox),
    midvals_(compute_midvals(bbox_))
{}

 
  
}; // end namespace TesselateUtils
