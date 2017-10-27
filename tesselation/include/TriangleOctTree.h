#ifndef _TRIANGLE_OCT_TREE_H
#define _TRIANGLE_OCT_TREE_H

#include <memory>

#include "common_defs.h"
#include "tesselate_utils.h"

namespace TesselateUtils {

class TriangleOctTree
{

public:
  
  TriangleOctTree(const Point3D* const pts,
                  const Triangle* const tris,
                  const uint num_tris,
                  const uint max_num_ixs = 24);

  // add a triangle to the oct-tree
  void addTriangle(const Triangle& t);

  // Get a list of indices to possible triangles that may intersect triangle
  // 't'.  The result is returned in 'candidate_ixs'.  Its length is equal to
  // the total number of triangles, and entries are 0 or 1, depending of whether
  // the corresponding triangle is an intersection candidate (meaning that an
  // intersection is possible) or not.  If 'clear' is set to 'true', the
  // contents of 'candidates' will be set to zero before searching for
  // candidates (typically what the end user wants).
  void getIntersectionCandidates(const Triangle& t,
                                 std::vector<uint>& candidates,
                                 bool clear = true) const;
  
private:
  const uint max_num_ixs_;
  const Point3D* const points_; // points referred to by triangles
  std::shared_ptr<std::vector<Triangle>> tris_; // triangles (shared with all children)

  // bounding box of the volume enclosed by this OctTree
  const std::array<double, 6> bbox_;
  const std::array<double, 3> midvals_;
  
  // indices to the triangles contained in this volume.  If the volume has been
  // subdivided, this vector will be empty, and children should be queried instead
  std::vector<uint> indices_;

  // if domain is subdivided, the children will be contained here
  std::array<std::shared_ptr<TriangleOctTree>, 8> children_; 

  void reorganize_if_necessary();
  bool is_subdivided() const;
  void include_last_triangle(); // the lastly added triangle should be included
                                // in the volume covered by this OctTree
  void determine_octs(const Triangle& t, std::array<bool, 8>& result) const;
  std::array<bool, 8> determine_octs(const Triangle& t) const;
  void test_integrity(); // for debugging
  
  TriangleOctTree(const Point3D* const pts,
                  std::shared_ptr<std::vector<Triangle>> tris,
                  const std::array<double, 6>& bbox,
                  const uint max_num_ixs);
  
}; // end class TriangleOctTree
  
}; // end namespace

#endif
