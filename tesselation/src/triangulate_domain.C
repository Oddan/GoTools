#include <algorithm>
#include <assert.h>
#include <stdexcept>
#include "triangulate_domain.h"
#include "tesselate_utils.h"

using namespace std;
using namespace TesselateUtils;

namespace {
  Triangle add_triangle(vector<Segment>& nsegs,   // input-output
			vector<Segment>& dsegs,   // input-output
			vector<uint>& unused_pts,    // input-output
			const Point2D* const points, // input only
			const double vdist);         // input only

  void find_candidate_points(const Segment& s,
			     const vector<Segment>& nsegs,
			     const vector<uint>& unused_pts,
			     const Point2D* const points,
			     const double vdist,
			     vector<uint>& cand_pts,         // input-output
			     vector<uint>& all_neigh_pts);    // input-output

  uint best_candidate_point(const vector<uint>& cand_pts, // nb: will be modified
			    const Segment& seg,
			    const Point2D* const points);

  bool is_delaunay(const Triangle& tri,
		   const vector<uint>& neigh_pts,
		   const Point2D* const points);

  bool add_or_remove_segment(const Segment& seg,
			     vector<Segment>& nsegs,
			     vector<Segment>& dsegs,
			     vector<uint>& unused_pts,
			     const bool delaunay,
			     const bool orientation); // 1: the new point is first in the segment
                                                      // 0: the new point is last in the segment 
  bool introduces_intersection(const Segment& s,
			       const uint pt_ix,
			       const vector<Segment>& nsegs,
			       const Point2D* const points,
			       const double tol);

  bool search_and_erase_segment(const Segment& seg, vector<Segment>& segvec);
  
};

namespace TesselateUtils {

// ============================================================================
// The algorithm used here is inspired by:
// S.H. Lo, "Delaunay Triangulation of Non-Convex Planar Domains" (1989)
vector<Triangle> triangulate_domain(const Point2D* const points,
				    const uint num_bpoints, 
				    const uint tot_num_points,
				    const double vdist)
// ============================================================================  
{

  // Preparing vectors keeping track of "non-delaunay" and delaunay segments.
  // The working front constists of both types of segments, but while adding
  // triangles we only need to check for intersections against the non-delaunay
  // segments.  At the beginning, the working front is equal to the boundary
  // polgygon, which is considered to consist solely of non-delaunay segments.
  // (If we were positive the boundary polygon was convex, we could probably
  // consider them delaunay segments instead).
  vector<Segment> nsegs(num_bpoints), dsegs; 

  // At start, all boundary segments are "non-delaunay".  We store them all.
  for (uint i = 0; i != num_bpoints; ++i)
    nsegs[i] = {i, (i+1)%num_bpoints};

  // We use this vector to keep track of indices to nodes that are still
  // interior or on the working front, which in the beginning means "all nodes".
  vector<uint> unused_pts(tot_num_points, 1);  // 1 for unused, 0 for used

  // Establish the result vector, and gradually fill it with rectangles until
  // there are no remaining points on the working front or in the interior.
  // NB: in the function call in the loop below, the three first arguments are
  // input-output (i.e. they will be modified by the function).
  vector<Triangle> result;

  while (nsegs.size() + dsegs.size() > 0)
    result.push_back(add_triangle(nsegs, dsegs, unused_pts, points, vdist));

  // sanity check: there should be no unused nodes left by now
  assert(accumulate(unused_pts.begin(), unused_pts.end(), 0) == 0);
  
  return result;
}

  
}; // end namespace TesselateUtils

namespace {

// ----------------------------------------------------------------------------
Triangle add_triangle(vector<Segment>& nsegs,   
		      vector<Segment>& dsegs,   
		      vector<uint>& unused_pts,    
		      const Point2D* const points, 
		      const double vdist)
// ----------------------------------------------------------------------------
{
  // We choose the next segment to work with.  We aim to deplete the
  // non-delaunay segments as fast as possible, since we do not have to do
  // intersection checks againt the delaunay segments.  We therefore pick a
  // segment from 'nsegs' as long as there are any left.
  const Segment cur_seg = (nsegs.size() > 0) ? nsegs.back() : dsegs.back();
  (nsegs.size() > 0) ? nsegs.pop_back() : dsegs.pop_back();
  
  // Identify all points within 'vdist' of the segment.  We separate these into
  // 'candidate' and 'non-candidate' points.  Candidate points must fulfill the 
  // additional criteria:
  // (1) They have to be on the 'left side' of the current segment
  // (2) No intersections with other segments on on the working front must occur
  //     when connecting the candidate point to either endpoint of the current
  //     segment (sufficient to check agaisnt the non-delaunay segments)
  // We will chose the third point of the triangle from the candidate points.
  // We must still check against all neighbor points to determine whether
  // our triangle is delaunay or not (semi-delaunay is still "not delaunay")
  vector<uint> cand_pts, all_neigh_pts;
  find_candidate_points(cur_seg, nsegs, unused_pts, points, vdist,
			cand_pts, all_neigh_pts);

  // Choosing the best candidate point, i.e. a candidate points for which no
  // other candidate points are within the circumscribing circle of the triangle
  // generated if using this point.  
  const uint chosen_pt = best_candidate_point(cand_pts, cur_seg, points);

  // Determine whether the generated triangle is a "pure" delaunay triangle, or
  // if it contains other points caused by the nonconvexity of the working front
  // (or if it is semi-delaunay by having other points exactly on its
  // circumscription).
  const bool is_del = is_delaunay({cur_seg[0], cur_seg[1], chosen_pt},
				  all_neigh_pts, points);
  
  // Modify working front and remaining nodes
  const bool removed1 = add_or_remove_segment({cur_seg[0], chosen_pt}, nsegs, dsegs,
					      unused_pts, is_del, false);
  const bool removed2 = add_or_remove_segment({chosen_pt, cur_seg[1]}, nsegs, dsegs,
					      unused_pts, is_del, true);
  if (removed1 && removed2) {
    // both segments already existed.  The triangle just filled in a triangular
    // hole in the mesh.  This means the chosen point will not be on the working
    // front and we should remove it from the active list as well
    unused_pts[chosen_pt] = 0;
  }
  return Triangle {cur_seg[0], cur_seg[1], chosen_pt};
}

// ----------------------------------------------------------------------------
void find_candidate_points(const Segment& s,
                           const vector<Segment>& nsegs,
                           const vector<uint>& unused_pts,
                           const Point2D* const points,
                           const double vdist,
                           vector<uint>& cand_pts,         // input-output
                           vector<uint>& all_neigh_pts)     // input-output
// ----------------------------------------------------------------------------
{
  const double TOL = 1.0e-6;  // @@ is this safe/general enough?
  const uint N = (uint)unused_pts.size(); // total number of points
  cand_pts.resize(N);     fill(cand_pts.begin(), cand_pts.end(), 0);
  all_neigh_pts.resize(N); fill(all_neigh_pts.begin(), all_neigh_pts.end(), 0);

  for(uint i = 0; i != N; ++i) {
    if ((!unused_pts[i]) || (i== s[0]) || (i == s[1])) {
      // this is either a non-active point or one of the segment endpoints; we
      // need not consider those
      continue;
    }
    if (point_on_line_segment(points[i], points[s[0]], points[s[1]], vdist)) {
      // point is close enough to be a neighbor point
      all_neigh_pts[i] = 1;
      // check if it is also a candidate point
      if ((projected_distance_to_line(points[i], points[s[0]], points[s[1]]) < 0) &&
	  (!introduces_intersection(s, i, nsegs, points, TOL))) {
	cand_pts[i] = 1; // point is a candidate point
      } 
    }
  }
}

// ----------------------------------------------------------------------------  
bool introduces_intersection(const Segment& s,
			     const uint pt_ix,
			     const vector<Segment>& nsegs,
			     const Point2D* const points,
			     const double tol)
// ----------------------------------------------------------------------------  
{
  // @@ NB: there are algorithms for efficiently intersecting one set of line
  //    segments against another in one go.  If the below, "one-by-one"
  //    algorithm becomes a bottleneck, consider switching to the more efficient
  //    set-based approach.
  for (const auto& ns: nsegs)
    if ((ns[0] == pt_ix) || ns[1] == pt_ix)
      // candidate point is one of the endpoints of this segment.  No
      // new intersection possible.
      continue;
    else if ((ns[0] == s[0] && ns[1] == s[1]) || (ns[0] == s[1] && ns[1] == s[0]))
      // the two segments are the same - no new intersection created
      continue;
    else if (segments_intersect(points[pt_ix], points[s[0]],
				points[ns[0]], points[ns[1]], -tol) ||
	     segments_intersect(points[pt_ix], points[s[1]],
				points[ns[0]], points[ns[1]], -tol))
      return true;
  return false;
}


// ----------------------------------------------------------------------------  
uint best_candidate_point(const vector<uint>& cand_pts,
                          const Segment& seg,
                          const Point2D* const points)
// ----------------------------------------------------------------------------
{
  // find indices of candidate points
  vector<uint> cand_ixs;
  for (uint i = 0; i != cand_pts.size(); ++i)
    if (cand_pts[i])
      cand_ixs.push_back(i);

  if (cand_ixs.empty()) 
    // @@ this ought to be handled more gracefully
    throw runtime_error("No candidate found, increase 'vdist' and try again.");

  // Searching for the candidate whose resulting triangle has a circumscription
  // that does not contain any other candidate
  Point2D center;  // center of the circumscribed circle
  double radius2;   // squared radius of the circumscribed circle
  const Point2D& p1 = points[seg[0]];
  const Point2D& p2 = points[seg[1]];

  // @@ It may be faster to eliminate neighbors as one go along (any neighbour
  // outside the circumscribed circle could be disregarded and removed from the
  // vector).  This introduces a bit more logic, but something to consider if
  // this part of the code proves to be a bottleneck.
  for (const auto ix : cand_ixs)  {
    if (circumscribe_triangle(points[ix], p1, p2, center, radius2)) {
      // check whether the circle contains any candidate point
      bool interior_point_found = false;
      for (const auto ix2 : cand_ixs) 
	if (dist2(center, points[ix2]) < radius2) {
	  interior_point_found = true;
	  break;
	}
      if (!interior_point_found)
	return ix;
    }
  }
  // we should never get here by design, so throw an error if it happens
  throw runtime_error("No suitable best candidate.  This should not happen.  "
		      "Please investigate.");
  return -1; // dummy, to keep compiler happy
}
  
// ----------------------------------------------------------------------------
bool is_delaunay(const Triangle& tri, 
                 const vector<uint>& neigh_pts,
                 const Point2D* const points)
// ----------------------------------------------------------------------------
{
  Point2D circ_center;
  double radius2;
  const double tol = 1e-5; //@@ OK?  In the worst case, we risk classifying a
			   //delaunay edge as non-delaunay, which would not do
			   //any harm to the final result, but slightly increase
			   //the computational cost (since there is now one more
			   //segment for which we need to check intersections).
  circumscribe_triangle(points[tri[0]], points[tri[1]], points[tri[2]],
			circ_center, radius2);

  radius2 *= (1+tol);  // slightly increase the radius to capture points on the
		      // edge of the circle.  If we capture points slightly
		      // outside, that should still be ok, see comment above.
  
  const uint N = (uint)neigh_pts.size();
  for (uint i = 0; i != N; ++i) 
    if ((neigh_pts[i]) &&                               // should be neigh poitn
	((i!=tri[0]) && (i!=tri[1]) && (i!=tri[2])) &&  // should not be corner of triangle
	(dist2(points[i], circ_center) < radius2))
      return false;
              
  return true;
}

// ----------------------------------------------------------------------------  
bool add_or_remove_segment(const Segment& seg,
                           vector<Segment>& nsegs,
                           vector<Segment>& dsegs,
                           vector<uint>& unused_pts,
                           const bool delaunay,
			   const bool orientation)
// ----------------------------------------------------------------------------
{
  // first, determine if segment already exists
  bool segment_found = search_and_erase_segment(seg, nsegs);
  if (!segment_found) {
    segment_found = search_and_erase_segment(seg, dsegs);
  }

  if (segment_found) {
    // segment already existed and has been removed (since both its triangles
    // has now been located).  We must remove the point that is no longer on the
    // working boundary
    unused_pts[orientation ? seg[1] : seg[0]] = 0;
    return true;
  } else {
    // this is a newly introduced segment.  Add it to the respective list of
    // segments
    if (delaunay)
      dsegs.push_back(seg);
    else
      nsegs.push_back(seg);
  }
  return false;
}

// ----------------------------------------------------------------------------
bool search_and_erase_segment(const Segment& seg, vector<Segment>& segvec)
// ----------------------------------------------------------------------------
{
  auto iter = find_if(segvec.begin(), segvec.end(), [&seg](const Segment& s) {
      return (((s[0]==seg[0]) && (s[1]==seg[1])) ||
	      ((s[0]==seg[1]) && (s[1]==seg[0])));
    });

  if (iter != segvec.end()) {
    segvec.erase(iter);
    return true;
  }
  return false; 
}
  
}; // end anonymous namespace

