#include <fstream> // @@ debug purposes
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

  Tet add_tet(vector<Triangle>& ntris, // input-output
              vector<Triangle>& dtris, // input-output
              vector<uint>& unused_pts,
              const Point3D* const points,
              const double vdist);
  
  void find_candidate_points(const Segment& s,
			     const vector<Segment>& nsegs,
			     const vector<uint>& unused_pts,
			     const Point2D* const points,
			     const double vdist,
			     vector<uint>& cand_pts,         // input-output
			     vector<uint>& all_neigh_pts);    // input-output

  void find_candidate_points(const Triangle& tri,
                             const vector<Triangle>& ntris,
                             const vector<uint>& unused_pts,
                             const Point3D* const points,
                             const double vdist,
			     vector<uint>& cand_pts,         // input-output
			     vector<uint>& all_neigh_pts);    // input-output

  uint best_candidate_point(const vector<uint>& cand_pts, // nb: will be modified
			    const Segment& seg,
			    const Point2D* const points);

  uint best_candidate_point(const vector<uint>& cand_pts, // nb: will be modified
                            const Triangle& cur_tri,
                            const Point3D* const points);

  bool is_delaunay(const Triangle& tri,
		   const vector<uint>& neigh_pts,
		   const Point2D* const points);

  bool is_delaunay(const Triangle& tri, const uint chosen_pt,
                   const vector<uint>& neigh_pts, const Point3D* const points);
  
  bool add_or_remove_segment(const Segment& seg,
			     vector<Segment>& nsegs,
			     vector<Segment>& dsegs,
			     vector<uint>& unused_pts,
			     const bool delaunay,
			     const bool orientation); // 1: the new point is first in the segment
                                                      // 0: the new point is last in the segment

  bool add_or_remove_face(const Triangle& tri, // third corner of 'tri' is the new point
                          vector<Triangle>& ntris,
                          vector<Triangle>& dtris,
                          vector<uint>& unused_pts,
                          const bool delaunay);
  
  bool introduces_intersection(const Segment& s,
			       const uint pt_ix,
			       const vector<Segment>& nsegs,
			       const Point2D* const points,
			       const double tol);

  bool introduces_intersection(const Triangle& tri,
                               const uint pt_ix,
                               const vector<Triangle>& ntris,
                               const Point3D* const points,
                               const double tol);

  bool search_and_erase_segment(const Segment& seg, vector<Segment>& segvec);
  bool search_and_erase_face(const Triangle& tri, vector<Triangle>& trivec);

  
  
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

// ============================================================================
// The algorithm used here is a 3D generalization og the algorithm used above,
// and which was inspired by S.H. Lo, "Delaunay Triangulation of Non-Convex
// Planar Domains" (1989)
vector<Tet> construct_tets(const Point3D* const points,
                           const uint tot_num_points,
                           const Triangle* btris,
                           const uint num_btris,
                           const double vdist)
// ============================================================================
{
  // preparing vectors keeping track of 'non-delaunay' and 'delaunay' faces
  // (terms used as the 3D generalization of terms used for triangulate_domain()).
  // At start, all boundary triangles are considered "non-delaunay".
  vector<Triangle> ntris(btris, btris + num_btris), dtris;

  // We use this vector to keep track of indices to nodes that are still
  // interior or on the working front, which in the beginning means "all nodes".
  vector<uint> unused_pts(tot_num_points, 1); // 1 for unused, 0 for used

  // Establish the result vector, and gradually fill it with tets until there
  // are no remaining points on the workign front nor in the interior.
  vector<Tet> result;

  while (ntris.size() + dtris.size() > 0)
    result.push_back(add_tet(ntris, dtris, unused_pts, points, vdist));

  // sanity check: there should be no unused nodes left by now
  assert(accumulate(unused_pts.begin(), unused_pts.end(), 0) == 0);

  return result;
}

  
}; // end namespace TesselateUtils

namespace {

// ----------------------------------------------------------------------------
Tet add_tet(vector<Triangle>& ntris,
            vector<Triangle>& dtris,
            vector<uint>& unused_pts,
            const Point3D* const points,
            const double vdist)
// ----------------------------------------------------------------------------
{
  // // @@ write current status to file
  // ofstream points_os("points.mat");
  // ofstream unused_os("unused.mat");
  // for (uint i = 0; i != (uint)unused_pts.size(); ++i) {
  //   points_os << points[i];
  //   unused_os << unused_pts[i] << " ";
  // }
  // points_os.close();
  // unused_os.close();
  // ofstream tris_os("triangles.mat");
  // for (uint i = 0; i != ntris.size(); ++i)
  //   tris_os << ntris[i] << '\n';
  // for (uint i = 0; i != dtris.size(); ++i)
  //   tris_os << dtris[i] << '\n';
  // tris_os.close();

  // Choose the next triangle to work with, aiming to deplete the non-delaunay
  // triangles as quickly as possible
  const Triangle cur_tri = (ntris.size() > 0) ? ntris.back() : dtris.back();
  (ntris.size() > 0) ? ntris.pop_back() : dtris.pop_back();

  // Identify all points within 'vdist' of triangle, and separate them into
  // 'candidate' and 'non-candidate' points, according to criteria analogous to
  // those mentioned in the comments to 'add_triangle'.  
  vector<uint> cand_pts, all_neigh_pts;
  find_candidate_points(cur_tri, ntris, unused_pts, points, vdist,
                        cand_pts, all_neigh_pts);

  // The 'best' of the candidate points will be chosen as the fourth corner of
  // the tet.  The best point is the one with no other candidate points within
  // the circumscribing sphere of the generated tet.
  const uint chosen_pt = best_candidate_point(cand_pts, cur_tri, points);

  // Determine whether teh generated tet is a "delaunay tet" or if it contains
  // other points caused by nonconvexity of the working front (or if it is
  // semi-delaunay by having other points exactly on its circumscription
  const bool is_del = is_delaunay(cur_tri, chosen_pt, all_neigh_pts, points);
  
  // modify working front and remaining nodes.  Make sure to add them in
  // _clockwise_ corner order, since we want them to be pointing outwards of the
  // active front.
  add_or_remove_face({cur_tri[0], cur_tri[1], chosen_pt}, ntris, dtris, unused_pts, is_del);
  add_or_remove_face({cur_tri[1], cur_tri[2], chosen_pt}, ntris, dtris, unused_pts, is_del);
  add_or_remove_face({cur_tri[2], cur_tri[0], chosen_pt}, ntris, dtris, unused_pts, is_del);

  // Verify the four corners of the inserted tet - how many of them remain on
  // the active front?  We set them to 0, and add back those that are still
  // found on the front.
  unused_pts[chosen_pt] = 0;
  for (uint i = 0; i != 3; ++i) unused_pts[cur_tri[i]] = 0;
  for(const auto t : ntris) { for (uint i = 0; i != 3; ++i) unused_pts[t[i]] = 1; }
  for(const auto t : dtris) { for (uint i = 0; i != 3; ++i) unused_pts[t[i]] = 1; }

  return Tet {cur_tri[0], cur_tri[1], cur_tri[2], chosen_pt};
}


  
// ----------------------------------------------------------------------------
uint best_candidate_point(const vector<uint>& cand_pts,
                          const Triangle& cur_tri,
                          const Point3D* const points)
// ----------------------------------------------------------------------------
{
  // find indices of candidate points
  vector<uint> cand_ixs;
  for (uint i = 0; i != (uint)cand_pts.size(); ++i)
    if (cand_pts[i])
      cand_ixs.push_back(i);

  if (cand_ixs.empty())
    // @@ this should be handled more gracefully
    throw runtime_error("No candidate found for new triangle point."
                        "  Increase 'vdist' and try again.");

  // searching for the candidate whose resulting tet has a minimal containing
  // sphere that does not contain any other candidate.
  Point3D center; // center of the minimal containing sphere
  double radius2; // squared radius of the containing sphere

  // @@ It may be faster to eliminate neighbors as one goes along (any neighbour
  // otuside the sphere could be immediately disregarded and removed from the
  // vector).  This would however introduce a bit of extra logic.
  for (const auto ix : cand_ixs) {
    if (fitting_sphere(points[ix], points[cur_tri[0]], points[cur_tri[1]],
                       points[cur_tri[2]], center, radius2)) {
      // check whether the sphere contains any candidate point
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
  // we should never get her by design, so throw an error if it happens
  throw runtime_error("No suitable best candidate.  This should not happen.");
  return -1; // dummy, to keep compiler happy
}
  
// ----------------------------------------------------------------------------  
void find_candidate_points(const Triangle& tri,
                           const vector<Triangle>& ntris,
                           const vector<uint>& unused_pts,
                           const Point3D* const points,
                           const double vdist,
                           vector<uint>& cand_pts,         // input-output
                           vector<uint>& all_neigh_pts)    // input-output
// ----------------------------------------------------------------------------
{
  const double TOL = 1.0e-6; // @@ safe/general enough?
  const uint N = (uint)unused_pts.size();
  const Point3D& p1 = points[tri[0]]; // p1, p2 and p3 references are for 
  const Point3D& p2 = points[tri[1]]; // convenience only
  const Point3D& p3 = points[tri[2]];
  
  cand_pts.resize(N); fill(cand_pts.begin(), cand_pts.end(), 0);
  all_neigh_pts.resize(N); fill(all_neigh_pts.begin(), all_neigh_pts.end(), 0);
  
  for (uint i = 0; i != N; ++i) {
    if ((!unused_pts[i]) || (i== tri[0]) || (i == tri[1]) || (i == tri[2])) {
      // this is either a non-active point or one of the triangle corners; we
      // need not consider those.
      continue;
    }
    if (point_on_triangle(points[i], p1, p2, p3, vdist)) {
      // point is close enough to triangle to be considered a neighbor point.
      all_neigh_pts[i] = 1;
      if ((point_on_inside_of_face(points[i], p1, p2, p3)) &&
          (!introduces_intersection(tri, i, ntris, points, TOL))) {
        cand_pts[i] = 1; // point is a candidate point
      }
    }
  }
}

// ----------------------------------------------------------------------------
bool introduces_intersection(const Triangle& tri,
                             const uint pt_ix,
                             const vector<Triangle>& ntris,
                             const Point3D* const points,
                             const double tol)
// ----------------------------------------------------------------------------
{
  // checking each non-delaunay triangle for intersection against the triangles
  // that would be introduced if the point indiced by 'pt_ix' were to form a tet
  // with the already-existing triangle 'tri'.
  for (const auto& nt : ntris) 
    if ((nt[0] == pt_ix) || (nt[1] == pt_ix) || (nt[2] == pt_ix))
      // candidate point is one of the corners of this triangle, no new
      // intersection possible
      continue;
    else if (nt == tri)
      // the two triangles are actually the same - no new intersection created
      continue;
    else if (triangles_intersect_3D(points[tri[0]], points[tri[1]], points[pt_ix],
                                    points[nt[0]], points[nt[1]], points[nt[2]], -tol) ||
             triangles_intersect_3D(points[tri[1]], points[tri[2]], points[pt_ix],
                                    points[nt[0]], points[nt[1]], points[nt[2]], -tol) ||
             triangles_intersect_3D(points[tri[2]], points[tri[0]], points[pt_ix],
                                    points[nt[0]], points[nt[1]], points[nt[2]], -tol))
      return true;
  return false;
}
  
  
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
    if (point_on_line_segment(points[i], points[s[0]], points[s[1]], vdist, false)) {
      // point is close enough to be a neighbor point
      all_neigh_pts[i] = 1;
      // check if it is also a candidate point
      if ((projected_distance_to_line_2D(points[i], points[s[0]], points[s[1]]) < 0) &&
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
    else if (segments_intersect_2D(points[pt_ix], points[s[0]],
				   points[ns[0]], points[ns[1]], -tol) ||
	     segments_intersect_2D(points[pt_ix], points[s[1]],
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
  for (uint i = 0; i != (uint)cand_pts.size(); ++i)
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

  // @@ It may be faster to eliminate neighbors as one goes along (any neighbour
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
    if ((neigh_pts[i]) &&                               // should be neigh point
	((i!=tri[0]) && (i!=tri[1]) && (i!=tri[2])) &&  // should not be corner of triangle
	(dist2(points[i], circ_center) < radius2))
      return false;
              
  return true;
}

// ----------------------------------------------------------------------------
bool is_delaunay(const Triangle& tri, const uint chosen_pt,
                 const vector<uint>& neigh_pts, const Point3D* const points)
// ----------------------------------------------------------------------------
{
  Point3D center;
  double radius2;
  const double tol = 1e-5; // same comment as for the 2D version of 'is_delaunay above
  fitting_sphere(points[tri[0]], points[tri[1]], points[tri[2]], points[chosen_pt],
                 center, radius2);
  radius2 *= (1+tol); // slightly increase radius to ensure points on the boundary are captured

  const uint N = (uint)neigh_pts.size();
  for (uint i = 0; i != N; ++i)
    if ((neigh_pts[i]) && // should be a neighbor point
        (i != chosen_pt) && // should not be the chosen point
        ((i!=tri[0]) && (i!=tri[1]) && (i!=tri[2])) && // should not be part of triangle
        (dist2(points[i], center) < radius2)) // should be within the radius
      return false;  // this point is inside the sphere, hence the proposed 'tet' is not delaunay

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
    // have now been located).  We must remove the point that is no longer on the
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

// ----------------------------------------------------------------------------
bool add_or_remove_face(const Triangle& tri, // third corner of 'tri' is the new point
                        vector<Triangle>& ntris,
                        vector<Triangle>& dtris,
                        vector<uint>& unused_pts,
                        const bool delaunay)
// ----------------------------------------------------------------------------
{
  const bool tri_found = search_and_erase_face(tri, ntris) ||
                         search_and_erase_face(tri, dtris);
  if (tri_found) 
    // triangle already existed and has been removed from the front.
    return true;

  // this is a newly introduced triangle.  Add it to the respective list of
  // triangles on the boundary
  (delaunay ? dtris : ntris).push_back(tri);
  return false;
}

// ----------------------------------------------------------------------------
bool search_and_erase_face(const Triangle& tri, vector<Triangle>& trivec)
// ----------------------------------------------------------------------------
{
  auto iter = find_if(trivec.begin(), trivec.end(), [&tri](Triangle t) {
      reverse(t.begin(), t.end());
      // flip triangle, rotate it three times, and check for a match each time
      for (uint i = 0; i!=3; ++i) {
        if (equal(t.begin(), t.end(), tri.begin()))
          return true;
        rotate(t.begin(), t.begin()+1,t.end());
      }
      return false;
    });
  if (iter != trivec.end()) {
    trivec.erase(iter);
    return true;
  }
  return false;
}

}; // end anonymous namespace

