#ifndef _INTERPOINT_DISTANCES_H
#define _INTERPOINT_DISTANCES_H

#include <vector>
#include "common_defs.h"

namespace TesselateUtils {

struct DistanceEntry {
  uint p1_ix; 
  uint p2_ix; 
  double dist;
  
  bool operator==(const DistanceEntry& rhs) const {
    return (p1_ix == rhs.p1_ix) &&
           (p2_ix == rhs.p2_ix) &&
           (dist  == rhs.dist );
  }

  bool operator<(const DistanceEntry& rhs) const {
    return (p1_ix < rhs.p1_ix) ? true :
           (p1_ix > rhs.p1_ix) ? false :
           (p2_ix < rhs.p2_ix) ? true : false;
  }
};

// Compute distances between points.  Only consider points whose distances are
// closer than R.  Each entry in the result vector contains the indices to a
// pair of points, as well as their mutual distances.  (This result can
// conceptually be considered as a sparse matrix.)
std::vector<DistanceEntry> interpoint_distances(const Point2D* points,
						const uint num_points,
						const double R,
						const bool bruteforce=false);
  
// Compute distances from points to a specific, separate point 'p'.  Only
// consider points whose distances to 'p' are shorter than R.  The second index
// in the DistanceEntry will be set to 0.
std::vector<DistanceEntry> interpoint_distances(const Point2D* points,
						const uint num_points,
						const Point2D& p,
						const double R);

// Compute interpoint distances between points in 'points' and points in
// 'other_points'.  Only consider distances shorter than R. If the same point is
// found in both sets (as identified by its memory location), it will ignore
// the distance between the point and itself.
std::vector<DistanceEntry> interpoint_distances(const Point2D* points,
						const uint num_points,
						const Point2D* other_points,
						const uint num_other_points,
						const double R);
  
}; //end namespace TesselateUtils


#endif
