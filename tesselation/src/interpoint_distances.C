#include <tuple>
#include <stdexcept>
#include "tesselate_utils.h"
#include "interpoint_distances.h"

using namespace std;
using namespace TesselateUtils;

namespace {
vector<DistanceEntry> 
interpoint_distances_bruteforce_impl(const Point2D* const points,
				     const uint num_points,
				     double R);
vector<DistanceEntry> 
interpoint_distances_smart_impl(const Point2D* const points,
				const uint num_points,
				double R);

const tuple<vector<uint>, vector<uint>> 
bin_points(const Point2D* const points, const uint num_points, double R);

};

namespace TesselateUtils {


// ============================================================================  
vector<DistanceEntry> interpoint_distances(const Point2D* const points,
					   const uint num_points,
					   double R, bool bruteforce)
// ============================================================================
{
  if (bruteforce) {
    return interpoint_distances_bruteforce_impl(points, num_points, R);
  }
  return interpoint_distances_smart_impl(points, num_points, R);
}

}; // end namespace TesselateUtils


namespace {
// ----------------------------------------------------------------------------  
vector<DistanceEntry> 
interpoint_distances_bruteforce_impl(const Point2D* const points, 
				     const uint num_points,
				     double R)
// ----------------------------------------------------------------------------
{
  const double R2 = R*R;
  vector<DistanceEntry> result;
  // brute force implementation
  // @@ This is an N^2 algorithm, so is not expected to scale wells

  const auto points_end = points + num_points;
  for (auto p = points; p != points_end; ++p)
    for (auto q = p+1; q != points_end; ++q)
      if (dist2(*p, *q) < R2)
  	result.push_back({uint(p-points), uint(q-points), dist(*p, *q)});
  
  return result;  
}

// ----------------------------------------------------------------------------
vector<DistanceEntry> 
interpoint_distances_smart_impl(const Point2D* const points, 
				const uint num_points, 
				double R)
// ----------------------------------------------------------------------------
{
  const double R2 = R*R;
  vector<DistanceEntry> result;

  const tuple<vector<uint>, vector<uint>> bins = bin_points(points, num_points, R);
  const auto& bin_x = get<0>(bins);
  const auto& bin_y = get<1>(bins);

  // loop over bins and identify neighbors.  Due to the bin size, neighbors can
  // only be found within a single bin or in adjacent bins.
  const uint maxbin_x = *std::max_element(bin_x.begin(), bin_x.end());
  const uint maxbin_y = *std::max_element(bin_y.begin(), bin_y.end());

  vector<uint> flag(num_points); // will be used to flag points in current bins
  vector<uint> cur_pts; // this vector will hold the indices of the flagged points
  
  for (int iy = 0; iy <= (int)maxbin_y; iy += 2) {
    for (int ix = 0; ix <= (int)maxbin_x; ix += 2) {

      // identify the points in the 3x3 set of bins under study
      transform(bin_x.begin(), bin_x.end(), bin_y.begin(), flag.begin(), 
		[ix, iy](int bx, int by) {
		  return (bx >= ix-1) && (bx <= ix+1) &&
		         (by >= iy-1) && (by <= iy+1);});
      cur_pts.clear();
      for (uint i = 0; i != flag.size(); ++i)
	if (flag[i]) 
	  cur_pts.push_back(i);

      // compare distances between flagged points
      for (uint i = 0; i != cur_pts.size(); ++i) {
	for (uint j = i+1; j != cur_pts.size(); ++j) {
	  auto& p1 = points[cur_pts[i]];
	  auto& p2 = points[cur_pts[j]];
	  if (dist2(p1, p2) < R2)
	    result.push_back( {cur_pts[i], cur_pts[j], dist(p1, p2)});
	}
      }
    }
  }
  sort(result.begin(), result.end());
  result.erase(unique(result.begin(), result.end()), result.end());

  return result;
}

// ----------------------------------------------------------------------------
const tuple<vector<uint>, vector<uint>> 
bin_points(const Point2D* const points, const uint num_points, double R)
// ----------------------------------------------------------------------------
{
  tuple<vector<uint>, vector<uint>> result {vector<uint>(num_points),
                                            vector<uint>(num_points)};
  auto& bin_x = get<0>(result);
  auto& bin_y = get<1>(result);

  const auto bbox = bounding_box(points, num_points);
  const double xmin = bbox[0];
  const double ymin = bbox[2];

  for (uint i = 0; i != num_points; ++i) {
    bin_x[i] = (uint)ceil((points[i][0] - xmin)/R);
    bin_y[i] = (uint)ceil((points[i][1] - ymin)/R);
  }
  return result;
}


};
