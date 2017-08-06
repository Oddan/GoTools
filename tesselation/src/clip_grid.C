#include "clip_grid.h"
#include "tesselate_utils.h"

using namespace std;
using namespace TesselateUtils;

namespace  {

void identify_intersected_domains(const array<double, 4>& bbox,
                                  const array<uint, 2>& res,
                                  const Point2D* const pcorners,
                                  const uint num_pcorners,
                                  vector<ClippedDomainType>& result);
                             
void classify_remaining_cells(const array<double, 4>& bbox,
                            const array<uint, 2>& res,
                            const Point2D* const pcorners,
                            const uint num_pcorners,
                            const double vdist,
                            vector<ClippedDomainType>& result);

void identify_crossings(Point2D p1,
                        Point2D p2,
                        array<double, 4> bbox,
                        array<uint, 2> res,
                        vector<ClippedDomainType>& result,
                        const uint dir);

  
}; //end anonymous namespace

namespace TesselateUtils {

//----------------------------------------------------------------------------
ClippedGrid<2> clip_grid_polygon_2D(const Point2D* const pcorners,
                                    const uint num_pcorners,
                                    const double vdist, 
                                    const uint res_x,
                                    const uint res_y)
//----------------------------------------------------------------------------
{
  const auto bbox = bounding_box_2D(pcorners, num_pcorners);
  ClippedGrid<2> result { bbox,
                          {res_x, res_y},
                          {(bbox[1] - bbox[0])/res_x, (bbox[3] - bbox[2])/res_y},
                          vector<ClippedDomainType>(res_x * res_y, UNDETERMINED)};
  
  // identifying all cells intersected by polygon
  identify_intersected_domains(result.bbox, result.res, pcorners,
                               num_pcorners, result.type);

  // classify non-intersected cells as either 'OUTSIDE', 'FAR_INSIDE' or
  // 'CLOSE_INSIDE'
  classify_remaining_cells(result.bbox, result.res, pcorners,
                           num_pcorners, vdist, result.type);
  
  return result;
};


}; //end namespace TesselateUtils

namespace {

//----------------------------------------------------------------------------  
void identify_intersected_domains(const array<double, 4>& bbox,
                                  const array<uint, 2>& res,
                                  const Point2D* const pcorners,
                                  const uint num_pcorners,
                                  vector<ClippedDomainType>& result)
//----------------------------------------------------------------------------  
{
  for (uint i = 0; i != num_pcorners; ++i) {
    const Point2D& p1 = pcorners[i];
    const Point2D& p2 = pcorners[(i+1) % num_pcorners];

    // identifying vertical and horizontal crossings
    identify_crossings(p1, p2, bbox, res, result, 0); // vertical
    identify_crossings(p1, p2, bbox, res, result, 1);  // horizontal
  }
}

//----------------------------------------------------------------------------
void classify_remaining_cells(const array<double, 4>& bbox,
                              const array<uint, 2>& res,
                              const Point2D* const pcorners,
                              const uint num_pcorners,
                              const double vdist,
                              vector<ClippedDomainType>& result)
//----------------------------------------------------------------------------
{
  const double dx = (bbox[1] - bbox[0]) / (double)res[0];
  const double dy = (bbox[3] - bbox[2]) / (double)res[1];
  for (uint y_ix = 0; y_ix != res[1]; ++y_ix) {
    for (uint x_ix = 0; x_ix != res[0]; ++x_ix) {
      const uint ix = y_ix * res[0] + x_ix;
      if (result[ix] == UNDETERMINED) {
        // This domain is not intersected by the polygon.  Check if it falls outside.
        const Point2D centroid { (0.5 + (double)x_ix) * dx, (0.5 + (double)y_ix) * dy};
        uint dummy = false; //we do not need this value, but required for function call below
        const double r = max(dx, dy);
        result[ix] =
          (!inpolygon(centroid, &pcorners[0], num_pcorners, 0.0, (bool&)dummy)) ? OUTSIDE:
          ((dist2(centroid,
                (closest_point_on_2D_boundary(centroid,
                                              pcorners,
                                              num_pcorners,
                                              dummy))) + r*r) < vdist*vdist) ? CLOSE_INSIDE :
                                                                               FAR_INSIDE;
      }
    }
  }
}
    
//----------------------------------------------------------------------------  
void identify_crossings(Point2D p1,
                        Point2D p2,
                        array<double, 4> bbox,
                        array<uint, 2> res,
                        vector<ClippedDomainType>& result,
                        const uint dir)
//----------------------------------------------------------------------------
{
  assert(dir==0 || dir == 1);
  // ensure the direction we will examine is the first one
  if (dir==1) {
    swap(p1[0], p1[1]);
    swap(p2[0], p2[1]);
    swap(bbox[0], bbox[2]); swap(bbox[1], bbox[3]);
    swap(res[0], res[1]);
  }

  const array<double, 2> cell_len = {(bbox[1] - bbox[0])/res[0],
                                     (bbox[3] - bbox[2])/res[1]};

  // Determine the range of cells along the chosen direction that are
  // intersected by the line segment from p1 to p2
  auto mm = minmax({p1[0], p2[0]});
  const double EPS = sqrt(numeric_limits<double>::epsilon());
  // add small tolerance to avoid missing edge cases
  const double tol = EPS * (mm.second - mm.first);
  mm.first -= tol;
  mm.second += tol;
  const pair<uint, uint> ix_range
     {min((uint)floor(max(mm.first  - bbox[0], 0.0)/cell_len[0]), res[0] - 1),
      min((uint)floor(max(mm.second - bbox[0], 0.0)/cell_len[0]), res[0] - 1)};
  
  // const auto ix_range = minmax({std::min((uint)floor((p1[0] - bbox[0])/cell_len[0]), res[0] - 1),
  //                               std::min((uint)floor((p2[0] - bbox[0])/cell_len[0]), res[0] - 1)});
  
  vector<array<uint,2>> ixs_pairs;
  for (uint i = ix_range.first; i != ix_range.second; ++i) {
    const double pval = (i+1) * cell_len[0];

    // determine at what point, in the other direction than 'dir', does the
    // intersection occur.
    const double t = (pval - p2[0]) / (p1[0] - p2[0]);
    
    const double pval_other = t * p1[1] + (1-t) * p2[1];

    // determine the second index of the cell
    const double pos = (pval_other - bbox[2])/cell_len[1];
    const uint other_ix = min((uint)floor(pos), res[1] - 1);
    
    ixs_pairs.push_back({i, other_ix});

    // handling edge cases
    if (floor(pos-EPS) < other_ix && other_ix > 0) {
      ixs_pairs.push_back( {i, (uint)floor(pos-EPS)});
    } else if (floor(pos+EPS) > other_ix && other_ix < res[1] - 1) {
      ixs_pairs.push_back( {i, (uint)floor(pos+EPS)});
    }
  }

  for (uint i = 0; i != ixs_pairs.size(); ++i) {
    array<uint, 2> ix_pair = ixs_pairs[i]; 
    uint ix = ix_pair[(dir+1)%2] * res[dir] + ix_pair[dir];
    
    result[ix] = INTERSECTED;
    
    ix_pair[0] += 1;
    ix = ix_pair[(dir+1)%2] * res[dir] + ix_pair[dir];

    result[ix] = INTERSECTED;
  }
  // for (uint i = 0; i != ixs_other.size(); ++i) {
  //   array<uint, 2> ix_pair {ix_range.first + i, ixs_other[i]};
  //   uint ix = ix_pair[(dir+1)%2] * res[dir] + ix_pair[dir];
    
  //   result[ix] = INTERSECTED;
    
  //   ix_pair[0] += 1;
  //   ix = ix_pair[(dir+1)%2] * res[dir] + ix_pair[dir];

  //   result[ix] = INTERSECTED;
  // }
}

  
}; // end anonymous namespace 
