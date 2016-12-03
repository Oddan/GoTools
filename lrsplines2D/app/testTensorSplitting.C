#include <locale>
#include <iostream>
#include <iterator> // for ostream_iterator
#include <vector>
//#include <tuple>
//#include <utility>
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/PlotUtils.h"

using namespace Go;
using namespace std;

using IntPair = pair<int, int>;
using IntVec  = vector<int>;

namespace {
  int insert_segment(Mesh2D& m,
		     Direction2D d,
		     double pval,
		     int start_ix,
		     int end_ix,
		     int mult);

  // determine the number of segments/multiplicities that are missing in order
  // to convert the specified subdomain to a full grid.   Does not consider the
  // boundary segments (first and last knot of each range)
  IntVec missing_tensorgrid_segments(const Mesh2D& m,
				     const Direction2D d,
				     const IntPair range_x,  // x-range
				     const IntPair range_y); // y-range

  // count total number of missing segments/multiplicities to make the specified
  // subdomain into a full grid (not counting boundaries).
  struct MissingSegInfo {int num; IntVec missing_x; IntVec missing_y;};
  MissingSegInfo total_missing_tensorgrid_segments(const Mesh2D&m,
						   const IntPair range_x,
						   const IntPair range_y);

  // Suggest a split of the specified domain into two subdomains with fewer
  // total missing segments/multiplicities.  'bnd_mult' specifies the
  // multiplicity each segment of the split itself should have.
  // Return values:
  // 1) Direction (whether the split occurs for XFIXED or YFIXED)
  // 2) Index of the line along which the split should occur (or -1 if no split
  //    was found)
  // 3) The total missing tensorgrid segments before split
  // 4) The cost of the split (the number of minisegments multiplicities to
  //    raise along the split line)
  // 5) The total missing tensorgrid segments of the first new subdomain
  // 6) The total missing tensorgrid segments of the second new subdomain
  // If a split is found, then the sum of the three last returned values should
  // be strictly lower than the value of the third returned value.  If this is
  // not possible, no split is proposed, and the second returned value will be -1.
  //tuple<Direction2D, int, int, int, int, int>
  struct SplitInfo {
    Direction2D d;
    int ix;
    int init_cost;
    int split_cost;
    int new_cost1;
    int new_cost2;
  };
  SplitInfo suggest_domain_split(const Mesh2D& m, const IntPair range_x,
				 const IntPair range_y, const int bnd_mult);
} // end anonymous namespace

// ============================================================================

int main(int argc, char* argv[])
{
  setlocale(LC_ALL, "en_US.UTF.8");
  
  // Defining initial knotvectors
  const double xvec[] {0, 1, 2, 3, 4, 5};
  const double yvec[] {0, 1, 2, 3, 4};

  // base, tensor-product mesh
  Mesh2D m = Mesh2D(xvec, xvec + 6, yvec, yvec + 5);

  // inserting modifications (no longer tensor product)
  insert_segment(m, YFIXED, 2.5, 2, 5, 1);
  insert_segment(m, XFIXED, 3.5, 1, 3, 1);

  plot_mesh(m);

  // TESTING
  IntVec n = missing_tensorgrid_segments(m, YFIXED, {0, 2}, {0, 5}); // (0,6), (0,5)

  copy(n.begin(), n.end(), ostream_iterator<int, wchar_t>(wcout, L" "));
  wcout<< endl;

  int total = total_missing_tensorgrid_segments(m, {0, 6}, {0, 5}).num;

  wcout << total << endl;
  
  
  //const auto knots = active_knots(m, XFIXED, {0, 6}, {3, 5});
  //const auto knots = active_knots(m, XFIXED, {4, 6}, {3, 5});
  //copy(knots.begin(), knots.end(), ostream_iterator<int, wchar_t>(wcout, L" \n"));
    
}

// ============================================================================

namespace {
  // --------------------------------------------------------------------------
  IntVec active_knots(const Mesh2D& m, Direction2D d, IntPair r1, IntPair r2)
  // --------------------------------------------------------------------------
  {
    // The first and last knot is automatically 'active' by merit to belong to
    // the boundary of the domain under consideration
    
    if (d==YFIXED) swap(r1, r2); // ensure that r1 correspond to the fixed direction
    IntVec result {r1.first}; // include first knot (see initial comment)

    // determine which of the "interior" knots should be considered active
    for (int i = r1.first+1; i < r1.second; ++i) {
      const vector<IntPair> segs = m.segments(d, i);
      if (any_of(segs.begin(), segs.end(),
		 [r2] (IntPair p) {return (p.first < r2.second) &
		                          (p.second > r2.first);}))
	result.push_back(i);
    }

    result.push_back(r1.second); // include last knot (see inital comment)

    return result;
  }

  // --------------------------------------------------------------------------
  IntVec make_segmult_map(const vector<GPos>& segs, int start_ix, int end_ix)
  // --------------------------------------------------------------------------
  {
    // determine the knot multiplicities of minisegments starting with a given knot
    IntVec result(end_ix - start_ix, segs[0].mult);

    for (int i = 1; i != segs.size(); ++i) {

      const int ix = segs[i].ix - start_ix;
      if (ix > 0 && ix < result.size())
	fill(result.begin()+ix, result.end(), segs[i].mult);
    }
    return result;
  }

  
  // --------------------------------------------------------------------------
  IntVec minisegment_multiplicities(const Mesh2D& m,
					 const Direction2D d,
					 const int seg_ix,
					 const IntVec knots)
  // --------------------------------------------------------------------------
  {
    assert(knots.size() >= 2);

    const auto all_segs = m.mrects(d, seg_ix); // all segments

    // segments intersecting with the intervals within the given knots

    // making map of multiplicities
    const IntVec segmult = make_segmult_map(all_segs, knots.front(), knots.back());

    // computing the multiplicities of each segment between two consecutive knots
    IntVec result;
    result.reserve(knots.size()-1);
    transform(knots.begin(), knots.end()-1, back_inserter(result),
	      [&] (int k) {return segmult[k - knots.front()];});
		  
    return result;
  }


  // ==========================================================================
  IntVec missing_tensorgrid_segments(const Mesh2D& m,
				     const Direction2D d,
				     const IntPair range_x, 
					  const IntPair range_y) 
  // ==========================================================================
  {
    // determine active knots in the opposite direction (not all zero in the segment)
    const auto d2_knots = active_knots(m, flip(d), range_x, range_y); 

    const IntPair r1  = (d == XFIXED) ? range_x : range_y; // d-direction
    const IntPair r2 = (d == XFIXED) ? range_y : range_x;  // other direction

    IntVec result;
    result.reserve(max(r1.second - r1.first - 1, 1));
    
    // we only consider interior.  Ignore first and last value of range
    for (int i = r1.first+1; i != r1.second; ++i) {

      // multiplicities for each "mini-segment" along this line, within the domain
      const auto ms_mult = minisegment_multiplicities(m, d, i, d2_knots);

      const int max_mult = *max_element(ms_mult.begin(), ms_mult.end());

      // compute number of missing minisegments, and add it to result vector
      result.push_back(max_mult * ms_mult.size() -
		       accumulate(ms_mult.begin(), ms_mult.end(), 0) );
    }
      
    return result;
  }

  // ==========================================================================
  MissingSegInfo total_missing_tensorgrid_segments(const Mesh2D&m,
						   const IntPair range_x,
						   const IntPair range_y)
  // ==========================================================================
  {
    const auto c1 = missing_tensorgrid_segments(m, XFIXED, range_x, range_y);
    const auto c2 = missing_tensorgrid_segments(m, YFIXED, range_x, range_y);

    return {
      accumulate(c1.begin(), c1.end(), 0) + accumulate(c2.begin(), c2.end(), 0),
      c1,
      c2};
  }
  

  // ==========================================================================
  int insert_segment(Mesh2D& m,
		    Direction2D d,
		    double pval,
		    int start_ix,
		    int end_ix,
		    int mult)
  // ==========================================================================
  {
    int l_ix = m.insertLine(d, pval);
    m.setMult(d, l_ix, start_ix, end_ix, mult);
    return l_ix;
  }

  // ==========================================================================
  SplitInfo suggest_domain_split(const Mesh2D& m,
				 const IntPair range_x,
				 const IntPair range_y,
				 const int bnd_mult)
  // ==========================================================================
  {
    //struct MissingSegInfo {int num; IntVec missing_x; IntVec missing_y;};
    const auto init = total_missing_tensorgrid_segments(m, range_x, range_y);
    const int missing_x = accumulate(init.missing_x.begin(), init.missing_x.end(), 0);
    const int missing_y = accumulate(init.missing_y.begin(), init.missing_y.end(), 0);


    // Choose direction of split (the direction with fewest missing multiplicities
    const Direction2D d = (missing_x < missing_y) ? XFIXED : YFIXED; 
    const IntPair range_d     = (d == XFIXED) ? range_x : range_y;
    const IntPair range_other = (d == XFIXED) ? range_y : range_x;

    // Identify candidates for split (only those associated with t-junctions in
    // the relevant range)
    IntVec cand(range_d.second - range_d.first + 1, 0);
    iota(cand.begin(), cand.end(), range_d.first); // consecutive values
    cand.erase(remove_if(cand.begin(), cand.end(),[&m, d, range_other](int ix)
			 {!has_t_junctions(m, d, ix, range_other);}),
	       cand.end());

    // Search for optimal candidate
    IntVec perf; // candidate performance
    transform(cand.begin(), cand.end(), back_inserter(perf), [&] (int ix) {
	__implement_me__;
      });

    // identify the winning candidate, and check if the corresponding split is
    // useful (positive performance)
    const auto max_perf = max_element(perf.begin(), perf.end());
    const int split_ix = (*max_perf > 0) ? cand(max_perf - perf.begin()) : -1;

    // TODO
    const int split_cost = 1;
    const int new_cost1 = 1;
    const int new_cost2 = 1;

    return {d, split_ix, init.num, split_cost, new_cost1, new_cost2};
  }
    
  
} // end anonymous namespace

  
