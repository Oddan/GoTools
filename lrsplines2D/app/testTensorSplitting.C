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
  vector<int> missing_tensorgrid_segments(const Mesh2D& m,
					  const Direction2D d,
					  const IntPair range_x,  // x-range
					  const IntPair range_y); // y-range
  
  // missing _inner_ refinements (i.e. excluding the boundary)
  int num_missing_increments(const Mesh2D& m,
			     Direction2D d, 
			     IntPair c1,  // lower left corner
			     IntPair c2); // upper right corner

  //vector<int> active_knots(const Mesh2D& m, Direction2D d, IntPair r1, IntPair r2);
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
  vector<int> n = missing_tensorgrid_segments(m, XFIXED, {0, 6}, {0, 4}); // (0,6), (0,5)

  copy(n.begin(), n.end(), ostream_iterator<int, wchar_t>(wcout, L" "));
  wcout<< endl;
  
  //const auto knots = active_knots(m, XFIXED, {0, 6}, {3, 5});
  //const auto knots = active_knots(m, XFIXED, {4, 6}, {3, 5});
  //copy(knots.begin(), knots.end(), ostream_iterator<int, wchar_t>(wcout, L" \n"));
    
}

// ============================================================================

namespace {
  // --------------------------------------------------------------------------
  vector<int> active_knots(const Mesh2D& m, Direction2D d, IntPair r1, IntPair r2)
  // --------------------------------------------------------------------------
  {
    // The first and last knot is automatically 'active' by merit to belong to
    // the boundary of the domain under consideration
    
    if (d==YFIXED) swap(r1, r2); // ensure that r1 correspond to the fixed direction
    vector<int> result {r1.first}; // include first knot (see initial comment)

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
  vector<int> make_segmult_map(const vector<GPos>& segs, int start_ix, int end_ix)
  // --------------------------------------------------------------------------
  {
    // determine the knot multiplicities of minisegments starting with a given knot
    vector<int> result(end_ix - start_ix, segs[0].mult);

    for (int i = 1; i != segs.size(); ++i) {

      const int ix = segs[i].ix - start_ix;
      if (ix > 0 && ix < result.size())
	fill(result.begin()+ix, result.end(), segs[i].mult);
    }
    return result;
  }

  
  // --------------------------------------------------------------------------
  vector<int> minisegment_multiplicities(const Mesh2D& m,
					 const Direction2D d,
					 const int seg_ix,
					 const vector<int> knots)
  // --------------------------------------------------------------------------
  {
    assert(knots.size() >= 2);

    const auto all_segs = m.mrects(d, seg_ix); // all segments

    // segments intersecting with the intervals within the given knots

    // making map of multiplicities
    const vector<int> segmult = make_segmult_map(all_segs, knots.front(), knots.back());

    // computing the multiplicities of each segment between two consecutive knots
    vector<int> result;
    result.reserve(knots.size()-1);
    transform(knots.begin(), knots.end()-1, back_inserter(result),
	      [&] (int k) {return segmult[k - knots.front()];});
		  
    return result;
  }


  // ==========================================================================
  vector<int> missing_tensorgrid_segments(const Mesh2D& m,
					  const Direction2D d,
					  const IntPair range_x, 
					  const IntPair range_y) 
  // ==========================================================================
  {
    // determine active knots in the opposite direction (not all zero in the segment)
    const auto d2_knots = active_knots(m, flip(d), range_x, range_y); 

    const IntPair r1  = (d == XFIXED) ? range_x : range_y; // d-direction
    const IntPair r2 = (d == XFIXED) ? range_y : range_x;  // other direction

    vector<int> result;
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
  int num_missing_increments(const Mesh2D& m,
			      Direction2D d, 
			      IntPair c1,  // lower left corner
			      IntPair c2)  // upper right corner
  // ==========================================================================
  {
    
    return {};
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
} // end anonymous namespace

  
