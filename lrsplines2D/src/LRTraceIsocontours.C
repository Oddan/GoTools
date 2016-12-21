#include <fstream> // @@ debug purpose
#include <functional>
#include "GoTools/lrsplines2D/LRTraceIsocontours.h"
#include "GoTools/lrsplines2D/SSurfTraceIsocontours.h"

using namespace std;
using namespace Go;

using LRSurfPtr = shared_ptr<LRSplineSurface>;
using IsectCurve = pair<CurvePtr, CurvePtr> ; // @@
namespace {

// ----------------------------------------------------------------------------
// Wrapper to simplify syntax when the shared pointer itself is not required
// except to manage the allocated memory
inline shared_ptr<SplineSurface> as_spline_surf(const shared_ptr<LRSplineSurface> lrs)
// ----------------------------------------------------------------------------
{
  return shared_ptr<SplineSurface>(lrs->asSplineSurface());
}

// ----------------------------------------------------------------------------  
vector<CurveVec> merge_isocontours(const vector<vector<CurveVec>>& curve_fragments,
				   const vector<LRSurfPtr>& surf_frags,
				   const double tol);
// ----------------------------------------------------------------------------

  
};// end anonymous namespace


namespace Go
{

// ============================================================================
vector<CurveVec> LRTraceIsocontours(const LRSplineSurface& lrs,
				    const std::vector<double>& isovals,
				    const double tol,
				    const bool include_3D_curves,
				    const bool use_sisl_marching)
// ============================================================================
{
  //break LR surface up into individual patches
  const vector<LRSurfPtr> surf_fragments = lrs.subdivideIntoSimpler();

  // Setting up function to compute isocontours for a surface fragment
  const function<vector<CurveVec>(LRSurfPtr)> compute_isovals {
    [&] (LRSurfPtr l){return SSurfTraceIsocontours(*as_spline_surf(l), isovals,
						   tol, include_3D_curves,
						   use_sisl_marching);}};

  // computing isocurves for each surface fragment (vector<vector<CurveVec>>)
  const auto curve_fragments = apply_transform(surf_fragments, compute_isovals);

  // // debug
  // ofstream os_surf("lrsurf.g2");
  // for (auto f : surf_fragments) {
  //   f->writeStandardHeader(os_surf);
  //   f->write(os_surf);
  // }
  // os_surf.close();

  // ofstream os_cv("curvefrags.g2");
  // for (auto frag : curve_fragments)
  //   for (auto ival : frag)
  //     for (auto c : ival) {
  // 	c.second->writeStandardHeader(os_cv);
  // 	c.second->write(os_cv);
  //     }
  // os_cv.close();
  
  // merge isocontours across patches and returning result
  return merge_isocontours(curve_fragments, surf_fragments,  tol);
}

}; // end namespace Go;

namespace {

// ----------------------------------------------------------------------------
  array<double, 4> parameter_domain(const LRSurfPtr& patch)
// ----------------------------------------------------------------------------
{
  return array<double, 4>{ { patch->paramMin(XFIXED),
	                     patch->paramMax(XFIXED), 
	                     patch->paramMin(YFIXED),
	                     patch->paramMax(YFIXED)} };
}

// ----------------------------------------------------------------------------
  int find_exit_edge(const IsectCurve& c, const array<double, 4>& domain,
		     const bool start, const double tol)
// ----------------------------------------------------------------------------
{
  const double t = start ? c.first->startparam() : c.first->endparam();
  Point uv;  // represent u and v parameters
  c.first->point(uv, t);  // evaluate (u, v) parameter pair at curve start or end

  const vector<double> dists = {fabs(domain[0] - uv[0]),
				fabs(domain[1] - uv[0]),
				fabs(domain[2] - uv[1]),
				fabs(domain[3] - uv[1])};
  auto min_it = min_element(dists.begin(), dists.end());
  if (*min_it > tol)
    return -1; // curve doesn't end on edge, but inside domain (either closed, or ends at
	       // singularity)

  return int(min_it - dists.begin()); // index of boundary edge where curve starts/ends
		 
}

// ----------------------------------------------------------------------------
void map_curve(const IsectCurve& c, const array<double, 4>& domain, const double tol,
	       const bool at_start, map<double, CurveVec>& u_map, map<double, CurveVec>& v_map)
// ----------------------------------------------------------------------------
{
  const int edge = find_exit_edge(c, domain, at_start, tol);

  if (edge >=0) {
    Point p1; c.first->point(p1, at_start ? c.first->startparam() : c.first->endparam());
    assert( (fabs(domain[edge] - p1[0]) < tol) | (fabs(domain[edge] - p1[1]) < tol) );
  }

  // if none of the cases below apply, curve does not terminate at edge.
  switch (edge) {
  case 0:
  case 1:
    u_map[domain[edge]].push_back(c);
    break;
  case 2:
  case 3:
    v_map[domain[edge]].push_back(c);
    break;
  }
}    

  // ----------------------------------------------------------------------------
void prepare_curvemaps(const array<double, 4>& domain, const double tol, const CurveVec& cvec,
		       map<double, CurveVec>& u_map, map<double, CurveVec>& v_map)
// ----------------------------------------------------------------------------
{
  for_each(cvec.begin(), cvec.end(), [&] (const IsectCurve& c) {
      map_curve(c, domain, tol, true, u_map, v_map);  // map curve according to start point
      map_curve(c, domain, tol, false, u_map, v_map); // map curve according to end point
  });
}

// ----------------------------------------------------------------------------
IsectCurve join_isectcurves(const IsectCurve& c1, const IsectCurve& c2,
			    bool c1_at_start, bool c2_at_start)
// ----------------------------------------------------------------------------
{
  shared_ptr<SplineCurve> pcurve1 = shared_ptr<SplineCurve>(c1.first->clone());
  shared_ptr<SplineCurve> scurve1 = shared_ptr<SplineCurve>(c1.second->clone());
  shared_ptr<SplineCurve> pcurve2 = shared_ptr<SplineCurve>(c2.first->clone());
  shared_ptr<SplineCurve> scurve2 = shared_ptr<SplineCurve>(c2.second->clone());

  if (c1_at_start) {
    pcurve1->reverseParameterDirection();
    scurve1->reverseParameterDirection();
  }

  if (!c2_at_start) {
    pcurve2->reverseParameterDirection();
    scurve2->reverseParameterDirection();
  }
  
  pcurve1->appendCurve(pcurve2.get());
  scurve1->appendCurve(scurve2.get());
  
  return IsectCurve { pcurve1, scurve1 };
}

// ----------------------------------------------------------------------------
void replace_segments(const IsectCurve& old1, const IsectCurve& old2,
		      const IsectCurve& updated, map<double, CurveVec>& target)
// ----------------------------------------------------------------------------
{
  for (auto& it_map : target) 
    for (auto& it_vec : it_map.second) 
      if ((it_vec.first == old1.first) || (it_vec.first == old2.first)) 
	it_vec = updated;
}


// ----------------------------------------------------------------------------
pair<double, double> identify_truncated_endpoints(const IsectCurve& c, Direction2D d,
						  double pval, const double tol)
// ----------------------------------------------------------------------------
{
  Point startpoint, endpoint;
  c.first->point(startpoint, c.first->startparam());
  c.first->point(endpoint  , c.first->endparam());

  const int ix = (d==YFIXED) ? 1 : 0;
  auto NaN = numeric_limits<double>::quiet_NaN();
  return pair<double, double> {
    (fabs(startpoint[ix] - pval) < tol) ? startpoint[(ix+1)%2] : NaN,
    (fabs(endpoint[ix] - pval) < tol) ? endpoint[(ix+1)%2] : NaN
  };
  
  // if (d==YFIXED) {
  //   swap(startpoint[0], startpoint[1]);
  //   swap(endpoint[0], endpoint[1]);
  // }
  // auto NaN = numeric_limits<double>::quiet_NaN();
  // return pair<double, double> {
  //   (fabs(startpoint[0] - pval) < tol ? startpoint[1] : NaN),
  //   (fabs(endpoint[0] - pval) < tol ? endpoint[1] : NaN)
  // };
}


// ----------------------------------------------------------------------------
void merge_segments(map<double, CurveVec>& mergemap, // map whose segments should be merged
		    map<double, CurveVec>& othermap, // map whose entries must also be updated
		    Direction2D d,                   // the fixed parameter direction
		    const double tol,               
		    CurveVec& bcurves,            
		    CurveVec& finished_curves) // insert newly merged finished curves here
// ----------------------------------------------------------------------------
{
  struct EndPoint {double pval; IsectCurve icurve;bool at_start;};

  while (!mergemap.empty()) {
    auto it = *mergemap.begin();   mergemap.erase(mergemap.begin());
    // sorting incident curve segments so that those that should be merged will lie right next
    // to each other
    vector<CurvePtr> encountered;
    vector<EndPoint> tp_vec;
    for (auto& ic : it.second) { // loop over intersection curve segments cut by this parameter line
      const auto ends = identify_truncated_endpoints(ic, d, it.first, tol);
      if (!isnan(ends.first + ends.second)) {
	// both endpoints of this curve are truncated by the parameter line.
	if (find(encountered.begin(), encountered.end(), ic.first) != encountered.end())
	  continue;
	encountered.push_back(ic.first); // ensure we will only register the curve once
      }
      if (!isnan(ends.first))
	tp_vec.push_back({ends.first, ic, true});
      if (!isnan(ends.second)) 
	tp_vec.push_back({ends.second, ic, false});
    }

    assert(tp_vec.size() % 2 == 0);  // we suppose that for each curve going in, there is one going out

    sort(tp_vec.begin(), tp_vec.end(), [](const EndPoint& t1, const EndPoint& t2) {return t1.pval < t2.pval;});

    for (size_t i = 0; i != tp_vec.size(); i += 2) {
      const auto entry1 = tp_vec[i];
      const auto entry2 = tp_vec[i+1];

      assert(fabs(entry1.pval - entry2.pval) < tol);
      if (entry1.icurve.first == entry2.icurve.first) {
	// this is a single curve whose endpoints meet across this edge.  There are no more
	// merges to be done.  The curve is finished, and can be returned.
	finished_curves.push_back(entry1.icurve);
      } else {
	auto new_curve = join_isectcurves(entry1.icurve, entry2.icurve, entry1.at_start, entry2.at_start);

	// replace references to the old curves with references to new_curve throughout
	replace_segments(entry1.icurve, entry2.icurve, new_curve, mergemap);
	replace_segments(entry1.icurve, entry2.icurve, new_curve, othermap);

	for (size_t j = i+2; j < tp_vec.size(); ++j) {
	  if ((tp_vec[j].icurve.first == entry1.icurve.first) |
	      (tp_vec[j].icurve.first == entry2.icurve.first)) {
	    tp_vec[j].icurve = new_curve;
	    tp_vec[j].at_start =
	      *(new_curve.first->coefs_begin() + (d==XFIXED ? 1 : 0)) == tp_vec[j].pval;
	  }
	}

	// updating boundary curve pointers if necessary
	transform(bcurves.begin(), bcurves.end(), bcurves.begin(), [&](const IsectCurve& c) {
	    return ((c.first == entry1.icurve.first) | (c.first == entry2.icurve.first)) ? new_curve : c;});
      }
    }
  }
}

// ----------------------------------------------------------------------------
array<double, 4> outer_boundary_pvals(const vector<LRSurfPtr>& patches)
// ----------------------------------------------------------------------------
{
  const auto start_u = [](LRSurfPtr l) {return l->startparam_u();};
  const auto end_u   = [](LRSurfPtr l) {return l->endparam_u();};
  const auto start_v = [](LRSurfPtr l) {return l->startparam_v();};
  const auto end_v   = [](LRSurfPtr l) {return l->endparam_v();};
  
  return array<double, 4>
    {(*min_element(patches.begin(), patches.end(),
	   [&start_u](LRSurfPtr l1, LRSurfPtr l2) {return start_u(l1) < start_u(l2);}))->startparam_u(),
     (*max_element(patches.begin(), patches.end(),
	   [&end_u](LRSurfPtr l1, LRSurfPtr l2) {return end_u(l1) < end_u(l2);}))->endparam_u(),
     (*min_element(patches.begin(), patches.end(),
	   [&start_v](LRSurfPtr l1, LRSurfPtr l2) {return start_v(l1) < start_v(l2);}))->startparam_v(),
     (*max_element(patches.begin(), patches.end(),
	   [&end_v](LRSurfPtr l1, LRSurfPtr l2) {return end_v(l1) < end_v(l2);}))->endparam_v()};
}

// ----------------------------------------------------------------------------
template<typename T>
void add_to_vec(vector<T>& target, const vector<T>& new_elements)
// ----------------------------------------------------------------------------
{
  target.insert(target.end(), new_elements.cbegin(), new_elements.cend());
}

// ----------------------------------------------------------------------------
template<typename T>
vector<T> expand_vec(const vector<vector<T>>& vv)
// ----------------------------------------------------------------------------
{
  vector<T> result;
  for (auto v : vv)
    add_to_vec(result, v);
  return result; 
}

// ----------------------------------------------------------------------------
template<typename K, typename T>
map<K, vector<T>> map_combine(const map<K, vector<T>>& m1,
			      const map<K, vector<T>>& m2)
// ----------------------------------------------------------------------------
{
  auto result = m1;
  result.insert(m2.begin(), m2.end());
  return result;
}
  
// ----------------------------------------------------------------------------
template<typename K, typename T>
vector<T> expand_map(const map<K, vector<T>>& m)
// ----------------------------------------------------------------------------
{
  vector<T> result;
  for (auto it : m)
    result.insert(result.end(), it.second.begin(), it.second.end());
  return result;
	    
}

// ----------------------------------------------------------------------------
template<typename C>
const C sort_container(C co)
// ----------------------------------------------------------------------------
{
  sort(co.begin(), co.end());
  return co;
}

// ----------------------------------------------------------------------------
template<typename T>
vector<T> remove_duplicates(vector<T>& c)
// ----------------------------------------------------------------------------
{
  auto tmp = sort_container(c);
  const auto it = unique(tmp.begin(), tmp.end());
  tmp.erase(it, tmp.end());
  return tmp;
}

// ----------------------------------------------------------------------------
CurveVec single_isocontour_merge(const vector<CurveVec>& curves,
				 const vector<LRSurfPtr>& surf_patches,
				 const double tol)
// ----------------------------------------------------------------------------
{
  map<double, CurveVec> u_map, v_map; // map curves exiting patch domains along u or v
					// parameter direction
  for (int i = 0; i != (int)surf_patches.size(); ++i) {
    prepare_curvemaps(parameter_domain(surf_patches[i]), tol, curves[i], u_map, v_map);
  }

  // identify intersection curves that do not need to be merged, and output them directly
  CurveVec result; // this will contain all isocontours (merged if necessary)
  const auto all_icurves = sort_container(expand_vec(curves));
  const auto mapped_icurves = sort_container(expand_map(map_combine(u_map, v_map)));
  set_difference(all_icurves.begin(), all_icurves.end(),
		 mapped_icurves.begin(), mapped_icurves.end(), back_inserter(result));

  // remove mapped entries corresponding with outer boundaries (no merge will take place across
  // these)
  const array<double, 4> outer_bnd = outer_boundary_pvals(surf_patches);
  CurveVec bcurves;
  auto it = u_map.find(outer_bnd[0]); if (it != u_map.end()) {add_to_vec(bcurves, it->second); u_map.erase(it);}
  it      = u_map.find(outer_bnd[1]); if (it != u_map.end()) {add_to_vec(bcurves, it->second); u_map.erase(it);}
  it      = v_map.find(outer_bnd[2]); if (it != v_map.end()) {add_to_vec(bcurves, it->second); v_map.erase(it);}
  it      = v_map.find(outer_bnd[3]); if (it != v_map.end()) {add_to_vec(bcurves, it->second); v_map.erase(it);}

  // looping through each parameter line, merging curve segments, and spitting out finished
  // curves
  CurveVec merged_segs;
  merge_segments(u_map, v_map, XFIXED, tol, bcurves, merged_segs);
  merge_segments(v_map, u_map, YFIXED, tol, bcurves, merged_segs);

  add_to_vec(result, remove_duplicates(bcurves));
  add_to_vec(result, merged_segs);
  
  return result;
}

  
// ----------------------------------------------------------------------------
vector<CurveVec> merge_isocontours(const vector<vector<CurveVec>>& curve_fragments,
				   const vector<LRSurfPtr>& surf_frags,
				   const double tol)
// ----------------------------------------------------------------------------
{
  // curve_fragments is indexed [surface patch][isovalue][set of curves]
  const int num_patches	    = (int)curve_fragments.size();  assert(num_patches > 0);
  const int num_isocontours = (int)curve_fragments[0].size();
  
  vector<CurveVec> result; // should contain one entry per isovalue

  for (int i = 0; i != num_isocontours; ++i) {

    // collect all curve fragments that belong to the set of isocurves for a
    // particular isovalue.
    const auto isocurves_to_merge =
      apply_transform(curve_fragments, function<CurveVec(const vector<CurveVec>&)> {
	  [i] (const vector<CurveVec>& v) {return v[i];}});

    // merge the fragments into complete curves, and store the resulting curves
    // in 'result'
    result.emplace_back(single_isocontour_merge(isocurves_to_merge, surf_frags, tol));
  }

  return result;
}


  
}; //end anonymous namespace
