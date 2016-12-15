#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <chrono>

#include "sislP.h"
#include "GoTools/geometry/SISLconversion.h"

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"

using namespace Go;
using namespace std;

using LRSurfPtr = shared_ptr<LRSplineSurface>;
using CurvePtr  = shared_ptr<const SplineCurve>;

struct IsectCurve {
  CurvePtr pcurve; // curve in parameter plane (2D)
  CurvePtr scurve; // space curve (3D)

  // IsectCurve(const IsectCurve& other)
  //   : pcurve(other.pcurve), scurve(other.scurve) {}
  bool operator<(const IsectCurve& rhs) const {
    return (pcurve.get() < rhs.pcurve.get());
  }

  bool operator==(const IsectCurve& rhs) const {
    return pcurve.get() == rhs.pcurve.get();
  }
};

using CurveVec  = vector<IsectCurve>;

// ----------------------------------------------------------------------------
namespace { // anonymous namespace
// ----------------------------------------------------------------------------  
  const string default_filename("data/64_lr_1d.g2");
  const string debug_file("result.g2");
  const int default_cvals = 10;

  LRSurfPtr read_surface(string fname);
  vector<double> contour_vals(LRSurfPtr lrs, int num_contours);
  vector<CurveVec> computeIsocurves(const LRSurfPtr lrs,
				    const vector<double>& isovals);
  vector<CurveVec> merge_isocontours(const vector<vector<CurveVec>>& curve_fragments,
				     const vector<LRSurfPtr>& surf_frags);
}; // end anonymous namespace


// ============================================================================
int main(int varnum, char* vararg[]) {
// ============================================================================
  // read LR-spline surface
  const auto lrs = read_surface(varnum > 1 ? vararg[1] : default_filename);
  
  // generate vectors of surface fragments
  auto t1 = chrono::high_resolution_clock::now();
  const vector<LRSurfPtr> frags = lrs->subdivideIntoSimpler();
  auto t2 = chrono::high_resolution_clock::now();
  cout << "** Subdivision took: " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " milliseconds. " << endl;
 
  // find isocontours on each surface fragment

  const vector<double> cvals = contour_vals(lrs, varnum > 2 ?
					    stoi(vararg[2]) :
					    default_cvals);

  t1 = chrono::high_resolution_clock::now();  
  vector<vector<CurveVec>> curve_fragments;
  transform(frags.begin(), frags.end(), back_inserter(curve_fragments),
	    [&cvals] (LRSurfPtr sp) { return computeIsocurves(sp, cvals); });
  t2 = chrono::high_resolution_clock::now();
  cout << "** Finding isocontours took: " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " milliseconds. " << endl;
  
  // merge isocontours
  t1 = chrono::high_resolution_clock::now();  
  const vector<CurveVec> curves = merge_isocontours(curve_fragments, frags);
  t2 = chrono::high_resolution_clock::now();
  cout << "** Merging segments took: " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " milliseconds. " << endl;  

  ofstream os(debug_file.c_str());
  for (auto c_level : curves)
    for (auto c : c_level) {
      c.scurve->writeStandardHeader(os);
      c.scurve->write(os);
    }
  os.close();

  return 0;
}

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
int find_exit_edge(const IsectCurve& c, const array<double, 4>& domain, bool start)
// ----------------------------------------------------------------------------
{
  const double epsge = 1e-6;
  const double t = start ? c.pcurve->startparam() : c.pcurve->endparam();
  Point uv;  // represent u and v parameters
  c.pcurve->point(uv, t);  // evaluate (u, v) parameter pair at curve start or end

  const vector<double> dists = {fabs(domain[0] - uv[0]),
				fabs(domain[1] - uv[0]),
				fabs(domain[2] - uv[1]),
				fabs(domain[3] - uv[1])};
  auto min_it = min_element(dists.begin(), dists.end());
  if (*min_it > epsge)
    return -1; // curve doesn't end on edge, but inside domain (either closed, or ends at
	       // singularity)

  return int(min_it - dists.begin()); // index of boundary edge where curve starts/ends
		 
}

// ----------------------------------------------------------------------------
void map_curve(const IsectCurve& c, const array<double, 4>& domain, bool at_start,
	       map<double, CurveVec>& u_map, map<double, CurveVec>& v_map)
// ----------------------------------------------------------------------------
{
  const int edge = find_exit_edge(c, domain, at_start);

  const double epsge = 1e-6;
  if (edge >=0) {
    Point p1; c.pcurve->point(p1, at_start ? c.pcurve->startparam() : c.pcurve->endparam());
    assert( (fabs(domain[edge] - p1[0]) < epsge) | (fabs(domain[edge] - p1[1]) < epsge) );
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

// void check_map_integrity(const map<double, CurveVec>& m, Direction2D d)
// {
//   for (auto it : m) {
//     const double pval = it.first;
//     for (auto v : it.second) {
//       Point p1, p2;
//       v.pcurve->point(p1, v.pcurve->startparam());
//       v.pcurve->point(p2, v.pcurve->endparam());
//       const int ix = (d==XFIXED) ? 0 : 1;
//       const double d1 = fabs(p1[ix]-pval);
//       const double d2 = fabs(p2[ix]-pval);

//       assert(min(d1, d2) < 1e-6);
//     }
//   }
// }
      
  
// ----------------------------------------------------------------------------
void prepare_curvemaps(const array<double, 4>& domain, const CurveVec& cvec,
		       map<double, CurveVec>& u_map, map<double, CurveVec>& v_map)
// ----------------------------------------------------------------------------
{
  for_each(cvec.begin(), cvec.end(), [&] (const IsectCurve& c) {
      map_curve(c, domain, true, u_map, v_map);  // map curve according to start point
      map_curve(c, domain, false, u_map, v_map); // map curve according to end point
  });
}

// ----------------------------------------------------------------------------
IsectCurve join_isectcurves(const IsectCurve& c1, const IsectCurve& c2,
			    bool c1_at_start, bool c2_at_start)
// ----------------------------------------------------------------------------
{
  SplineCurve pcurve1 = SplineCurve(*c1.pcurve);
  SplineCurve scurve1 = SplineCurve(*c1.scurve);
  SplineCurve pcurve2 = SplineCurve(*c2.pcurve);
  SplineCurve scurve2 = SplineCurve(*c2.scurve);

  if (c1_at_start) {
    pcurve1.reverseParameterDirection();
    scurve1.reverseParameterDirection();
  }

  if (!c2_at_start) {
    pcurve2.reverseParameterDirection();
    scurve2.reverseParameterDirection();
  }

  pcurve1.appendCurve(&pcurve2);
  scurve1.appendCurve(&scurve2);
  
  return IsectCurve { CurvePtr(pcurve1.clone()), CurvePtr(scurve1.clone())};
}

// ----------------------------------------------------------------------------
void replace_segments(const IsectCurve& old1, const IsectCurve& old2,
		      const IsectCurve& updated, map<double, CurveVec>& target,
		      Direction2D d) //@@ d is here only for debugging.  Can be removed
// ----------------------------------------------------------------------------
{
  for (auto& it_map : target) 
    for (auto& it_vec : it_map.second) 
      if ((it_vec.pcurve == old1.pcurve) || (it_vec.pcurve == old2.pcurve)) 
	it_vec = updated;
}


// ----------------------------------------------------------------------------
pair<double, double> identify_truncated_endpoints(const IsectCurve& c, Direction2D d, double pval)
// ----------------------------------------------------------------------------
{
  const double epsge = 1e-6;
  Point startpoint, endpoint;
  c.pcurve->point(startpoint, c.pcurve->startparam());
  c.pcurve->point(endpoint  , c.pcurve->endparam());
  if (d==YFIXED) {
    swap(startpoint[0], startpoint[1]);
    swap(endpoint[0], endpoint[1]);
  }
  auto NaN = numeric_limits<double>::quiet_NaN();
  return pair<double, double> {
    (fabs(startpoint[0] - pval) < epsge ? startpoint[1] : NaN),
    (fabs(endpoint[0] - pval) < epsge ? endpoint[1] : NaN)
  };
}


// ----------------------------------------------------------------------------
void merge_segments(map<double, CurveVec>& mergemap, // map whose segments should be merged
		    map<double, CurveVec>& othermap, // map whose entries must also be updated
		    Direction2D d,                   // the fixed parameter direction
		    CurveVec& bcurves,            
		    CurveVec& finished_curves) // insert newly merged finished curves here
// ----------------------------------------------------------------------------
{
  const double epsge = 1e-6;
  struct EndPoint {double pval; IsectCurve icurve;bool at_start;};

  while (!mergemap.empty()) {
    auto it = *mergemap.begin();   mergemap.erase(mergemap.begin());
    // sorting incident curve segments so that those that should be merged will lie right next
    // to each other
    vector<CurvePtr> encountered;
    vector<EndPoint> tp_vec;
    for (auto& ic : it.second) { // loop over intersection curve segments cut by this parameter line
      const auto ends = identify_truncated_endpoints(ic, d, it.first);
      if (!isnan(ends.first + ends.second)) {
	// both endpoints of this curve are truncated by the parameter line.
	if (find(encountered.begin(), encountered.end(), ic.pcurve) != encountered.end())
	  continue;
	encountered.push_back(ic.pcurve); // ensure we will only register the curve once
      }
      if (!isnan(ends.first))  tp_vec.push_back({ends.first, ic, true});
      if (!isnan(ends.second)) tp_vec.push_back({ends.second, ic, false});
    }

    assert(tp_vec.size() % 2 == 0);  // we suppose that for each curve going in, there is one going out

    sort(tp_vec.begin(), tp_vec.end(), [](const EndPoint& t1, const EndPoint& t2) {return t1.pval < t2.pval;});

    for (size_t i = 0; i != tp_vec.size(); i += 2) {
      const auto entry1 = tp_vec[i];
      const auto entry2 = tp_vec[i+1];

      assert(fabs(entry1.pval - entry2.pval) < epsge);
      if (entry1.icurve.pcurve == entry2.icurve.pcurve) {
	// this is a single curve whose endpoints meet across this edge.  There are no more
	// merges to be done.  The curve is finished, and can be returned.
	finished_curves.push_back(entry1.icurve);
      } else {
	auto new_curve = join_isectcurves(entry1.icurve, entry2.icurve, entry1.at_start, entry2.at_start);

	// replace references to the old curves with references to new_curve throughout
	replace_segments(entry1.icurve, entry2.icurve, new_curve, mergemap, d);
	replace_segments(entry1.icurve, entry2.icurve, new_curve, othermap, flip(d));

	for (size_t j = i+2; j < tp_vec.size(); ++j) 
	  if ((tp_vec[j].icurve.pcurve == entry1.icurve.pcurve) |
	      (tp_vec[j].icurve.pcurve == entry2.icurve.pcurve)) {
	    tp_vec[j].icurve = new_curve;
	    const auto ends = identify_truncated_endpoints(new_curve, d, it.first);
	    tp_vec[j].at_start = (!isnan(ends.first));
	    // if ends.second is also truncated by this line (i.e. !isnan(ends.second)), we
	    // have just closed a loop.  It will be picked up in one of the upcoming
	    // iterations.
	  }

	// updating boundary curve pointers if necessary
	transform(bcurves.begin(), bcurves.end(), bcurves.begin(), [&](const IsectCurve& c) {
	    return ((c.pcurve == entry1.icurve.pcurve) | (c.pcurve == entry2.icurve.pcurve)) ? new_curve : c;});
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
				 const vector<LRSurfPtr>& surf_patches)
// ----------------------------------------------------------------------------
{
  map<double, CurveVec> u_map, v_map; // map curves exiting patch domains along u or v
					// parameter direction
  for (int i = 0; i != (int)surf_patches.size(); ++i) {
    prepare_curvemaps(parameter_domain(surf_patches[i]), curves[i], u_map, v_map);
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
  merge_segments(u_map, v_map, XFIXED, bcurves, merged_segs);
  merge_segments(v_map, u_map, YFIXED, bcurves, merged_segs);

  add_to_vec(result, remove_duplicates(bcurves));
  add_to_vec(result, merged_segs);
  
  return result;
}
  
// =============================================================================
vector<CurveVec> merge_isocontours(const vector<vector<CurveVec>>& curve_fragments,
				   const vector<LRSurfPtr>& surf_frags)
// =============================================================================
{
  // curve_fragments is indexed [surface patch][isovalue][set of curves]
  const int num_patches	    = (int)curve_fragments.size();  assert(num_patches > 0);
  const int num_isocontours = (int)curve_fragments[0].size();
  
  vector<CurveVec> result; // should contain one entry1 per isovalue

  for (int i = 0; i != num_isocontours; ++i) {
    vector<CurveVec> isocurves_to_merge; // isocurve(s) for each surface patch
    transform(curve_fragments.begin(),
	      curve_fragments.end(),
	      back_inserter(isocurves_to_merge),
	      [i] (const vector<CurveVec>& v) {return v[i];});

    result.emplace_back(single_isocontour_merge(isocurves_to_merge, surf_frags));
  }

  return result;
}

// =============================================================================
LRSurfPtr read_surface(string fname)
// =============================================================================
{
  ifstream is(fname.c_str());

  if (!is) throw runtime_error("could not open input file.");

  ObjectHeader header;  header.read(is);  // removing header info

  return LRSurfPtr(new LRSplineSurface(is));
}

// =============================================================================
vector<double> contour_vals(LRSurfPtr lrs, int num_contours)
// =============================================================================
{
  const auto mM =
    minmax_element(lrs->basisFunctionsBegin(),
		   lrs->basisFunctionsEnd(),
		   [] (const pair<const LRSplineSurface::BSKey, unique_ptr<LRBSpline2D>>& p1,
		       const pair<const LRSplineSurface::BSKey, unique_ptr<LRBSpline2D>>& p2)
		   {return p1.second->Coef()[0] < p2.second->Coef()[0];});

  const double minval = mM.first->second->Coef()[0];
  const double maxval = mM.second->second->Coef()[0];
  
  vector<double> result(num_contours, 0);
  for (int i = 0; i != num_contours; ++i)
    result[i] = minval + ((maxval - minval)/(num_contours-1)) * i;

  return result;
}

// ----------------------------------------------------------------------------
// Two first members of tuple represents parameter values of intersection point.
// Last member expresses the nature of associated intersection curve.
pair<SISLIntcurve**, int> get_isocontour_topology(SISLSurf* s, double isoval)
// ----------------------------------------------------------------------------
{
  // This function does the equivalent of SISL sh1851, but assumes that the
  //provided surface is 1D (bivariate spline function), noncyclic and k-regular.

  
  // Finding intersections, using SISL function sh1761.
  int jpt = 0; // number of single intersection point
  int jcrv = 0; // number of intersection curves
  int jsurf = 0; 
  double* gpar = SISL_NULL; // parameter values of the single intersection points
  double* spar = SISL_NULL; // dummy array
  int *pretop=SISL_NULL;
  SISLIntcurve** wcurve; // array containing descriptions of the intersection curves
  SISLIntsurf** wsurf; 

  auto qp  = newPoint(&isoval, 1, 1);
  auto qo1 = newObject(SISLSURFACE);
  qo1->s1 = s;
  qo1->o1 = qo1;
  auto qo2 = newObject(SISLPOINT);
  qo2->p1 = qp;
  
  SISLIntdat* qintdat = SISL_NULL; // intersection result
  auto freeall = [&] (bool skip_wcurves = false) {
    if (gpar)                         free(gpar);
    if (spar)                         free(spar);
    if (pretop)                       free(pretop);
    if (qintdat)                      freeIntdat(qintdat);
    if ((bool)wcurve & (jcrv > 0) & !skip_wcurves) freeIntcrvlist(wcurve, jcrv);
    for (int i = 0; i < jsurf; ++i)   freeIntsurf(wsurf[i]);
    if ((bool)wsurf & (jsurf > 0))          free(wsurf);
  };
  
  auto cleanup_and_throw = [&freeall] (string s) {
    freeall();
    throw runtime_error(s.c_str());
  };
  
  const double epsge = 1e-6;
  int kstat = 0;
  
  // find intersections
  sh1761(qo1, qo2, epsge, &qintdat, &kstat);
  if (kstat < 0) cleanup_and_throw("SISL error in sh1761.");

  const int kdeg = 1;
    
  // express intersections on output format (SISL)
  if (qintdat)
    hp_s1880 (qo1, qo1, kdeg, 2, 0, qintdat, &jpt, &gpar, &spar, &pretop,
	      &jcrv, &wcurve, &jsurf, &wsurf, &kstat);
  
  freeall(true);
  
  return pair<SISLIntcurve**, int>(wcurve, jcrv);
}

// ----------------------------------------------------------------------------
IsectCurve trace_isoval_curve(SISLSurf* s, SISLIntcurve* ic, double isoval)
// ----------------------------------------------------------------------------
{
  double pnt[] = {0, 0, isoval};
  double nrm[] = {0, 0, 1};
  const double epsge = 1e-6;
  const double epsco = 1e-15; // Not used
  const double maxstep = 0.0;
  const int dim = 3;
  const int makecurv = 2; // make both geometric and parametric curves
  int stat;
  //@@UNIMPLEMENTED
  s1314(s, pnt, nrm, dim, epsco, epsge, maxstep, ic, makecurv, 0, &stat);

  SISLCurve* sc = ic->pgeom;
  SISLCurve* sp = ic->ppar1;
  if (sc == 0) {
    MESSAGE("s1314 returned code: " << stat << ", returning.");
    return IsectCurve {};
  }
  
  assert(sc->rcoef == 0); // otherwise, slight modification of the below is necessary
  const bool rational = false;
  return IsectCurve
  {
    CurvePtr(new SplineCurve(sp->in, sp->ik, sp->et, sp->ecoef, 2, rational)), // pcurve
    CurvePtr(new SplineCurve(sc->in, sc->ik, sc->et, sc->ecoef, 3, rational)) // scurve
  };
  
}

  
// =============================================================================
vector<CurveVec> computeIsocurves(const LRSurfPtr lrs, const vector<double>& isovals)
// =============================================================================
{
  auto lrs_working_copy = shared_ptr<LRSplineSurface>(lrs->clone());

  // make 1D spline function 
  const auto ss   = shared_ptr<SplineSurface>(lrs_working_copy->asSplineSurface());
  assert(ss->basis_u().isKreg() && ss->basis_v().isKreg());
  SISLSurf* sislsurf1D = GoSurf2SISL(*ss, false);

  // make 3D spline function
  lrs_working_copy->to3D();
  const auto ss3D = shared_ptr<SplineSurface>(lrs_working_copy->asSplineSurface());
  SISLSurf* sislsurf3D = GoSurf2SISL(*ss3D, false);
  
  vector<CurveVec> result;
  pair<SISLIntcurve**, int> curves;

  try {
    for (auto ival : isovals) {
      CurveVec isocurves;
      // determine topology (simplified version of 1851)
      curves = get_isocontour_topology(sislsurf1D, ival);
      
      // march out curves
      transform(curves.first, curves.first + curves.second, back_inserter(isocurves),
		[&] (SISLIntcurve* ic) {
		  return trace_isoval_curve(sislsurf3D, ic, ival);
		});
      result.emplace_back(isocurves);
    }
    if (curves.second > 0) freeIntcrvlist(curves.first, curves.second);
  } catch (exception& e) {
    if (curves.second > 0) freeIntcrvlist(curves.first, curves.second);
    freeSurf(sislsurf1D);
    freeSurf(sislsurf3D);
    throw e;
  }

  freeSurf(sislsurf3D);
  freeSurf(sislsurf1D);

  return result;
};

}; // end anonymous namespace
