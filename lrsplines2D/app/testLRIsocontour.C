#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>

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
using CurvePtr  = shared_ptr<SplineCurve>;

struct IsectCurve {
  CurvePtr pcurve; // curve in parameter plane (2D)
  CurvePtr scurve; // space curve (3D)
};

using CurveVec  = vector<IsectCurve>;

namespace { // anonymous namespace
  
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

int main(int varnum, char* vararg[]) {

  // read LR-spline surface
  const auto lrs = read_surface(varnum > 1 ? vararg[1] : default_filename);
  
  // generate vectors of surface fragments
  const vector<LRSurfPtr> frags = lrs->subdivideIntoSimpler();
  
  // find isocontours on each surface fragment
  const vector<double> cvals = contour_vals(lrs, varnum > 2 ?
					    stoi(vararg[2]) :
					    default_cvals);

  vector<vector<CurveVec>> curve_fragments;
  transform(frags.begin(), frags.end(), back_inserter(curve_fragments),
	    [&cvals] (LRSurfPtr sp) { return computeIsocurves(sp, cvals); });


  // ofstream os(debug_file.c_str());

  // for (size_t i = 0; i != curve_fragments.size(); ++i) {
  //   const int cnum = accumulate(curve_fragments[i].begin(),
  // 				curve_fragments[i].end(),
  // 				0,
  // 				[] (const int& cur, const CurveVec& cv) {return cur + cv.size();});
    
  //   cout << "Patch " << i << " has " << cnum << " curves." << endl;
    
  //   frags[i]->writeStandardHeader(os);
  //   frags[i]->write(os);
  //   for (size_t j = 0; j != curve_fragments[i].size(); ++j) {
  //     for (size_t k = 0; k != curve_fragments[i][j].size(); ++k){
  // 	const IsectCurve& cv = curve_fragments[i][j][k];
  // 	cv.scurve->writeStandardHeader(os);
  // 	cv.scurve->write(os);
  //     }
  //   }
  // }

  // os.close();
  
  // merge isocontours
  const vector<CurveVec> curves = merge_isocontours(curve_fragments, frags);

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
inline double endpoint_on_parameter_line(CurvePtr c, double pval, Direction2D d)
// ----------------------------------------------------------------------------
{
  Point p1, p2;
  c->point(p1, c->startparam());
  c->point(p2, c->endparam());
  if (d==YFIXED) {
    swap(p1[0], p1[1]);
    swap(p2[0], p2[1]);
  }
  return fabs(p1[0] - pval) < fabs(p2[0] - pval) ? p1[1] : p2[1];
}

// ----------------------------------------------------------------------------
inline double endpoint_dist(CurvePtr c1, bool c1_at_start, CurvePtr c2, bool c2_at_start)
// ----------------------------------------------------------------------------
{
  Point p1, p2;
  c1->point(p1, c1_at_start ? c1->startparam() : c1->endparam());
  c2->point(p2, c2_at_start ? c2->startparam() : c2->endparam());
  return p1.dist2(p2);
}

// ----------------------------------------------------------------------------
IsectCurve join_splinecurves(IsectCurve& c1, IsectCurve& c2, Direction2D d, double pval)
// ----------------------------------------------------------------------------
{
  // ensure parameter directions are OK (the curves should be oriented so that the start of c2
  // can be joined to the end of c1).
  
  const array<double, 4> dists { endpoint_dist(c1.pcurve, true, c2.pcurve, true), // c1 start, c2 start
                                 endpoint_dist(c1.pcurve, true, c2.pcurve, false), // c1 start, c2 end
                                 endpoint_dist(c1.pcurve, false, c2.pcurve, true), // c1 end, c2 start
                                 endpoint_dist(c1.pcurve, false, c2.pcurve, false)}; // c1 end, c2 end

  const auto it = min_element(dists.begin(), dists.end());
  const double epsge = 1e-6; //@@
  assert(fabs(*it) < epsge);
  switch (it - dists.begin()) {
  case 0:
    c1.pcurve->reverseParameterDirection();
    c1.scurve->reverseParameterDirection();
    break;
  case 1:
    swap(c1, c2);
    break;
  case 2:
    // perfect, do nothing
    break;
  case 3:
    c2.pcurve->reverseParameterDirection();
    c2.scurve->reverseParameterDirection();
  };
  
  // joining curves
  IsectCurve result {CurvePtr {c1.pcurve->clone()}, CurvePtr {c1.scurve->clone()}};

  result.pcurve->appendCurve(c2.pcurve.get());
  result.scurve->appendCurve(c2.scurve.get());
  
  return result;
}

// ----------------------------------------------------------------------------
void replace_segments(const IsectCurve& old1, const IsectCurve& old2,
		      const IsectCurve& updated, map<double, CurveVec>& target)
// ----------------------------------------------------------------------------
{
  for (auto& it_map : target) 
    for (auto& it_vec : it_map.second)
      if ((it_vec.pcurve == old1.pcurve) || (it_vec.pcurve == old2.pcurve))
	it_vec = updated;
}

// ----------------------------------------------------------------------------
void merge_segments(map<double, CurveVec>& mergemap, // map whose segments should be merged
		    map<double, CurveVec>& othermap, // map whose entries must also be updated
		    Direction2D d,                   // the fixed parameter direction
		    CurveVec& finished_curves) // insert finished curves here
// ----------------------------------------------------------------------------
{
  for (auto it = mergemap.begin(); it != mergemap.end(); ++it) {
    // sorting incident curve segments so that those that should be merged will lie right next
    // to each other
    CurveVec& incident_curves = it->second;

    ofstream os("dilldall.g2");
    for (int i = 0; i != incident_curves.size(); ++i) {
      incident_curves[i].scurve->writeStandardHeader(os);
      incident_curves[i].scurve->write(os);
    }
    os.close();
    
    assert(incident_curves.size() % 2 == 0);  // we suppose that for each curve going in, there is one going out
    sort(incident_curves.begin(), incident_curves.end(),
	 [&it, d](const IsectCurve& c1, const IsectCurve& c2) {
	   return endpoint_on_parameter_line(c1.pcurve, it->first, d) < 
   	          endpoint_on_parameter_line(c2.pcurve, it->first, d);
	 });
    
    while (!incident_curves.empty()) {
      // pop the two last elements
      auto cp1 = incident_curves.back(); incident_curves.pop_back(); 
      auto cp2 = incident_curves.back(); incident_curves.pop_back();
      if (cp1.pcurve == cp2.pcurve) {
	// this is a single curve whose endpoints meet across this edge.  There are no more
	// merges to be done.  This curve is finished and can be returned.
	finished_curves.push_back(cp1);
      } else {
	auto cp_new = join_splinecurves(cp1, cp2, d, it->first);
	// replace references to cp1 and cp2 with references to cp_new throughout
	replace_segments(cp1, cp2, cp_new, mergemap);
	replace_segments(cp1, cp2, cp_new, othermap);
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
CurveVec single_isocontour_merge(const vector<CurveVec>& curves,
				 const vector<LRSurfPtr>& surf_patches)
// ----------------------------------------------------------------------------
{
  map<double, CurveVec> u_map, v_map; // map curves exiting patch domains along u or v
					// parameter direction
  ofstream os("krull.g2");
  for (auto cv : curves)
    for (auto c : cv) {
      c.scurve->writeStandardHeader(os);
      c.scurve->write(os);
    }
  os.close();
  
  for (int i = 0; i != (int)surf_patches.size(); ++i) {
    prepare_curvemaps(parameter_domain(surf_patches[i]), curves[i], u_map, v_map);
  }

  // remove mapped entries corresponding with outer boundaries (no merge will take place across
  // these)
  const array<double, 4> outer_bnd = outer_boundary_pvals(surf_patches);
  u_map.erase(outer_bnd[0]);
  u_map.erase(outer_bnd[1]);
  v_map.erase(outer_bnd[2]);
  v_map.erase(outer_bnd[3]);
  
  // looping through each parameter line, merging curve segments, and spitting out finished
  // curves
  CurveVec result;
  merge_segments(u_map, v_map, XFIXED, result);
  merge_segments(v_map, u_map, YFIXED, result);

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
  
  vector<CurveVec> result; // should contain one entry per isovalue

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
