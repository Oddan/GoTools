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


  ofstream os(debug_file.c_str());

  for (size_t i = 0; i != curve_fragments.size(); ++i) {
    const int cnum = accumulate(curve_fragments[i].begin(),
				curve_fragments[i].end(),
				0,
				[] (const int& cur, const CurveVec& cv) {return cur + cv.size();});
    
    cout << "Patch " << i << " has " << cnum << " curves." << endl;
    
    frags[i]->writeStandardHeader(os);
    frags[i]->write(os);
    for (size_t j = 0; j != curve_fragments[i].size(); ++j) {
      for (size_t k = 0; k != curve_fragments[i][j].size(); ++k){
	const IsectCurve& cv = curve_fragments[i][j][k];
	cv.scurve->writeStandardHeader(os);
	cv.scurve->write(os);
      }
    }
  }

  os.close();
  
  
  // merge isocontours
	      
  

  return 0;
}

namespace {
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
