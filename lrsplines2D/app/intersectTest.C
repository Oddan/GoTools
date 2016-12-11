#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include "sislP.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SISLconversion.h"
//#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/IntResultsSfModel.h"

using namespace Go;
using namespace std;

using SCurvePtr =  shared_ptr<SplineCurve>;
using GPoint    = tuple<double, double, int>;

const string filename("data/64_lr_1d.g2");
const string savefile("result.g2");

void simple_SurfaceModel_test(LRSplineSurface& lrs, ofstream& os);
void simple_isocontour_test(LRSplineSurface& lrs, ofstream& os, double isoval);
vector<pair<SCurvePtr, SCurvePtr>> compute_isocontour(const SplineSurface& ss, double isoval);
  
int main()
{
  ifstream is(filename.c_str());
  if (!is) {
    cout << "File not found." << endl;
    return 0;
  }
  ObjectHeader header;
  header.read(is);
  LRSplineSurface lrs(is);
  is.close();


  // determine z bounds
  const auto minmax = minmax_element(lrs.basisFunctionsBegin(), lrs.basisFunctionsEnd(),
				     [] (const pair<const LRSplineSurface::BSKey, unique_ptr<LRBSpline2D>>& p1,
					 const pair<const LRSplineSurface::BSKey, unique_ptr<LRBSpline2D>>& p2)
				     {return p1.second->Coef()[0] < p2.second->Coef()[0];});

  const double minval = minmax.first->second->Coef()[0];
  const double maxval = minmax.second->second->Coef()[0];

  cout << "Min/max: " << minval << " " << maxval << endl;

  ofstream os(savefile.c_str());
  
  // Simple test of one intersection curve, using SurfaceModel as currently
  // defined
  //simple_SurfaceModel_test(lrs, os);

  // Test function computing a isocontour for a 1D spline function
  simple_isocontour_test(lrs, os, 1600);
  
  os.close();
  
  //lrs.expandToFullTensorProduct();
  
  
  // splineSurf->writeStandardHeader(os);
  // splineSurf->write(os);
  // os.close();
  
};

// ----------------------------------------------------------------------------
void simple_isocontour_test(LRSplineSurface& lrs, ofstream& os, double isoval)
// ----------------------------------------------------------------------------
{
  // make spline function (1D)
  const auto splfun = shared_ptr<SplineSurface>(lrs.asSplineSurface());

  auto curves = compute_isocontour(*splfun, isoval);
}

// ----------------------------------------------------------------------------
// Two first members of tuple represents parameter values of intersection point.
// Last member expresses the nature of associated intersection curve.
vector<SISLIntcurve*> get_isocontour_topology(SISLSurf* s, double isoval)
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
  auto freeall = [&] () {
    if (gpar)    free(gpar);
    if (spar)    free(spar);
    if (pretop)  free(pretop);
    if (qintdat) freeIntdat(qintdat);
    if (wcurve)  freeIntcrvlist(wcurve, jcrv);
    for (int i = 0; i < jsurf; ++i) freeIntsurf(wsurf[i]);
    if (wsurf) free(wsurf);
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

  // // represent degenerated intersection curves as one point
  // sh6degen(qo1, qo1, &qintdat, epsge, &kstat);
  // if (kstat < 0) cleanup_and_throw("SISL error in sh6degen.");

  // join periodic curves
  const int kdeg = 1;
  double eimpli = 1;
  // int_join_per(&qintdat, qo1, qo1, &eimpli, kdeg, epsge, &kstat);
  // if (kstat < 0) cleanup_and_throw("SISL error in int_join_per.");
  
  // express intersections on output format (SISL)
  if (qintdat)
    hp_s1880 (qo1, qo1, kdeg, 2, 0, qintdat, &jpt, &gpar, &spar, &pretop,
	      &jcrv, &wcurve, &jsurf, &wsurf, &kstat);
  
  // @@ still unimplemented
  vector<SISLIntcurve*> result;
  for (int i = 0; i != jcrv; ++i)
    result.push_back(wcurve[i]);
  
  freeall();
  
  return result;
}

// // ----------------------------------------------------------------------------
// pair<SCurvePtr, SCurvePtr>
// trace_isoval_curve(shared_ptr<SISLSurf> s, double p1, double p2)
// // ----------------------------------------------------------------------------
// {
//   //@@UNIMPLEMENTED
//   return pair<SCurvePtr, SCurvePtr>();
// }

// ----------------------------------------------------------------------------
// First returned value is the curve in the parametric plane.  Second returned
// value is the curve as a (1D) spline function.  (This function should be of
// constant value).
vector<pair<SCurvePtr, SCurvePtr>> 
compute_isocontour(const SplineSurface& ss, double isoval)
// ----------------------------------------------------------------------------
{
  // The following asserts spefifies prerequisites for calling this function.
  assert(ss.dimension() == 1);
  assert(ss.basis_u().isKreg() && ss.basis_v().isKreg());

  SISLSurf* sislsurf = GoSurf2SISL(ss, false);
  vector<pair<SCurvePtr, SCurvePtr>> result;
  
  try {
    // determine topology (simplified version of 1851)
    auto curves = get_isocontour_topology(sislsurf, isoval);
    result = vector<pair<SCurvePtr, SCurvePtr>>(curves.size());
    // march out curves
    
    
    // transform(guidepoints.begin(), guidepoints.end(), result.begin(),
    // 	      [&] (const GPoint& p) {
    // 		return trace_isoval_curve(sislsurf, get<0>(p), get<1>(p));
    // 	      });
  } catch (exception& e) {
    freeSurf(sislsurf);
    throw e;
  }

  freeSurf(sislsurf);

  return result;
}
  
// ----------------------------------------------------------------------------
void simple_SurfaceModel_test(LRSplineSurface& lrs, ofstream& os)
// ----------------------------------------------------------------------------
{
  lrs.to3D();
  auto splineSurf = shared_ptr<SplineSurface>(lrs.asSplineSurface());
  const int id = 1;
  auto  ftSurf = shared_ptr<ftSurface>(new ftSurface(splineSurf, id));

  auto surfvec = vector<shared_ptr<ftSurface>>(1, ftSurf);
  
  SurfaceModel smod(1e-5, // approx tol
  		    1e-5, // gap between surfaces
  		    1e-3, // threshold for whether surfaces are adjacent
  		    1e-3, // kink
  		    1e-3, // bend
  		    surfvec);
  
  //const auto plane = ftPlane( {double(0),0,1}, {double(0),0,0}); // normal, point
  const auto plane = ftPlane({double(0), 0, 1}, {double(0), 0, 1600});

  auto result = smod.intersect_plane(plane);

  cout << "Curve segments: " << result->nmbCurveSegments() << endl;
  cout << "Points: " << result->nmbIntPoints() << endl;

  auto icurve = dynamic_pointer_cast<IntResultsSfModel>(result)->getIntersectionCurves();

  surfvec[0]->surface()->writeStandardHeader(os);
  surfvec[0]->surface()->write(os);

  cout << "Num segments: " << icurve.numSegments() << endl;
  for (int i = 0; i != icurve.numSegments(); ++i) {
    icurve.segment(i).spaceCurve()->writeStandardHeader(os);
    icurve.segment(i).spaceCurve()->write(os);
  }
  
}
