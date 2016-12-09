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
  simple_SurfaceModel_test(lrs, os);

  // Test function computing a isocontour for a 1D spline function
  //simple_isocontour_test(lrs, os, 1600);
  
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
vector<GPoint> get_isocontour_topology(shared_ptr<SISLSurf> s, double isoval)
// ----------------------------------------------------------------------------
{
  // This function does the equivalent of SISL sh1851, but assumes that the
  //provided surface is 1D (bivariate spline function), noncyclic and k-regular.

  
  // Finding intersections, using SISL function sh1761.
  SISLIntdat* qintdat = SISL_NULL; // intersection result
  auto qp = shared_ptr<SISLPoint>(newPoint(&isoval, 1, 1));
  auto qo1 = shared_ptr<SISLObject>(newObject(SISLSURFACE));
  qo1->s1 = s.get();
  qo1->o1 = qo1.get();
  auto qo2 = shared_ptr<SISLObject>(newObject(SISLPOINT));
  qo2->p1 = qp.get();
  const double epsge = 1e-6;
  int kstat = 0;
  sh1761(qo1.get(), qo2.get(), epsge, &qintdat, &kstat);
  if (kstat < 0) {
    if (qintdat) freeIntdat(qintdat);
    throw runtime_error("SISL error in sh1761");
  }
  
  // @@ still unimplemented
  if (qintdat) freeIntdat(qintdat);
  return vector<GPoint>();
}

// ----------------------------------------------------------------------------
pair<SCurvePtr, SCurvePtr>
trace_isoval_curve(shared_ptr<SISLSurf> s, double p1, double p2)
// ----------------------------------------------------------------------------
{
  //@@UNIMPLEMENTED
  return pair<SCurvePtr, SCurvePtr>();
}

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
  
  auto sislsurf = shared_ptr<SISLSurf>(GoSurf2SISL(ss, false));
                                       
  // determine topology (simplified version of 1851)
  auto guidepoints = get_isocontour_topology(sislsurf, isoval);
  
  // march out curves
  vector<pair<SCurvePtr, SCurvePtr>> result(guidepoints.size());

  transform(guidepoints.begin(), guidepoints.end(), result.begin(),
	    [&] (const GPoint& p) {
	      return trace_isoval_curve(sislsurf, get<0>(p), get<1>(p));
	    });

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
