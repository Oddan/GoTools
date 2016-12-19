#include <fstream>
#include <iostream>

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRTraceIsocontours.h"


using namespace std;
using namespace Go;

namespace {
  vector<double> contour_vals(const LRSplineSurface& lrs, int num_contours);
}// end anonymous namespace


int main()
{
  // Establishing test LR spline surface
  ifstream is("data/64_lr_1d.g2");
  ObjectHeader header;
  header.read(is);
  LRSplineSurface lrsurf(is);
  is.close();

  // Computing isocontours
  
  const bool include_3D = true;
  const bool use_sisl_marching = false;
  const auto isovals = contour_vals(lrsurf, 4);
  const vector<CurveVec> curves = LRTraceIsocontours(lrsurf,
						     isovals,
						     include_3D,
						     use_sisl_marching);

  
  cout << "Number of curves found: " << endl;
  for (size_t i = 0; i != curves.size(); ++i) 
    cout << "At height " << isovals[i] << ": " << curves[i].size() << " curves." << endl;
  
  return 0;
}



namespace {


// =============================================================================
vector<double> contour_vals(const LRSplineSurface& lrs, int num_contours)
// =============================================================================
{
  const auto mM =
    minmax_element(lrs.basisFunctionsBegin(),
		   lrs.basisFunctionsEnd(),
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


} // end anonymous namespace
