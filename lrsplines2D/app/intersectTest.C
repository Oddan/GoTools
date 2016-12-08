#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"


using namespace Go;
using namespace std;

const string filename("data/64_lr_1d.g2");
const string savefile("result.g2");

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
				     [] (std::map<LRSplineSurface::BSKey, std::unique_ptr<LRBSpline2D> >::iterator it1,
					 std::map<LRSplineSurface::BSKey, std::unique_ptr<LRBSpline2D> >::iterator it2)
				     {return it1->second->coefTimesGamma()[0] <
				             it1->second->coefTimesGamma()[0];});
  const double minval = minmax.first->second->coefTimesGamma()[0];
  const double maxval = minmax.second->second->coefTimesGamma()[0];
  
  
  lrs.to3D();
  lrs.expandToFullTensorProduct();
  auto splineSurf = shared_ptr<SplineSurface>(lrs.asSplineSurface());


  


  
  ofstream os(savefile.c_str());
  splineSurf->writeStandardHeader(os);
  splineSurf->write(os);
  os.close();
  
};
