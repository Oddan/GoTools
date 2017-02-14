#include <vector>
#include <iostream>
#include <fstream>
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"

using namespace std;
using namespace Go;

int main(int varnum, char* vararg[]) {

  const int cpoints = 5;
  
  vector<double> coefs;
  for (size_t i = 0; i != cpoints; ++i) {
    for (size_t j = 0; j != cpoints; ++j) {
      coefs.push_back(double(i));
      coefs.push_back(double(j));
      coefs.push_back(0);
      if ((i==3) & (j == 2))
	coefs.back() = 3;
      else if ((i==1) & (j==1))
	coefs.back() = -1;
    }
  }
  // ---------------------- Make underlying spline surface ----------------------
  shared_ptr<SplineSurface>
    ss(new SplineSurface(cpoints, cpoints, 
			 4, 4,  // order in either parameter direction
			 &(vector<double> {0, 0, 0, 0, 0.5, 1, 1, 1, 1})[0],
			 &(vector<double> {0, 0, 0, 0, 0.2, 1, 1, 1, 1})[0],
			 &coefs[0],
			 3, // dimension
			 false)); // rationals


  // --------------------------- make parameter-curves ---------------------------

  // straight line segment 1
  
  shared_ptr<SplineCurve>
    sc1(new SplineCurve(2, 2, // # of coefs, and order
			&(vector<double> {0, 0, 1, 1})[0],   // knots
			&(vector<double> {4/5.0, 3/5.0, 4/5.0, 4.5/5.0})[0], // coefs
			2)); // dimension (parameter domain has dim 2)

  shared_ptr<CurveOnSurface> cos1(new CurveOnSurface(ss, sc1, true));

  shared_ptr<SplineCurve>
    sc2(new SplineCurve(7, 2,
			&(vector<double> {0, 0, 1/6.0, 2/6.0, 3/6.0, 4/6.0, 5/6.0, 6/6.0, 6/6.0})[0],
			&(vector<double> {4/5.0, 4.5/5.0,
			      2/5.0,   4.5/5.0,
			      0.5/5.0, 4.5/5.0,
			      0.5/5.0, 2.5/5.0,
			      0.5/5.0, 0.5/5.0,
			      2/5.0,   0.5/5.0,
			      4/5.0,   0.5/5.0})[0],
			2));
			
  shared_ptr<CurveOnSurface> cos2(new CurveOnSurface(ss, sc2,true));
  
  shared_ptr<SplineCurve> tmp(cos1->geometryCurve());
  shared_ptr<SplineCurve> tmp2(cos2->geometryCurve());
  
  
  // ------------------------------- Write result -------------------------------
  
  ofstream os("bsurf.g2");
  ss->writeStandardHeader(os);
  ss->write(os);
  cos1->writeStandardHeader(os);
  cos1->write(os);
  cos2->writeStandardHeader(os);
  cos2->write(os);
  
  return 0;
};
