#include <fstream>
#include <iostream>
#include <string>
#include "tesselate_curve.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"

#include "tesselate_debug.h"

using namespace std;
using namespace Go;

namespace {
  void read_away_header(ifstream& is) {
    for (int i=0, tmp=0;i != 4; is >> tmp, ++i); // read away header  
  }
  
}; // anonymous namespace 

int main(int varnum, char* vararg[]) {

  // ------------------------------ loading object ------------------------------
  const string filename = (varnum > 1) ? vararg[1] : "curve_and_surface.g2";

  ifstream is(filename.c_str());

  read_away_header(is);
  shared_ptr<ParamSurface> ss {new SplineSurface};
  ss->read(is);

  // read four curves

  GoTools::init();
  vector<shared_ptr<CurveOnSurface>> curve_loop;
  for (size_t i = 0; i != 4; ++i) {
    read_away_header(is);
    shared_ptr<CurveOnSurface> c {new CurveOnSurface};
    c->read(is);
    c->setUnderlyingSurface(ss);
    curve_loop.emplace_back(c);
  }

  // ----------------------------- tesselate curves -----------------------------

  ofstream os("krull2.g2");
  for (auto c : curve_loop) {
    auto krull = c->geometryCurve();
    krull->writeStandardHeader(os);
    krull->write(os);
    vector<double> pvec = tesselate_curve(*(c->geometryCurve()), 20);

    store_points_and_curve(*krull, &pvec[0], pvec.size(), "dill.g2");
  }
  os.close();
  
  return 0;
}
  
  
