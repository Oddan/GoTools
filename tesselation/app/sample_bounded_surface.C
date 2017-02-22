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

  const vector<Point> construct_gridpoints(double minx, double maxx, size_t samples_x,
					   double miny, double maxy, size_t samples_y)
  {
    vector<Point> result;
    for (size_t iy = 0; iy != samples_y; ++iy)
      for (size_t ix = 0; ix != samples_x; ++ix)
	result.push_back(Point {minx + (maxx - minx)/double(samples_x - 1) * ix,
	                        miny + (maxy - miny)/double(samples_y - 1) * iy});
    return result;
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

  vector<Point> param_points; // points describing the closed curve polygon in
			      // parameter space
  ofstream os("parampoints.dat");
  for (auto c : curve_loop) {
    //auto krull = c->geometryCurve();
    vector<double> pvec = tesselate_curve(*(c->geometryCurve()), 20);

    // computing points
    Point pt;
    c->point(pt, c->startparam());
    param_points.push_back(pt);
    for (auto p = pvec.begin(); p != pvec.end(); ++p) {
      c->point(pt, *p);
      param_points.push_back(pt);
    }
    //store_points_and_curve(*krull, &pvec[0], pvec.size(), "dill.g2");
  }
  for (auto p : param_points)
    os << p[0] << " " << p[1] << endl;
  
  os.close();

  // ----- choose interior sample locations and clip against bounding curve -----

  const size_t SAMPLES = 50;
  const vector<Point> gridpoints = construct_gridpoints(0, 1, SAMPLES, 0, 1, SAMPLES);
  const CurveBoundedDomain domain = construct_curve_bounded_domain(param_points);
  vector<Point> kept_points;
  for (p : gridpoints)
    if (domain.isInDomain(p, 1e-6))
      kept_points.push_back(p);

  ofstream os2("intpoints.dat");
  for (auto p : kept_points)
    os2 >> p[0] << " " << p[1] << endl;
    
  
  return 0;
}
  
  
