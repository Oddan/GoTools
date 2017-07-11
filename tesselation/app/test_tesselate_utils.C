#include <fstream>
#include <iostream>
#include <vector>
#include <iterator>
#include "tesselate_utils.h"
#include "GoTools/utils/Point.h"
#include "interpoint_distances.h"

using namespace std;
using namespace TesselateUtils;

namespace Test { 
  void test_poly_area();
  void test_bounding_box();
  void test_generate_grid();
  void test_inpolygon();
  void test_interpoint_distances();
};

// ============================================================================
int main() {

  // cout << "Testing polygon area: " << endl;
  // Test::test_poly_area();

  // cout << "Testing bounding box: " << endl;
  // Test::test_bounding_box();

  // cout << "Testing generate_grid: " << endl;
  // Test::test_generate_grid();

  // cout << "Testing inpolygon: " << endl;
  // Test::test_inpolygon();
  
  cout << "Testing interpoint distance computation: " << endl;
  Test::test_interpoint_distances();

  return 0;
};

// ============================================================================

namespace Test {
  
// ----------------------------------------------------------------------------
void test_poly_area()
// ----------------------------------------------------------------------------
{
  vector<Go::Point> poly { {0.0, 0.0}, {1.0, 0.0}, {0.5, 0.5}, {0.0, 1.0}};

  cout << "Area is: " << polygon_area(&poly[0], (int)poly.size()) << endl;
}

// ----------------------------------------------------------------------------
void test_bounding_box()
// ----------------------------------------------------------------------------
{
  vector<Go::Point> poly { {0.0, 0.0}, {1.0, 0.0}, {0.5, 0.5}, {0.0, 1.0}};
  
  const auto box = bounding_box(&poly[0], (int)poly.size());
  cout << "Bounding box is: \n";
  cout << "xmin: " << box[0] << ", xmax: " << box[1] << '\n';
  cout << "ymin: " << box[2] << ", ymax: " << box[3] << endl;
}

// ----------------------------------------------------------------------------  
void test_generate_grid()
// ----------------------------------------------------------------------------  
{
  const Go::Point p1 {0.0, 0.0};
  const Go::Point p2 {1.0, 2.0};
  const int nx = 4;
  const int ny = 6;

  const auto result = generate_grid_2D(p1, p2, nx, ny);

  cout << "Gridpoints are: \n";
  copy(result.begin(), result.end(), std::ostream_iterator<Go::Point>(cout, "\n"));
  cout << endl;
}

// ----------------------------------------------------------------------------
void test_inpolygon()
// ----------------------------------------------------------------------------  
{
  const vector<Go::Point> poly = { {0.0, 0.0}, {1.0, 0.0}, 
				   {1.0, 2.0}, {0.5, 1.0}, {0.0, 1.0}};
  const auto bbox = bounding_box(&poly[0], (int)poly.size());
  const int nx = 30;
  const int ny = 30;
  const double tol = 1e-6;
  const auto candidates = generate_grid_2D(Go::Point {bbox[0], bbox[2]}, 
					   Go::Point {bbox[1], bbox[3]}, nx, ny);
  const auto result = inpolygon(&candidates[0], (unsigned int)candidates.size(),
				&poly[0], (unsigned int)poly.size(), tol);

  cout << "Polygon is: \n";
  copy(poly.begin(), poly.end(), std::ostream_iterator<Go::Point>(cout, "\n"));

  cout << "Candidate points were: \n";
  copy(candidates.begin(), candidates.end(), std::ostream_iterator<Go::Point>(cout, "\n"));

  cout << "Interior points are: \n" << endl;
  copy(result.begin(), result.end(), std::ostream_iterator<Go::Point>(cout, "\n"));

}

// ----------------------------------------------------------------------------  
void test_interpoint_distances()
// ----------------------------------------------------------------------------
{
  // Define a random set of test points
  const int N = 100;
  const double R = 0.05;
  vector<Point2D> points(N);
  generate(points.begin(), points.end(), 
	   [](){return Point2D {random_uniform(0, 1), random_uniform(0, 1)};});

  const auto result = interpoint_distances(&points[0], N, R, false);

  ofstream os("distances.mat");
  for (auto e : result)
    os << e.p1_ix << " " << e.p2_ix << " " << e.dist << '\n';
  os.close();

  ofstream os2("points.mat");
  for (auto p : points)
    os2 << p[0] << " " << p[1] << '\n';
  os2.close();

  cout << "Number of relations found: " << result.size() << endl;
    
}

};
