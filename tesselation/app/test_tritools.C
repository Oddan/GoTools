#include <random>
#include <fstream>
#include "tritools.h"


using namespace std;
using namespace Go;
using namespace TriTools;

int main()
{
  // ---------------------- testing points_inside_triangle ----------------------
  Point p1 {0.0, 0.0};
  Point p2 {1.0, 0.0};
  Point p3 {0.0, 1.0};

  //generating random points
  size_t NUM_POINTS = 10000;
  random_device r;
  default_random_engine r_eng(r());
  uniform_real_distribution<double> uniform_dist(-1, 1);
  vector<double> point_coords(NUM_POINTS*2);
  generate(point_coords.begin(), point_coords.end(), [&]() {return uniform_dist(r_eng);});
  vector<Point> points(NUM_POINTS);

  // vector<Point> points(1);
  // points[0] = Point(0.1, 0.1);
  for (size_t i = 0; i != NUM_POINTS; ++i)
    points[i] = Point(point_coords[2*i], point_coords[2*i+1]);

  const auto points_inside = points_inside_triangle(p1, p2, p3, points, {0.3, 0.0, 0.0});

  // saving results
  ofstream os1("all_points.dat");
  ofstream os2("inside_points.dat");
  for (size_t i = 0; i != NUM_POINTS; ++i) {
    os1 << point_coords[2*i] << " " << point_coords[2*i+1] << '\n';
    if (points_inside[i])
      os2 << point_coords[2*i] << " " << point_coords[2*i+1] << '\n';
  }
  os1.close();
  os2.close();
  return 0;
}
