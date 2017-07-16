#include <fstream>
#include <iostream>
#include <vector>
#include <iterator>
#include "tesselate_utils.h"
#include "GoTools/utils/Point.h"
#include "interpoint_distances.h"
#include "polyhedral_energies.h"
#include "triangulate_domain.h"
#include "SimplePolyhedronTesselation.h"

using namespace std;
using namespace TesselateUtils;


namespace Test { 
  void test_poly_p();
  void test_bounding_box();
  void test_generate_grid();
  void test_inpolygon();
  void test_interpoint_distances();
  void test_energy();
  void test_circumscribe_triangle();
  void test_segment_intersection();
  void test_triangulation();
  void test_volume_tesselation();
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
  
  // cout << "Testing interpoint distance computation: " << endl;
  // Test::test_interpoint_distances();

  // cout << "Testing energy computation: " << endl;
  // Test::test_energy();

  // cout << "testing circumscribe triangle: " << endl;
  // Test::test_circumscribe_triangle();

  // cout << "testing segment intersection: " << endl;
  // Test::test_segment_intersection();
  
  // cout << "testing triangulation generation: " << endl;
  // Test::test_triangulation();

  cout << "testing volume tesselation: " << endl;
  Test::test_volume_tesselation();
  
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
  const int N = 20000;
  const double R = 0.1;
  vector<Point2D> points(N);
  generate(points.begin(), points.end(), 
	   [](){return Point2D {random_uniform(0, 1), random_uniform(0, 1)};});

  const auto result = interpoint_distances(&points[0], N, R, false);

  // ofstream os("distances.mat");
  // for (auto e : result)
  //   os << e.p1_ix << " " << e.p2_ix << " " << e.dist << '\n';
  // os.close();

  // ofstream os2("points.mat");
  // for (auto p : points)
  //   os2 << p[0] << " " << p[1] << '\n';
  // os2.close();

  cout << "Number of relations found: " << result.size() << endl;
    
}

// ----------------------------------------------------------------------------
void test_energy()
// ----------------------------------------------------------------------------
{
  const vector<Point2D> bpoints { {0,0}, {1,0}, {1,1}, {0,1} };

  const vector<Point2D> points { {0.25, 0.5}, {0.75, 0.5}};
  const double R = 1;

  const auto result = polygon_energy(&bpoints[0], (uint)bpoints.size(),
				     &points[0],  (uint)points.size(),
				     R);

  cout << "Energy is: " << result.val << '\n';

  cout << "Energy derivatives are: \n";
  for (auto d : result.der) {
    cout << d[0] << " " << d[1] << '\n';
  }
    
}
// ----------------------------------------------------------------------------
void test_circumscribe_triangle()
// ----------------------------------------------------------------------------
{
  const Point2D p1 {0, 0};
  const Point2D p2 {1, 0};
  const Point2D p3 {1, 1};

  Point2D center;
  double radius2;
  circumscribe_triangle(p1, p2, p3, center, radius2);

  cout << "Points are: \n";
  cout << "(" << p1[0] << ", " << p1[1] << "), ";
  cout << "(" << p2[0] << ", " << p2[1] << "), ";
  cout << "(" << p3[0] << ", " << p3[1] << ")\n";
  cout << "Center is: \n";
  cout << "(" << center[0] << ", " << center[1] << ")\n";
  cout << "Radius is: " << sqrt(radius2) << endl << endl;

  cout << "Computed distances to each point are: \n";
  cout << "For p1 -- " << dist(p1, center) << endl;
  cout << "For p2 -- " << dist(p2, center) << endl;
  cout << "For p3 -- " << dist(p3, center) << endl;
  cout << endl;
}

// ----------------------------------------------------------------------------
void test_segment_intersection()
// ----------------------------------------------------------------------------
{
  // obvious case
  const Point2D p1a {1, 1}, p1b {4, 3};
  const Point2D p2a {3, 1}, p2b {2, 3};

  // endpoints-meet case
  const Point2D q1a {1, 1}, q1b {1, 4};
  const Point2D q2a {1, 4}, q2b {2, 3};

  // endpoint-touch-segment case
  const Point2D r1a {1, 1}, r1b {3, 1};
  const Point2D r2a {2, 1}, r2b {2, 2};

  const double tol = 1e-6;
  cout << "Obvious case, positive tol: " << segments_intersect_2D(p1a, p1b, p2a, p2b,  tol) << endl;
  cout << "Obvious case, negative tol: " << segments_intersect_2D(p1a, p1b, p2a, p2b, -tol) << endl;
  cout << endl;
  cout << "Endpoint-meet, positive tol: " << segments_intersect_2D(q1a, q1b, q2a, q2b,  tol) << endl;
  cout << "Endpoint-meet, negative tol: " << segments_intersect_2D(q1a, q1b, q2a, q2b, -tol) << endl;
  cout << endl;
  cout << "Point-touch-line, pos. tol: " << segments_intersect_2D(r1a, r1b, r2a, r2b,  tol) << endl;
  cout << "Point-touch-line, pos. tol: " << segments_intersect_2D(r1a, r1b, r2a, r2b, -tol) << endl;
  cout << endl;
  
  
}

// ----------------------------------------------------------------------------
void test_triangulation()
// ----------------------------------------------------------------------------
{
  vector<Point2D> bpoints { {0, 0}, {1, 0}, {1, 1}, {0,1}};
  vector<Point2D> ipoints { {0.25, 0.25}, {0.75, 0.25}, {0.5, 0.5}, {0.25, 0.75}, {0.75, 0.75}};

  vector<Point2D> all_points(bpoints);
  all_points.insert(all_points.end(), ipoints.begin(), ipoints.end());

  const double vdist = 0.9;
  vector<Triangle> tris;
  // vector<Triangle> tris = triangulate_domain(&all_points[0], (uint)bpoints.size(),
  // 					     (uint)all_points.size(), vdist);

  // cout << "Points: " << endl;
  // for (const auto p : all_points)
  //   cout << p[0] << ", " << p[1] << '\n';
  // cout << endl;
  // cout << "Tris: " << endl;
  // for (const auto t : tris)
  //   cout << t[0] << " " << t[1] << " " << t[2] << '\n';
  // cout << endl << endl;


  // non-convex case
  bpoints = vector<Point2D> { {0,0}, {2,0}, {2,2}, {1,2}, {1,1}, {0,1} };
  ipoints = vector<Point2D> { {0.5, 0.5}, {1, 0.5}, {1.5, 0.5}, {1.5, 1}, {1.5, 1.5}};
  all_points = bpoints;
  all_points.insert(all_points.end(), ipoints.begin(), ipoints.end());

  tris = triangulate_domain(&all_points[0], (uint)bpoints.size(),
			    (uint) all_points.size(), vdist);

  cout << "--- Non-convex case: ---" << endl<<endl;
  cout << "Points: " << endl;
  for (const auto p : all_points)
    cout << p[0] << ", " << p[1] << '\n';
  cout << endl;
  cout << "Tris: " << endl;
  for (const auto t : tris)
    cout << t[0] << " " << t[1] << " " << t[2] << '\n';
  cout << endl << endl;
}

// ----------------------------------------------------------------------------
void test_volume_tesselation()
// ----------------------------------------------------------------------------
{
  // Tesselate regular prism
  const vector<Point3D> corners { {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, 
                                {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1} };
  const vector<Segment> edges   { {0, 1}, {1, 2}, {2, 3}, {3, 0},  // bottom face
                                {4, 5}, {5, 6}, {6, 7}, {7, 4},  // top face
                                {0, 4}, {1, 5}, {2, 6}, {3, 7}};  // "vertical" edges
  const vector<FaceLoop> faces { { {0, 1, 2, 3}, false}, // bottom (zmin) face
                                 { {4, 5, 6, 7}, true}, // top (zmax) face
                                 { {0, 9, 4, 8}, true}, // front (ymin) face
                                 { {2, 11, 6, 10}, true}, // back (ymax) face
                                 { {8, 7, 11, 3}, true}, // left (xmin) face
                                 { {1, 10, 5, 9}, true}}; // right (xmax) face
  
  const SimpleVolumeType v {};
  SimplePolyhedron spoly {corners, edges, faces, v};

  cout << "Polyhedron:\n" << endl;
  cout << spoly << endl << endl;

  const double vdist = 0.3;
  spoly.tesselate(vdist);

  cout << "Writing wirefame: " << endl;
  spoly.writeTesselatedOutline(cout);
}
};
