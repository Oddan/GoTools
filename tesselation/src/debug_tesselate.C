#include "debug_tesselate.h"
#include "tesselate_utils.h"
#include "polyhedral_energies.h"
#include <fstream>

using namespace TesselateUtils;
using namespace std;

// ----------------------------------------------------------------------------
void sample_polyhedron_energy(double xmin, double xmax, int num_x,
                              double ymin, double ymax, int num_y,
                              double zmin, double zmax, int num_z,
                              const Point3D* const bpoints,
                              const unsigned int num_bpoints,
                              const Triangle* const btris,
                              const unsigned int num_btris,
                              const double vdist,
                              string filename)
// ----------------------------------------------------------------------------
{
  vector<Point3D> sample_points = generate_grid_3D(Point3D {xmin, ymin, zmin},
                                                   Point3D {xmax, ymax, zmax},
                                                   (uint) num_x,
                                                   (uint) num_y,
                                                   (uint) num_z,
                                                   false);

  // computing values
  vector<ValAndDer<Point3D>> vals;
  transform(sample_points.begin(), sample_points.end(), back_inserter(vals),
            [&] (const Point3D& p) { return polyhedron_energy(bpoints, num_bpoints,
                                                              btris, num_btris,
                                                              &p, 1, vdist); });
  
  // saving results
  ofstream os(filename.c_str());

  os << "A - samples: " << endl;
  for (auto p : sample_points)
    os << p[0] << " " << p[1] << " " << p[2] << '\n';
  
  os << '\n' << "B - values: " << endl;
  for (auto v : vals)
    os << v.val << ' ',

  os << '\n' << "C - derivatives: " << endl;
  for (auto v : vals) {
    for (auto d : v.der)
      os << d[0] << ' ' << d[1] << ' ' << d[2] << '\n';
    os << '\n\n';
  }
  
  os.close();
}
