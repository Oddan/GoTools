#include "polyhedral_energies.h"
#include "interpoint_distances.h"
#include "tesselate_utils.h"
#include "common_defs.h"

using namespace std;
using namespace TesselateUtils;

namespace {

// first entry in result array is the energy, the second the derivative
array<double, 2> energy(double dist, double R);

ValAndDer internal_energy(const Point2D* const ipoints,
			  const unsigned int num_ipoints,
			  double vdist);
ValAndDer boundary_energy(const Point2D* const bpoints,
			  const unsigned int num_bpoints,
			  const Point2D* const ipoints,
			  const unsigned int num_ipoints,
			  double vdist);
ValAndDer mirror_energy(const Point2D* const bpoints,
			const unsigned int num_bpoints,
			const Point2D* const ipoints,
			const unsigned int num_ipoints,
			double vdist);

}; // end anonymous namespace 

namespace TesselateUtils {

// ----------------------------------------------------------------------------
ValAndDer polygon_energy(const Point2D* const bpoints,
			 const unsigned int num_bpoints,
			 const Point2D* const ipoints,
			 const unsigned int num_ipoints,
			 double vdist)
// ----------------------------------------------------------------------------
{
  
  // compute internal energy (potential energy between internal points)
  const ValAndDer E_int = internal_energy(ipoints, num_ipoints, vdist);

  // compute boundary energy (energy from interaction between internal points
  // and boundary points)
  const ValAndDer E_bnd = boundary_energy(bpoints, num_bpoints,
					  ipoints, num_ipoints, vdist);

  // compute mirroring energy (energy from interaction between internal points
  // and their virtual neighbors across boundary)
  const ValAndDer E_mir = mirror_energy(bpoints, num_bpoints, 
					ipoints, num_ipoints, vdist);

  // Adding up components and returning results
  ValAndDer E_tot = E_int;

  E_tot.val += E_bnd.val; 
  E_tot.der += E_bnd.der;
  
  E_tot.val += E_mir.val;
  E_tot.der += E_mir.der;

  return E_tot;
}

}; // end namespace Go


namespace {

// ----------------------------------------------------------------------------
ValAndDer internal_energy(const Point2D* const ipoints,
			  const unsigned int num_ipoints,
			  double vdist)
// ----------------------------------------------------------------------------
{
  const auto dists = interpoint_distances(ipoints, num_ipoints, vdist, false);

  vector<array<double,2>> energies(dists.size());
  
  transform(dists.begin(), dists.end(), energies.begin(), 
	    [vdist](const DistanceEntry& d) {return energy(d.dist, vdist);});

  do_stuff_here;


  return ValAndDer();
}
// ----------------------------------------------------------------------------
ValAndDer boundary_energy(const Point2D* const bpoints,
			  const unsigned int num_bpoints,
			  const Point2D* const ipoints,
			  const unsigned int num_ipoints,
			  double vdist)
// ----------------------------------------------------------------------------
{
  return ValAndDer();
}

// ----------------------------------------------------------------------------
ValAndDer mirror_energy(const Point2D* const bpoints,
			const unsigned int num_bpoints,
			const Point2D* const ipoints,
			const unsigned int num_ipoints,
			double vdist)
// ----------------------------------------------------------------------------
{
  return ValAndDer();
}

// ----------------------------------------------------------------------------
array<double, 2> energy(double dist, double R)
// ----------------------------------------------------------------------------
{
  const double tmp = max(R-dist, double(0));
  const double tmp2 = tmp*tmp;
  return {tmp2 * tmp2, -4 * tmp*tmp*tmp}; // energy and derivative
}

}; // end namespace Go
