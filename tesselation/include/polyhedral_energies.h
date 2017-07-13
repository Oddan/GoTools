#ifndef _POLYHEDRAL_ENERGIES_H
#define _POLYHEDRAL_ENERGIES_H

#include <utility>
#include "common_defs.h"

namespace TesselateUtils {

  ValAndDer polygon_energy(const Point2D* const bpoints,
			   const unsigned int num_bpoints,
			   const Point2D* const ipoints,
			   const unsigned int num_ipoints,
			   const double vdist);

}; // end namespace Go

#endif
