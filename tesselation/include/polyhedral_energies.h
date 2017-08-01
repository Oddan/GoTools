#ifndef _POLYHEDRAL_ENERGIES_H
#define _POLYHEDRAL_ENERGIES_H

#include <utility>
#include "common_defs.h"

namespace TesselateUtils {

// ============================================================================
ValAndDer<Point2D> polygon_energy(const Point2D* const bpoints,
                                  const unsigned int num_bpoints,
                                  const Point2D* const ipoints,
                                  const unsigned int num_ipoints,
                                  const double vdist);
// ============================================================================


// ============================================================================
ValAndDer<Point3D> polyhedron_energy(const Point3D* const bpoints,
                                     const unsigned int num_bpoints,
                                     const Triangle* const btris,
                                     const unsigned int num_btris,
                                     const Point3D* const ipoints,
                                     const unsigned int num_ipoints,
                                     const double vdist);
// ============================================================================  
  

}; // end namespace Go

#endif
