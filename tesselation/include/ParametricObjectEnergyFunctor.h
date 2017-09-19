#ifndef _PARAMETRIC_OBJECT_ENERGY_FUNCTOR_H
#define _PARAMETRIC_OBJECT_ENERGY_FUNCTOR_H

#include <vector>
#include "common_defs.h"

namespace TesselateUtils {

// ============================================================================
template<typename ParamObj, int dim> // can be ParamCurve (dim=1), 
class ParametricObjectEnergyFunctor  // ParamSurface (dim=2) or ParamVolume (dim=3)
// ============================================================================
{
public:
  ParametricObjectEnergyFunctor(const ParamObj& pobj,
                                double radius,
                                unsigned int num_points)
    : pobj_(pobj), radius_(radius), np_(num_points) {}

  double operator() (const double* const arg) const;
  void grad(const double* arg, double* grad) const;
  double minPar(int i) const;
  double maxPar(int i) const;
  
private:
  const ParamObj& pobj_; // parametric object
  const double radius_; // radius of energy function
  const np_; // number of unknown points

  mutable ValAndDer<dim> cached_result_;
  mutable std::vector<double> cached_arg_;

  bool use_cached(const double* const arg) const;
  void update_cache(const double* const arg) const;
  
}; // end ParametricObjectEnergyFunctor

  
} // end namespace TesselateUtils


#endif
