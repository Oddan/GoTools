#include "ParametricObjectEnergyFunctor.h"

namespace TesselateUtils {

// ----------------------------------------------------------------------------
template<typename ParamObj, int dim> double
ParametricObjectEnergyFunctor<ParamObj, dim>::operator()(const double* const arg) const 
// ----------------------------------------------------------------------------
{
  if (!use_cached(arg))
    update_cache(arg);
  return cached_result_.val;
}

// ----------------------------------------------------------------------------
template<typename ParamObj, int dim> double
ParametricObjectEnergyFunctor<ParamObj, dim>::grad(const double* arg,
                                                   double* grad) const
// ----------------------------------------------------------------------------
{
  if (!use_cached(arg))
    update_cache(arg);
  const double* const dp = (const double* const)&cached_result_.der[0];
  copy(dp, dp + dim * n_, grad);
}

// ----------------------------------------------------------------------------
template<> double
ParametricObjectEnergyFunctor<Go::ParamCurve, 1>::minPar(int i) const
// ----------------------------------------------------------------------------  
{
  return pobj_.startparam();
}

// ----------------------------------------------------------------------------
template<> double 
ParametricObjectEnergyFunctor<Go::ParamCurve, 1>::maxPar(int i) const
// ----------------------------------------------------------------------------  
{
  return pobj_.endparam();
}

// ----------------------------------------------------------------------------    
  template<> bool ParametricObjectEnergyFunctor<Go::ParamCurveuse_cached(const double* const arg) const
// ----------------------------------------------------------------------------
{

}
  
// ----------------------------------------------------------------------------    
template<> void update_cache(const double* const arg) const
// ----------------------------------------------------------------------------  
{
  
}

  
}; // end namespace TesselateUtils



