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
template<typename ParamObj, int dim> bool
ParametricObjectEnergyFunctor<ParamObj, dim>::use_cached(const double* const arg) const
// ----------------------------------------------------------------------------
{
  return (cached_arg_.size() > 0) &&
         (std::equal(arg, arg + dim * np_, &cached_arg_[0]))
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
template<> void
ParametricObjectEnergyFunctor<Go::ParamCurve, 1>::update_cache(const double* const arg) const
// ----------------------------------------------------------------------------  
{
  if (cached_arg_.empty())
    cached_arg_.resize(dim * np_);

  std::copy(arg, arg + dim * np_, &cached_arg_[0]);

  cached_result_ = 


template<typename PointXD>  
struct ValAndDer {
  double val;
  std::vector<PointXD> der;

  void reset(uint num_der);
};


}
  
}; // end namespace TesselateUtils



