namespace TesselateUtils {

// ----------------------------------------------------------------------------
template<typename Traits> double
ParametricObjectEnergyFunctor<Traits>::operator()(const double* const arg) const 
// ----------------------------------------------------------------------------
{
  if (!use_cached(arg))
    update_cache(arg);
  return cached_result_.val;
}

// ----------------------------------------------------------------------------
template<typename Traits> void
ParametricObjectEnergyFunctor<Traits>::grad(const double* arg, double* grad) const
// ----------------------------------------------------------------------------
{
  if (!use_cached(arg))
    update_cache(arg);
  const double* const dp = (const double* const)&cached_result_.der[0];
  std::copy(dp, dp + Dim * np_, grad);
}

// ----------------------------------------------------------------------------  
template<typename Traits> bool
ParametricObjectEnergyFunctor<Traits>::use_cached(const double* const arg) const
// ----------------------------------------------------------------------------
{
  return (cached_arg_.size() > 0) &&
         (std::equal(arg, arg + Dim * np_, &cached_arg_[0]));
}

// ----------------------------------------------------------------------------
template<typename Traits> double
ParametricObjectEnergyFunctor<Traits>::minPar(int n) const
// ----------------------------------------------------------------------------  
{
  return bbox_[(n % 2) * 2];
}
  
// ----------------------------------------------------------------------------
template<typename Traits> double
ParametricObjectEnergyFunctor<Traits>::maxPar(int n) const
// ----------------------------------------------------------------------------
{
  return bbox_[(n % 2) * 2 + 1];
}

// ----------------------------------------------------------------------------    
template<> void
ParamCurveEnergyFunctor::update_cache(const double* const arg) const
// ----------------------------------------------------------------------------  
{
  if (cached_arg_.empty())
    cached_arg_.resize(1 * np_);

  std::copy(arg, arg + 1 * np_, &cached_arg_[0]);

  cached_result_ = parametric_curve_energy(pobj_,
                                           minPar(0),
                                           maxPar(0),
                                           &cached_arg_[0],
                                           np_,
                                           radius_);
}
  
}; // end namespace TesselateUtils
