#include "ParamSurfEnergyFunctor.h"

namespace TesselateUtils {

// ----------------------------------------------------------------------------  
ParamSurfEnergyFunctor::
ParamSurfEnergyFunctor(const shared_ptr<const ParamSurface> surf,
                       const Point2D* bpoints,
                       const uint num_bpoints,
                       const uint num_ipoints,
                       const double radius)
// ----------------------------------------------------------------------------  
  : surf_(surf), poly_(bpoints), nc_(num_bpoints),
    ni_(num_ipoints), r_(radius), bbox_(bounding_box_2D(bpoints, num_bpoints)),
    cgrid_(clip_grid_polygon_2D(bpoints, num_bpoints, radius, 80, 80, surf)
           // @@ 80 hard-coded above
{}

// ----------------------------------------------------------------------------  
double ParamSurfEnergyFunctor::operator() (const double* const arg) const
// ----------------------------------------------------------------------------  
{
  if (!use_cached(arg))
    update_cache(arg);
  return cached_result_.val;
}
  
// ----------------------------------------------------------------------------
void ParamSurfEnergyFunctor::grad(const double* arg, double* grad) const
// ----------------------------------------------------------------------------
{
  if (!use_cached(arg))
    update_cache(arg);
  const double* const dp = (const double* const)&cached_result_der[0];
  copy(dp, dp + 2 * ni_, grad);
}

// ----------------------------------------------------------------------------
double ParamSurfEnergyFunctor::minPar(int n) const
// ----------------------------------------------------------------------------  
{
  return bbox_[(n % 2) * 2];
}
  
// ----------------------------------------------------------------------------  
double ParamSurfEnergyFunctor::maxPar(int n) const
// ----------------------------------------------------------------------------
{
  return bbox_[(n % 2) * 2 + 1];
}

// ----------------------------------------------------------------------------  
bool ParamSurfEnergyFunctor::use_cached(const double* const arg) const
// ----------------------------------------------------------------------------
{
  return (cached_arg_.size() > 0) && (std::equal(arg,
                                                 arg + 2 * ni_,
                                                 &cached_arg_[0]));
}

// ----------------------------------------------------------------------------  
void ParamSurfEnergyFunctor::update_cache(const double* const arg) const
// ----------------------------------------------------------------------------  
{
  if (cached_arg_.empty())
    cached_arg_.resize(2 * ni_);

  copy(arg, arg + 2 * ni_, &cached_arg_[0]);

  cached_result_ = parametric_surf_energy(surf_, poly_, nc_, 
                                          (const Point2D* const) &cached_arg_[0],
                                          ni_, r_, &cgrid_);
}
  
}; // end namespace TesselateUtils
