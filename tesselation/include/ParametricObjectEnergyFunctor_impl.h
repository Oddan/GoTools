#include <iostream> // @@ for debugging only
namespace TesselateUtils {

// ----------------------------------------------------------------------------
ClippedGrid<1> inline
ParamCurveEnergyFunctionTraits::compute_cgrid(const ParamObj pobj,
                                              const BoundaryPolytope& boundary,
                                              const double radius)
// ----------------------------------------------------------------------------  
{
  // ClippedGrid not used in the curve case, so we just return a dummy one.
  return ClippedGrid<1>();
}

// ----------------------------------------------------------------------------
ClippedGrid<2> inline
ParamSurfaceEnergyFunctionTraits::compute_cgrid(const ParamObj pobj,
                                                const BoundaryPolytope& boundary,
                                                const double radius)
// ----------------------------------------------------------------------------  
{
  std::cout << "Warning in ParamSurfaceEnergyFunctionTraits::compute_grid: "
    "dummy implementation." << std::endl;

  // @@@ dummy implementation for now
  auto krull = clip_grid_polygon_2D(boundary.bpoints,
                                    boundary.num_bpoints,
                                    radius,
                                    80, 80); // @@ 80 hardcoded here
  fill(krull.type.begin(), krull.type.end(), UNDETERMINED);
  return krull;
}

// ----------------------------------------------------------------------------
ClippedGrid<3> inline
ParamVolumeEnergyFunctionTraits::compute_cgrid(const ParamObj pobj,
                                               const BoundaryPolytope& boundary,
                                               const double radius)
// ----------------------------------------------------------------------------  
{
  assert(false);
  return ClippedGrid<3>();
}
  
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
template<> void inline
ParamCurveEnergyFunctor::update_cache(const double* const arg) const
// ----------------------------------------------------------------------------  
{
  if (cached_arg_.empty())
    cached_arg_.resize(Dim * np_);

  std::copy(arg, arg + Dim * np_, &cached_arg_[0]);

  cached_result_ = parametric_curve_energy(pobj_,
                                           minPar(0),
                                           maxPar(0),
                                           &cached_arg_[0],
                                           np_,
                                           radius_);
}

// ----------------------------------------------------------------------------    
template<> void inline
ParamSurfaceEnergyFunctor::update_cache(const double* const arg) const
// ----------------------------------------------------------------------------  
{
  if (cached_arg_.empty())
    cached_arg_.resize(Dim * np_);

  std::copy(arg, arg + Dim * np_, &cached_arg_[0]);

  cached_result_ = parametric_surf_energy(pobj_,
                                          boundary_.bpoints,
                                          boundary_.num_bpoints,
                                          (const Point2D* const) &cached_arg_[0],
                                          np_,
                                          radius_,
                                          &cgrid_);
}
  
}; // end namespace TesselateUtils
