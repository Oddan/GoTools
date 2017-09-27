#ifndef _PARAM_SURF_ENERGY_FUNCTOR_H
#define _PARAM_SURF_ENERGY_FUNCTOR_H

#include "common_defs.h"

namespace TesselateUtils
{

class ParamSurfEnergyFunctor
{
public:
  ParamSurfEnergyFunctor(const shared_ptr<const ParamSurface> surf,
                         const Point2D* bpoints,
                         const uint num_bpoints,
                         const uint num_ipoints,
                         const double radius);

  double operator() (const double* const arg) const;
  void grad(const double* arg, double* grad) const;
  double minPar(int n) const;
  double maxPar(int n) const;
  
private:
  bool use_cached(const double* const arg) const;
  void update_cache(const double* const arg) const;

  const shared_ptr<const ParamSurface> surf_,
  const Point2D* const poly_;
  const uint nc_; // num polygon corners
  const uint ni_; // num interior points
  const double r_; // radius
  const std::array<double, 4> bbox_; // bounding box
  const ClippedGrid<2> cgrid_; // precomputed classifiation of subdivided domain
                               // parts, according to their relationship with
                               // the polygon boundary (inside, outside, etc.).
                               // used to improve computational efficiency.

  mutable ValAndDer<Point2D> cached_result_;
  mutable vector<double> cached_arg_;
  
};
  
  
}; // end namespace TesselateUtils

#endif
