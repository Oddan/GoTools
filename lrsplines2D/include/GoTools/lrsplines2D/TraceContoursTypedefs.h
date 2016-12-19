#ifndef _TRACECONTOURSTYPEDEFS_H
#define _TRACECONTOURSTYPEDEFS_H

#include <vector>
#include <functional>
#include <algorithm>
#include "GoTools/geometry/SplineCurve.h"

namespace Go
{
  // Common typedefs used both by SSurfTraceIsocontours and LRTraceIsocontours
  
  using CurvePtr = std::shared_ptr<const SplineCurve>;

  // The first curve pointer of the pair represents a 2D curve in the parameter
  // domain of the investigated LR spline function.  The second curve pointer
  // represents the corresponding 3D curve, if requested.
  using CurveVec = std::vector<std::pair<CurvePtr, CurvePtr>>;


  // ----------------------------------------------------------------------------
  // The following utility template function ought not to have any overhead in
  // optimized mode, due to the 'move' semantics of STL vectors.
  template<typename T, typename R>
  std::vector<R> apply_transform(const std::vector<T>& input,
				 const std::function<R(T)>& trans_fun) 
  // ----------------------------------------------------------------------------
  {
    std::vector<R> tmp(input.size());
    std::transform(input.begin(), input.end(), tmp.begin(), trans_fun);
    return tmp;
  }

};

#endif
