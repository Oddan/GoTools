#ifndef _CURVESPACINGFUN_H
#define _CURVESPACINGFUN_H

#include <utility>
#include <functional>

#include "GoTools/geometry/ParamCurve.h"
#include "find_root.h"

namespace {
  typedef std::function<double(const double* const, int)> RnToRFunction;
};

namespace Go {

std::tuple<RnToRFunction, RnToRnFunction> make_curve_spacing_fun(const Go::ParamCurve& c);

};

#endif // _CURVESPACINGFUN_H
