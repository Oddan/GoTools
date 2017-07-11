#ifndef _COMMON_DEFS_H
#define _COMMON_DEFS_H

#include <array>
#include <vector>

namespace TesselateUtils {

  using uint = unsigned int;
  using Point2D = std::array<double,2>;
  struct ValAndDer {double val; std::vector<Point2D> der;};
}

#endif
