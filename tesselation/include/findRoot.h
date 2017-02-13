#ifndef _FINDROOT_H
#include _FINDROOT_H

#include <vector>
#include <functional>

namespace Go {

// ----------------------------------------------------------------------------

// Class representing the value and jacobian of a multi-parameter, multi-valued function
struct ValAndJac
{
  std::vector<double> value; // of size N
  std::vector<double> jacobian; // of size NxN (column-wise storage)
};
  
// ----------------------------------------------------------------------------

// Function from Rn to Rn, returning a ValAndJac.  First argument is a pointer
// to the function arguments, second argument is the number of function arguments
typedef std::function<ValAndJac(const double* const, unsigned int)> RnToRnFunction;

// ----------------------------------------------------------------------------

// Function searching for root
std::vector<double> FindRoot(const RnToRnFunction& fun,
			     const double* const start,
			     unsigned int dim);
  
}; //end namespace Go

#endif // _FINDROOT_H
