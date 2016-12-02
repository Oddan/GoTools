#include <locale>
#include <iostream>
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/PlotUtils.h"

using namespace Go;

namespace {
  int insertSegment(Go::Mesh2D& m,
		    Direction2D d,
		    double pval,
		    int start_ix,
		    int end_ix,
		    int mult);
} // end anonymous namespace

int main(int argc, char* argv[])
{
  setlocale(LC_ALL, "en_US.UTF.8");
  
  // Defining initial knotvectors
  const double xvec[] = {0, 1, 2, 3, 4, 5};
  const double yvec[] = {0, 1, 2, 3, 4};

  // base, tensor-product mesh
  Mesh2D m = Mesh2D(xvec, xvec + 6, yvec, yvec + 5);

  // inserting modifications (no longer tensor product)
  insertSegment(m, YFIXED, 2.5, 2, 5, 1);
  insertSegment(m, XFIXED, 3.5, 1, 3, 1);
  
  plot_mesh(m);

}

// ----------------------------------------------------------------------------
namespace {
  int insertSegment(Go::Mesh2D& m,
		    Direction2D d,
		    double pval,
		    int start_ix,
		    int end_ix,
		    int mult)
  {
    int l_ix = m.insertLine(d, pval);
    m.setMult(d, l_ix, start_ix, end_ix, mult);
    return l_ix;
  }
} // end anonymous namespace

  
