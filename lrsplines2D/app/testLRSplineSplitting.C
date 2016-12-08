#include <iostream>
#include <fstream>
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/PlotUtils.h"


using namespace Go;
using namespace std;

// data/bicubic_7_7_ref_u_lr.g2 - very simple candidate


int main(int argc, char* argv[]) {
  setlocale(LC_ALL, "en_US.UTF.8");
  ifstream is(argv[1]);
  if (!is) {
    wcout << L"File not found.\n";
    return 0;
  }
  ObjectHeader header;
  header.read(is);
  
  LRSplineSurface lrs(is);

  plot_mesh(lrs.mesh());
  
  return 0;
}
  
