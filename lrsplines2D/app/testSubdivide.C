#include <iostream>
#include <fstream>
#include <locale>
#include <iostream>
#include <iterator> // for ostream_iterator
#include <vector>
#include <chrono>
//#include <tuple>
//#include <utility>
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/PlotUtils.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"

using namespace Go;
using namespace std;


int main(int argc, char* argv[])
{
  setlocale(LC_ALL, "en_US.UTF.8");

  ifstream is(argv[1]);
  if (!is) {
    wcout << L"File not found.\n";
    return 0;
  }
  ObjectHeader header;
  header.read(is);
  
  LRSplineSurface lrs(is);

  //plot_mesh(lrs.mesh());

  auto result = lrs.subdivideIntoSimpler();
  
  ofstream os("result.g2");
  for (size_t i = 0; i != result.size(); ++i) {
    (*result[i]).writeStandardHeader(os);
    (*result[i]).write(os);
  }
  os.close();
  
}; // main


