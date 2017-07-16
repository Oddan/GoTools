#include <iostream> // for debugging
#include <array>
#include "SimplePolyhedronTesselation.h"

using namespace std;

namespace TesselateUtils
{
// ----------------------------------------------------------------------------
template<>
void SimplePolyhedron::tesselate()
// ----------------------------------------------------------------------------
{
  cout << "Tesselate() called." << endl;
}

// ----------------------------------------------------------------------------
template<> 
void SimplePolyhedron::compute_tesselation(const std::array<Point3D, 2>& boundary,
                                           const Segment& edge,
                                           std::vector<Point3D>& ipoints)
// ----------------------------------------------------------------------------
{
  
}

// ----------------------------------------------------------------------------
template<> 
void SimplePolyhedron::compute_tesselation(const std::vector<Point3D>& boundary,
                                           const FaceLoop& face,
                                           std::vector<Point3D>& ipoints,
                                           std::vector<Triangle>& triangles)
// ----------------------------------------------------------------------------
{

}

// ----------------------------------------------------------------------------  
template<> 
void SimplePolyhedron::compute_tesselation(const std::vector<Point3D>& bpoints,
                                           const std::vector<Triangle>& btris,
                                           const VolumeType& volume,
                                           std::vector<Point3D>& ipoints,
                                           std::vector<Tet>& tets)
// ----------------------------------------------------------------------------
{

}

// ----------------------------------------------------------------------------  
template<> array<Point3D, 2> SimplePolyhedron::edgeCorners(uint edge_ix) const
// ----------------------------------------------------------------------------
{
  
}
  
// ----------------------------------------------------------------------------
template<> vector<Point3D> SimplePolyhedron::faceBoundaryPoints(uint face_ix) const
// ----------------------------------------------------------------------------
{

}

  
  
}; // end namespace TesselateUtils
