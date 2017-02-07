/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#ifndef _CREATETRIMVOLUME_H
#define _CREATETRIMVOLUME_H

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/utils/DirectionCone.h"

namespace Go
{
  class ftSurface;
  class ftVolume;
  class ParamVolume;
  class ParamSurface;

  class CreateTrimVolume
  {
  public:
    /// Constructor
    CreateTrimVolume(shared_ptr<SurfaceModel> model);

    /// Destructor
    ~CreateTrimVolume();

    shared_ptr<ftVolume> fetchOneTrimVol();

  private:
    shared_ptr<SurfaceModel> model_;

    vector<vector<shared_ptr<ftSurface> > > face_grp_;
    vector<shared_ptr<ParamSurface> > under_sf_;
    vector<BoundingBox> bbox_;
    vector<DirectionCone> cone_;
    vector<double> sfsize_;

    void identifyBoundaryFaces(std::vector<std::pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> > >& side_sfs);

    void extractMaxSet(std::vector<shared_ptr<ftSurface> >& bd_faces,
		       std::vector<shared_ptr<ftSurface> >& trim_faces);

    shared_ptr<ftVolume> 
      createTrimVolume(shared_ptr<ParamVolume> vol, 
		       std::vector<std::pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> > >& side_sfs);
    
    void identifyInnerTrim(std::vector<shared_ptr<ftSurface> >& bd_faces,
			   std::vector<shared_ptr<ftSurface> >& trim_faces);

    void computeGroupInfo();

    void findSideSfs(double tol, double angtol,
		     std::vector<std::pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> > >& side_sfs);

    int
      checkCandPair(Point vec, shared_ptr<ParamSurface> sf1, int bd_type1, 
		    BoundingBox& box1, shared_ptr<ParamSurface> sf2,
		    int bd_type2, BoundingBox& box2, 
		    double tol);

    void
      oneSideSf(int bd_type, std::vector<int>& face_group_ix, 
		Point bd_vec, Point dir, double tol, double angtol,
		std::pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> >& side_sf);
    void extendSurfaces(std::vector<std::pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> > >& side_sfs);

    void trimSideSurfaces(std::vector<std::pair<shared_ptr<ftSurface>, shared_ptr<ParamSurface> > >& side_sfs);
    void refineInSharpEdges(shared_ptr<ParamVolume>& vol);

    bool checkIsoPar(shared_ptr<ParamSurface> surf,
		     shared_ptr<ParamVolume> vol,
		     int pardir, double parval, double tol);
  };
} // namespace Go


#endif // _CREATETRIMVOLUME_H
