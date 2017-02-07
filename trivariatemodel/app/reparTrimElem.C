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
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/ftVolumeTools.h"
#include "GoTools/trivariatemodel/VolumeModelFileHandler.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

using namespace Go;
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

//===========================================================================
//
/// Description: 
///
///              
///  
/// 
///
/// Input/Output: 
///
///               
/// 
/// Note:       
///
///             
///
//   
//===========================================================================

int main( int argc, char* argv[] )
{
  if (argc != 3)
    {
      cout << "Usage: " << "<infile> " << " <file type> "<< endl;
      exit(-1);
    }
  std::string infile(argv[1]);
  int file_type = atoi(argv[2]);

  shared_ptr<ftVolume> curr_vol;
  int degree = 3;
  if (file_type == 1)
    {
      // The tolerances must be set according to the properties of the model.
      // The neighbour tolerance must be smaller than the smallest entity in the
      // model, but larger than the largest gap.
      // The gap tolerance must be smaller than the neighbour tolerance
      double gap = 0.0001; 
      double neighbour = 0.001; 
      double kink = 0.01;
      double approxtol = 0.001;

      CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

      ifstream is(infile);
      CompositeModel *model = factory.createFromG2(is);

      shared_ptr<SurfaceModel> sfmodel = 
	shared_ptr<SurfaceModel>(dynamic_cast<SurfaceModel*>(model));
      if (!sfmodel.get())
	{
	  std::cout << "No input model read" << std::endl;
	  exit(-1);
	}
 
      if (sfmodel->nmbBoundaries() > 0)
	{
	  std::cout << "Not a brep solid. Consider increasing the neighbour tolerance" << std::endl;
	  exit(-1);
	}
      
      bool isOK = sfmodel->checkShellTopology();
      std::cout << "Shell topology: " << isOK << std::endl;

 
      // Select the first volume and pass through all elements and check if
      // they intersect the non-boundary trimming surface
      curr_vol = shared_ptr<ftVolume>(new ftVolume(sfmodel));

      std::ofstream out_file("volmodel.g22");
      VolumeModelFileHandler filehandler;
      filehandler.writeStart(out_file);
      filehandler.writeHeader("Test ftVolume", out_file);
      filehandler.writeVolume(curr_vol, out_file);
      filehandler.writeEnd(out_file);

      // Number of elements in underlying volume
      SplineVolume* curr_under = curr_vol->getVolume()->asSplineVolume();
      curr_under->raiseOrder(2, 2, 2);
      int nn = 8;
      for (int dir=0; dir<3; ++dir)
	{
	  double tstart = curr_under->startparam(dir);
	  double tend = curr_under->endparam(dir);
	  double tdel = (tend - tstart)/(double)(nn+1);
	  vector<double> newknots(nn);
	  for (int ka=0; ka<nn; ++ka)
	    newknots[ka] = tstart + (ka+1)*tdel;
	  curr_under->insertKnot(dir, newknots);
	}
    }
  else
    {
      VolumeModelFileHandler filehandler;
      curr_vol = filehandler.readVolume(infile.c_str());
      SplineVolume* curr_under = curr_vol->getVolume()->asSplineVolume();
      for (int dir=0; dir<3; ++dir)
	{
	  int num_el = curr_under->numElem(dir);
	  if (num_el < 10)
	    {
	      double tstart =  curr_under->startparam(dir);
	      double tend =  curr_under->endparam(dir);
	      double del1 = (tend-tstart)/10.0;
	      vector<double> knots;
	      curr_under->basis(dir).knotsSimple(knots);
	      vector<double> newknots;
	      for (size_t kj=1; kj<knots.size(); ++kj)
		{
		  double del2 = (knots[kj] - knots[kj-1]);
		  if (del2 < del1)
		    continue;
		  int nmb = std::max(1, (int)(del2/del1));
		  double del3 = del2/(double)(nmb+1);
		  int kr;
		  double par;
		  for (kr=0, par=knots[kj-1]+del3; kr<nmb; ++kr, par+=del3)
		    newknots.push_back(par);
		}

	      curr_under->insertKnot(dir, newknots);
	    }
	}
     }

  SplineVolume* under = curr_vol->getVolume()->asSplineVolume();
  double gap = curr_vol->getTolerances().gap;
  
  int nmb_elem = under->numElem();
  std::cout << "No of elements: " << nmb_elem << std::endl;

  std::ofstream of5("tmp5.g2");
  std::ofstream of6("tmp6.g2");
  std::ofstream of7("tmp7.g2");
  for (int ki=0; ki<nmb_elem; ++ki)
    {
      int elem_stat = curr_vol->ElementBoundaryStatus(ki);
      std::cout << "Boundary status, element " << ki+1 << ": " << elem_stat << std::endl;

      if (elem_stat == 1)
	{
	  // Element intersects trimming surface
	  // Split element with trimming shell
	  vector<shared_ptr<ftVolume> > sub_elem;
	  vector<int> is_inside; // Equal 1 if the sub element is inside
	  // the trimmed volume
	  curr_vol->splitElementByTrimSfs(ki, gap, sub_elem, is_inside);

	  std::ofstream of4("tmp4.g2");
	  for (size_t kj=0; kj<sub_elem.size(); ++kj)
	    {
	      shared_ptr<SurfaceModel> mod = sub_elem[kj]->getOuterShell();
	      int nmb = mod->nmbEntities();
	      for (int kr=0; kr<nmb; ++kr)
		{
		  shared_ptr<ParamSurface> sf = mod->getSurface(kr);
		  sf->writeStandardHeader(of4);
		  sf->write(of4);
		}
	    }
	  int stop_break = 1;

	  // if (sub_elem.size() < 2)
	  //   continue;

 	  // Check if the remaining element is hexagonal
	  for (size_t kj=0; kj<sub_elem.size(); ++kj)
	    {
	      if (!is_inside[kj])
		continue;

	      bool regular = sub_elem[kj]->isRegularized();
	      std::cout << "Sub element nr " << kj+1 << ": " << regular << std::endl;
	      if (regular)
		{
		  if (false)
		    {
		  // Create non-trimmed parameter element
		  shared_ptr<ParamVolume> reg_vol = 
		    sub_elem[kj]->getRegParVol(degree);
		  if (reg_vol.get())
		    {
		      reg_vol->writeStandardHeader(of5);
		      reg_vol->write(of5);
		    }

		  // Create non-trimmed element
		  sub_elem[kj]->untrimRegular(degree);
		  shared_ptr<ParamVolume> tmp_vol = sub_elem[kj]->getVolume();
		  tmp_vol->writeStandardHeader(of6);
		  tmp_vol->write(of6);
		    }
		}
	      else
		{
		  std::ofstream of4_2("tmp4_2.g2");
		  shared_ptr<SurfaceModel> mod = sub_elem[kj]->getOuterShell();
		  int nmb = mod->nmbEntities();
		  for (int kr=0; kr<nmb; ++kr)
		    {
		      shared_ptr<ParamSurface> sf = mod->getSurface(kr);
		      sf->writeStandardHeader(of4_2);
		      sf->write(of4_2);
		      sf->writeStandardHeader(of7);
		      sf->write(of7);
		    }

		  std::ofstream of_mod("tmp_mod.g22");
		  VolumeModelFileHandler filewrite;
		  filewrite.writeStart(of_mod);
		  filewrite.writeHeader("Irregular volume", of_mod);
		  filewrite.writeVolume(sub_elem[kj], of_mod);
		  filewrite.writeEnd(of_mod);
		  

		  std::cout << "Number of surfaces: " << sub_elem[kj]->getOuterShell()->nmbEntities() << std::endl;
		}
	    }
	}
      else if (elem_stat == 2)
	{
	  SplineVolume *vol = curr_vol->getVolume()->asSplineVolume();
	  if (vol)
	    {
	      // Fetch parameter values surrounding the specified element
	      double elem_par[6];
	      vol->getElementBdPar(ki, elem_par);
	      
	      // Create an ftVolume entity corresponding to the element. 
	      // First create underlying SplineVolume
	      shared_ptr<ParamVolume> elem_vol(vol->subVolume(elem_par[0], elem_par[2],
							      elem_par[4], elem_par[1],
							      elem_par[3], elem_par[5]));
	      elem_vol->writeStandardHeader(of6);
	      elem_vol->write(of6);
	    }

	}
    }
}
