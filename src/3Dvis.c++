/************************************************************************

   File: 3Dvis.c++

   Author: Jennifer Campbell

   Created:
   Revisions:

   Description:

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.

***********************************************************************/

#include <vtkLookupTable.h>
#include <vtkHedgeHog.h>
#include <vtkConeSource.h>


#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <float.h>

#include "jc_tools.h"

/*
NOTE: if set plot_fan to true, must also set it to true in ODFPlotter constructor. currently hardcoded.  

-noresample is forced for vis of ProbLabels because I'm really not confident doing anything but nearest neighbour resampling of the ProbLabels.  This could be revisited. 

note that the map_tract_point_scalars option works only if the scalars are between 0 and 1, and high values will be bluer, low redder.  The output of mincFibreTrack, at this point, is the real scalars, not the scaled ones, so the user needs to convert the file if this feature is desired (impractical for high number of seed points/iterations).


could do: all tracts in a scene should have the SAME lut.  they currently don't, as the lut gets initalized based on the range of scalars.  check this.


sigma_theta as a file input is only implemented for ODF, with maxima (not eigenvector, although conceivably could want a way of attaching it to vector imagedatas too)


DO: add -radius input for -tracts




for opacity <1, can't tell what is in front of what very well (shading is wrong: a vtk bug?)

note map_tract_point_scalars and interactive_level_set and movie are currently exclusive (could be combined in one Style).


**************************
**************************

on resampling: check all possible input combinations and see whether the resampling/non is OK:

DO: add resampling of eg ODF and ODFmask, vector and vectormask, etc

NOTE: directional information is assumed to be in the same space (same direction cosines) as anything it may need to be resampled -like: resampling of directions to another space is not handled here


*********************
*********************



 */

void PlotODFs(char *ODF_filename, char *mask_filename, DISPLAY_MODE display_mode, DisplayWindow *display_window, double ODF_opacity, float colour[]) //ilana vtk5 SetOpacity takes double
{

  int counter;
  ODFPlotter **ODFplotter_array;
  int x,y,z;
  gsl_vector_float *voxel_scale_factors;
  float minimum_scale_factor;
  int i;
  float scale_factor;


  //create ImageData for each input:

  ImageData *ODF_image_data = new ImageData(ODF_filename);
  ImageData *mask_image_data = new ImageData(mask_filename);



  //allocate max size this array may have to be (wasteful: could use realloc):

  ODFplotter_array=(ODFPlotter **)malloc(ODF_image_data->GetXSize()*ODF_image_data->GetYSize()*ODF_image_data->GetZSize()*sizeof(ODFPlotter *));

  ODFPlotter *ODFplotter;


  //loop through all voxels of ODF: if mask>0.5; dynamically allocate a new ODFPlotter and prepare it.  keep track of their best scale factors in a gsl_vector_float array.  keep track of how many plotting voxels there are.

  int plot_fan=false;
  ImageData *polarity;
  ImageData *full_pdf;

  if (plot_fan)
    {
  
      //for CI, read in fullpdf to plot it instead of histo_most_likely.  must match.
  
      
      full_pdf=new ImageData("full_pdf.mnc");
      polarity=new ImageData("polarity.mnc");
    }


		  

  counter=0;
  minimum_scale_factor=FLT_MAX;

  for (x=0;x<ODF_image_data->GetXSize();x++)
    {
	for (y=0;y<ODF_image_data->GetYSize();y++)
	{
	  for (z=0;z<ODF_image_data->GetZSize();z++)
	    {

	      if (mask_image_data->GetValue(x,y,z)>0.5)
		{
		  //optionally normalize the ODF first: a visualization choice:
		  //ODF_image_data->NormalizeODF(x,y,z);

		  if (plot_fan)
		    {
		      ODFplotter = new ODFPlotter(ODF_image_data,polarity,full_pdf);
		    }
		  else
		    {
		      ODFplotter = new ODFPlotter(ODF_image_data);
		    }
		  ODFplotter->SetSmoothOn();
		  ODFplotter->SetVoxel(x,y,z);
		  ODFplotter->SetOpacity(ODF_opacity);

		  if (colour[0]==-1)
		    {
		      ODFplotter->SetColourScheme(RGB);
		    }
		  else
		    {
		      ODFplotter->SetColour(colour[0],colour[1],colour[2]);
		    }


		  switch (display_mode)
		    {
		    case SubtractMinimum:
		      ODFplotter->SetDisplayMode(SubtractMinimum);
		      break;
		    case Stretch:
		      ODFplotter->SetDisplayMode(Stretch);
		      break;
		    }

		  if (counter>0)
		    {
		      ODFplotter->SetInterpolator(ODFplotter_array[0]->GetInterpolator());
		      if (plot_fan) ODFplotter->SetInterpolator2(ODFplotter_array[0]->GetInterpolator2());
		    }
		  else
		    {
		      ODFplotter->GetInterpolator()->CalculateLookupTable();
		      if (plot_fan) ODFplotter->GetInterpolator2()->CalculateLookupTable();
		    }


		  ODFplotter->CalculatePlottingPoints();

		  scale_factor=float_abs(ODFplotter->GetBestScaleFactorForVoxelFit());

		  if(scale_factor<minimum_scale_factor)
		    {
		      minimum_scale_factor=scale_factor;
		    }

		  ODFplotter_array[counter]=ODFplotter;
		  counter++;

		}

	    }

	}

    }

  //cout << "done prep\n";


  for (i=0;i<counter;i++)
    {
      
      ODFplotter_array[i]->SetScaleFactor(minimum_scale_factor);

      display_window->AddObject(ODFplotter_array[i]->CreateActor());

    }

  //DO: delete all ODFPlotters in array, then:
  free(ODFplotter_array);

  ODF_image_data->Release();
}


void Plot3DImage(char *image_filename, DisplayWindow *display_window,ThreeOrthogonalPlanesPlotter *three_planes_plotter)
{

  display_window->AddObject((vtk3DWidget **)three_planes_plotter->GetWidgets(),3);



}


void PlotTracts(char **tract_set_filename_array, int number_of_tract_sets, DisplayWindow *display_window, bool plot_tubes, bool map_tract_point_scalars, double *tract_colour, bool colour_by_direction, char *tract_scalars_name, float scalar_low_thresh, float scalar_high_thresh)
{

//DO: could make ColourByDirection be an option for many systems as well.
//ColourByDirection currently only works for ascii tracts file (in which case, it is the default for one system, if multiple they will each have one solid colour) otherwise will be one colour
//DO: get the image size for a reasonable tube radius. hardcoded to 0.3 here, which is nice for 2x2x2


  int i;
  FibreTractPlotter **fibre_tract_plotter_array;
  fibre_tract_plotter_array=(FibreTractPlotter **) malloc(number_of_tract_sets*sizeof(FibreTractPlotter *));
  


  bool free_tract_colour=false;

  
  if (number_of_tract_sets>1)
    {
      
      colour_by_direction=false; //currently overriding input
    }
  
  if (tract_colour!=NULL)
    {
      colour_by_direction=false;
    }

  if (!colour_by_direction)


    {
      if (tract_colour==NULL)
	{
	  free_tract_colour=true;
	  tract_colour=(double *) malloc(3*sizeof(double));
	}
      

   }

  //check (not nec)
    /*
  if (tract_scalars_name!=NULL)
    {
      free_tract_colour=false;
      colour_by_direction=false;
    }
    */
  


  vtkLookupTable *colour_lookup_table=vtkLookupTable::New();
  colour_lookup_table->SetTableRange(0,number_of_tract_sets+1);
  colour_lookup_table->SetHueRange(0,1);
  colour_lookup_table->SetSaturationRange(1, 1);
  colour_lookup_table->SetValueRange(1, 1);
  colour_lookup_table->Build();

  
  for (i=0;i<number_of_tract_sets;i++)

    {
      fibre_tract_plotter_array[i]= new FibreTractPlotter(tract_set_filename_array[i]);
      if (tract_scalars_name!=NULL)
	{
	  
	  fibre_tract_plotter_array[i]->ColourByScalarFile(tract_scalars_name, scalar_low_thresh, scalar_high_thresh);
	  

	}
      if ( !colour_by_direction)
      
	{
	  if (free_tract_colour)
	    {
	      colour_lookup_table->GetColor((double)i,tract_colour);/*VTK5 expects doubles ilana*/
	      fibre_tract_plotter_array[i]->SetColour(tract_colour[0],tract_colour[1],tract_colour[2]);
	    }
	  //not sure:
	  //display_window->AddInteractiveTracts(fibre_tract_plotter_array[i]->CreateActor());
	  else
	    {
	      fibre_tract_plotter_array[i]->SetColour(tract_colour[3*i+0],tract_colour[3*i+1],tract_colour[3*i+2]);
	    }
	}
      else if (free_tract_colour)
	{
	  colour_lookup_table->GetColor((double)i,tract_colour);
	  //tract_colour[0]=1;
	  //tract_colour[1]=0;
	  //tract_colour[2]=0;
	  fibre_tract_plotter_array[i]->SetColour(tract_colour[0],tract_colour[1],tract_colour[2]);
	}
      else if (colour_by_direction && !map_tract_point_scalars)
	{
	 
	  fibre_tract_plotter_array[i]->ColourByDirection();
	}
      //fibre_tract_plotter_array[i]->SetColour(tract_colour[0],tract_colour[1],tract_colour[2]);
      if (plot_tubes)
	{
	  fibre_tract_plotter_array[i]->SetTubesOn(0.3);
	}
  
      if (map_tract_point_scalars)
	{
	  
	  fibre_tract_plotter_array[i]->MapPointScalarsOn();
	  //temporarily disabled interaction: using color.  change in FibreTractPlotter and here to go to opacity-windowing
	  //cout << "Opacity interaction disabled" << endl;



	  //display_window->AddInteractiveTracts(fibre_tract_plotter_array[i]->CreateActor());
	  
	  
	  display_window->AddObject(fibre_tract_plotter_array[i]->CreateActor());
	  
	}
      else if  (tract_scalars_name!=NULL)
	{
	  //HERE InteractorStyleFibreTract() (via display_window->AddInteractiveTracts) has to somehow know it is being invoked for windowing of colourbar, not alpha. set a case.  make an enum.  change above. could just have it be a different mouse click and let the user use the one that makes sense. no error msg..
	  //display_window->AddObject2D(fibre_tract_plotter_array[i]->ColourBar());
	  display_window->AddInteractiveTracts(fibre_tract_plotter_array[i]->CreateActor(),fibre_tract_plotter_array[i]->ColourBar());

	  //for no intereaction:
	  //display_window->AddObject(fibre_tract_plotter_array[i]->CreateActor());
	  //display_window->AddObject2D(fibre_tract_plotter_array[i]->ColourBar());

	}
      else
	{
	  
	  display_window->AddObject(fibre_tract_plotter_array[i]->CreateActor());
	  
	}


  }

  //DO: loop though all tract sets and free(fibre_tract_plotter[i]): not sure if can do this until the end, in which case need this outside...because interactor is part of the FibreTractPlotter...

  if (free_tract_colour)
    {
      free(tract_colour);
    }

}

void PlotROIs(char **roi_filename_array, int number_of_rois, DisplayWindow *display_window, float *low_thresh, float *high_thresh, double *roi_colour, double *roi_opacities, bool interactive_level_set) //ilana vtk5 SetOpacity takes double
{
  int i;
  LevelSetPlotter **level_set_plotter_array;
  level_set_plotter_array=(LevelSetPlotter **) malloc(number_of_rois*sizeof(LevelSetPlotter *));
  bool use_colour_table=false;
  vtkLookupTable *colour_lookup_table;
  double *roi_colour_double=(double *) malloc(3*sizeof(double));

  if (roi_colour==NULL)
    {
      use_colour_table=true;
      colour_lookup_table=vtkLookupTable::New();
      colour_lookup_table->SetTableRange(0,number_of_rois+1);
      colour_lookup_table->SetHueRange(0,1);
      colour_lookup_table->SetSaturationRange(1, 1);
      colour_lookup_table->SetValueRange(1, 1);
      colour_lookup_table->Build();
      roi_colour=(double *) malloc(3*sizeof(double));/*VTK5 expects doubles ilana*/

    }

  for (i=0;i<number_of_rois;i++)
    {
      level_set_plotter_array[i]=new LevelSetPlotter(roi_filename_array[i],low_thresh[i],high_thresh[i]);
      if (use_colour_table)
	{
	  //5.0
	  colour_lookup_table->GetColor((double)i,roi_colour);/*VTK5 expects doubles ilana*/

	  //test 4.4:
	  //roi_colour_double=colour_lookup_table->GetColor((double) i);
	  //roi_colour[0]=(float) roi_colour_double[0];
	  //roi_colour[1]=(float) roi_colour_double[1];
	  //roi_colour[2]=(float) roi_colour_double[2];
	}

      
      /*
      if (i==1)
	{
	  //red
	  level_set_plotter_array[i]->SetColour(3.0,0,0.0);
	}

      if (i==1)
	{
	  //green
	  level_set_plotter_array[i]->SetColour(0,3.0,0.0);
	}

      if (i==2)
	{
	  //blue
	  level_set_plotter_array[i]->SetColour(0,0,3.0);
	}

      else if (i==3)
	{
	  //turquoisish
	  level_set_plotter_array[i]->SetColour(3.0*0.1943,3.0*0.3951,3.0*0.4106);

	}

	else if (i==4)
	{
	  //yellow with a hint of green
	  level_set_plotter_array[i]->SetColour(1.4,1.6,0.0);
	}
      else if (i==5)
	{
	  //orangish
	  level_set_plotter_array[i]->SetColour(3.0*0.5021,3.0*0.3080,3.0*0.1899);
	}
      */

	/*if (i==1)
	{
	  //yellow
	  level_set_plotter_array[i]->SetColour(1.5,1.5,0.0);
	}
	else
	  {
	  level_set_plotter_array[i]->SetColour(1.5*roi_colour[0],1.5*roi_colour[1],1.5*roi_colour[2]);
	  }
	*/ //IRL added to be able to specify individual ROI colour
      if (use_colour_table)
      {
    	  level_set_plotter_array[i]->SetColour(1.5*roi_colour[0],1.5*roi_colour[1],1.5*roi_colour[2]);
      }
      else
      {
       level_set_plotter_array[i]->SetColour(1.5*roi_colour[i*3],1.5*roi_colour[i*3+1],1.5*roi_colour[i*3+2]);
      }
      //}

      level_set_plotter_array[i]->SetOpacity(roi_opacities[i]);

      if (!interactive_level_set)
	{
	  display_window->AddObject(level_set_plotter_array[i]->CreateActor());
	}
      else
	{
	  display_window->AddInteractiveLevelSet(level_set_plotter_array[i]);

	}


    }


  free(roi_colour);
  free(roi_colour_double);
}


void PlotVectors(char *file1, char *file2, char *file3, char *mask, char *scale_name,DisplayWindow *display_window,char *ODF_maxima_filename, float sigma_theta, char *sigma_theta_name, bool prob, Label geom, bool mostlikely, bool greyscale)

{

  

  //sigma_theta is in degrees  

  bool ODF;
  bool scale;
  int x,y,z,i,j;
  vtkPolyData *polydata=vtkPolyData::New();
  vtkHedgeHog *hedgehog=vtkHedgeHog::New();
  vtkImagePlaneWidget *image_plane_widget=vtkImagePlaneWidget::New();
  vtkPoints *points=vtkPoints::New();
  vtkFloatArray *scalars=vtkFloatArray::New(); 
  vtkFloatArray *vectors=vtkFloatArray::New();
  ImageData *image_data1;
  ImageData *image_data2;
  ImageData *image_data3;
  ImageData *image_data4;
  ImageData *image_data5;
  ImageData *image_data6;
  ImageData *mask_image_data;
  ImageData *scale_image_data;
  ImageData *mostlikelygeom;

  int counter;
  double v[3]; //JC changed this to double; check behaviour
  Directions *directions;
  POINT3D point;
  VECTOR3D vector;
  VECTOR3D vector_tmp;
  VECTOR3D vector0,vector1,vector2;

  ImageData *tmp;
  Directions *dirs;
  SIGMA_THETA sigma_theta_vox;

  float singleprob,doubleprob,tripleprob,fanprob;

  vtkConeSource *cone_source1;
  vtkConeSource *cone_source2;
  vtkPolyDataMapper *mapper;
  vtkActor *c_actor;

  bool therewasanonfan=false;


  vtkLookupTable *lut= vtkLookupTable::New();
  lut->SetRange(0.0, 1.0);
  lut->Build();   

  
  float nullv[3];
  nullv[0]=0;
  nullv[1]=0;
  nullv[2]=0;

  float minstep;


  //test
  float ODF_mag;
  int max_index;

  if (file2==NULL)
    {
      ODF=true;
    }
  else 
    {
      ODF=false;
    }

  if (scale_name==NULL)
    {
      scale=false;
    }
  else 
    { 
      scale=true;
    }



  if (prob)
    {
      //create full filenames and init ImageDatas for six ImageDatas:


	for (i=0; i<6; i++)
	{
	  char inttostring[5];
	  char OutFile[200];
		memset(OutFile, '\0', 200);
		strcpy(OutFile, file1);
		sprintf(inttostring,"%i",i);
		strcat(OutFile, inttostring);
		strcat(OutFile, ".mnc");
		switch (i)
		  {
		  case 0:
		    image_data1=new ImageData(OutFile);
		    break;
		  case 1:
		    image_data2=new ImageData(OutFile);
		    break;
		  case 2:
		    image_data3=new ImageData(OutFile);
		    break;
		  case 3:
		    image_data4=new ImageData(OutFile);
		    break;
		  case 4:
		    image_data5=new ImageData(OutFile);
		    break;
		  case 5:
		    image_data6=new ImageData(OutFile);
		    break;
		  }
	}

      
    }
  else //!prob
    {
    image_data1=new ImageData(file1);
    }

  mask_image_data=new ImageData(mask);


  if (prob)//we set the maxima and sigma_theta in a new dummy ODF image_data1 
    
    {

  
      //temporary routine to set the correct values in a standard ProbLabels ImageData from new format size files, and to set a volume of most likely labels if that is what was requested
      //note the loop structure 2x is not necessary
      if (mostlikely)
	{
	  mostlikelygeom=new ImageData(Image3D_vanilla,image_data1->GetXSize(),image_data1->GetYSize(),image_data1->GetZSize(),image_data1->GetXStart(),image_data1->GetYStart(),image_data1->GetZStart(),image_data1->GetXStep(),image_data1->GetYStep(),image_data1->GetZStep(),image_data1->GetXDirectionCosine(),image_data1->GetXDirectionCosine(),image_data1->GetXDirectionCosine());
	  //mostlikelygeom=new ImageData(ODFdata,image_data1->GetImageGeom());
	  tmp=new ImageData(ODFdata,image_data1->GetImageGeom());
	  for (x=0;x<image_data1->GetXSize();x++)
	    {
	      for (y=0;y<image_data1->GetYSize();y++)
		{
		  for (z=0;z<image_data1->GetZSize();z++)
		    {
		      if (mask_image_data->GetValue(x,y,z)>FLT_EPSILON)
			{
			  
			  //we transfer the values:
		
			  singleprob=image_data1->GetValue(x,y,z,0);
			  //cout << "singleprob " << singleprob << endl;
			  tmp->SetValue(x,y,z,1,image_data1->GetValue(x,y,z,2));//angle
			  tmp->SetValue(x,y,z,2,image_data1->GetValue(x,y,z,3));//mean|
			  tmp->SetValue(x,y,z,3,image_data1->GetValue(x,y,z,4));
			  tmp->SetValue(x,y,z,4,image_data1->GetValue(x,y,z,5));

			  doubleprob=image_data2->GetValue(x,y,z,0);
			  //cout << "doubleprob " << doubleprob << endl;
			  tmp->SetValue(x,y,z,6,image_data2->GetValue(x,y,z,2));
			  tmp->SetValue(x,y,z,7,image_data2->GetValue(x,y,z,3));
			  tmp->SetValue(x,y,z,8,image_data2->GetValue(x,y,z,4));
			  tmp->SetValue(x,y,z,9,image_data2->GetValue(x,y,z,5));

			  tmp->SetValue(x,y,z,11,image_data3->GetValue(x,y,z,2));
			  tmp->SetValue(x,y,z,12,image_data3->GetValue(x,y,z,3));
			  tmp->SetValue(x,y,z,13,image_data3->GetValue(x,y,z,4));
			  tmp->SetValue(x,y,z,14,image_data3->GetValue(x,y,z,5));

			  tripleprob=image_data4->GetValue(x,y,z,0);
			  //cout << "tripleprob " << tripleprob << endl;

			  tmp->SetValue(x,y,z,16,image_data4->GetValue(x,y,z,2));
			  tmp->SetValue(x,y,z,17,image_data4->GetValue(x,y,z,3));
			  tmp->SetValue(x,y,z,18,image_data4->GetValue(x,y,z,4));
			  tmp->SetValue(x,y,z,19,image_data4->GetValue(x,y,z,5));

			  tmp->SetValue(x,y,z,21,image_data5->GetValue(x,y,z,2));
			  tmp->SetValue(x,y,z,22,image_data5->GetValue(x,y,z,3));
			  tmp->SetValue(x,y,z,23,image_data5->GetValue(x,y,z,4));
			  tmp->SetValue(x,y,z,24,image_data5->GetValue(x,y,z,5));

			  tmp->SetValue(x,y,z,26,image_data6->GetValue(x,y,z,2));
			  tmp->SetValue(x,y,z,27,image_data6->GetValue(x,y,z,3));
			  tmp->SetValue(x,y,z,28,image_data6->GetValue(x,y,z,4));
			  tmp->SetValue(x,y,z,29,image_data6->GetValue(x,y,z,5));







			  //DO: redo the rest when figure out the fanning input

			  if (image_data1->GetNSize()>36)
			    {

			  fanprob=image_data1->GetValue(x,y,z,37);
			  //cout << "fanprob " << fanprob << endl;
			  tmp->SetValue(x,y,z,31,image_data1->GetValue(x,y,z,38));
			  tmp->SetValue(x,y,z,32,image_data1->GetValue(x,y,z,39));
			  tmp->SetValue(x,y,z,33,image_data1->GetValue(x,y,z,40));
			  tmp->SetValue(x,y,z,34,image_data1->GetValue(x,y,z,41));

			  tmp->SetValue(x,y,z,36,image_data1->GetValue(x,y,z,44));
			  tmp->SetValue(x,y,z,37,image_data1->GetValue(x,y,z,45));
			  tmp->SetValue(x,y,z,38,image_data1->GetValue(x,y,z,46));
			  tmp->SetValue(x,y,z,39,image_data1->GetValue(x,y,z,47));

			  tmp->SetValue(x,y,z,41,image_data1->GetValue(x,y,z,50));
			  tmp->SetValue(x,y,z,42,image_data1->GetValue(x,y,z,51));
			  tmp->SetValue(x,y,z,43,image_data1->GetValue(x,y,z,52));
			  tmp->SetValue(x,y,z,44,image_data1->GetValue(x,y,z,53));
			    }
			  
			  //we set geom:

			  if (singleprob>=doubleprob && singleprob>=tripleprob && singleprob>=fanprob)
			    {
			      geom=single_curve;
			      //cout << "single_curve\n";
			    }
			  else if (doubleprob>=singleprob && doubleprob>=tripleprob && doubleprob>=fanprob)
			    {
			      geom=doublecross;
			      //cout << "doublecross\n";
			    }
			  else if (tripleprob>=singleprob && tripleprob>=doubleprob && tripleprob>=fanprob)
			    {
			      geom=triplecross;
			      //cout << "triplecross\n";
			    }
			  else
			    {
			      geom=fan;
			      //cout << "fan\n";
			    }


			  mostlikelygeom->SetValue(x,y,z,(float) geom);

			}
		    }
		}
	    }
	  image_data1->Release();
	  image_data1=tmp;

	  
	}//mostlikely


      tmp=new ImageData(ODFdata,image_data1->GetImageGeom());
      for (x=0;x<image_data1->GetXSize();x++)
	{
	  for (y=0;y<image_data1->GetYSize();y++)
	    {
	      for (z=0;z<image_data1->GetZSize();z++)
		{
		  if (mask_image_data->GetValue(x,y,z)>FLT_EPSILON)
		    {

		      if (mostlikely)
			{
			  geom=(Label) mostlikelygeom->GetValue(x,y,z);
			}

		      //cout << "geom " << geom << endl;

		      switch (geom)
			{
			case single_curve:
			  //cout << "single_curve\n";
			  dirs=new Directions(1);
			  dirs->SetXComponent(0,image_data1->GetValue(x,y,z,2));
			  dirs->SetYComponent(0,image_data1->GetValue(x,y,z,3));
			  dirs->SetZComponent(0,image_data1->GetValue(x,y,z,4));
			  sigma_theta_vox.sigma_theta=(float *) malloc(1*sizeof(float));
			  sigma_theta_vox.sigma_theta[0]=image_data1->GetValue(x,y,z,1);
			  sigma_theta_vox.number=1;
			  sigma_theta_vox.calculated=true;

			  break;
			case doublecross:
			  //cout << "doublecross\n";
			  dirs=new Directions(2);
			  dirs->SetXComponent(0,image_data1->GetValue(x,y,z,7));
			  dirs->SetYComponent(0,image_data1->GetValue(x,y,z,8));
			  dirs->SetZComponent(0,image_data1->GetValue(x,y,z,9));
			  dirs->SetXComponent(1,image_data1->GetValue(x,y,z,12));
			  dirs->SetYComponent(1,image_data1->GetValue(x,y,z,13));
			  dirs->SetZComponent(1,image_data1->GetValue(x,y,z,14));
			  sigma_theta_vox.sigma_theta=(float *) malloc(2*sizeof(float));
			  sigma_theta_vox.sigma_theta[0]=image_data1->GetValue(x,y,z,6);
			  sigma_theta_vox.sigma_theta[1]=image_data1->GetValue(x,y,z,11);
			  sigma_theta_vox.number=2;
			  sigma_theta_vox.calculated=true;
			  break;
			case triplecross:
			  //cout << "triplecross\n";
			  dirs=new Directions(3);
			  dirs->SetXComponent(0,image_data1->GetValue(x,y,z,17));
			  dirs->SetYComponent(0,image_data1->GetValue(x,y,z,18));
			  dirs->SetZComponent(0,image_data1->GetValue(x,y,z,19));
			  dirs->SetXComponent(1,image_data1->GetValue(x,y,z,22));
			  dirs->SetYComponent(1,image_data1->GetValue(x,y,z,23));
			  dirs->SetZComponent(1,image_data1->GetValue(x,y,z,24));
			  dirs->SetXComponent(2,image_data1->GetValue(x,y,z,27));
			  dirs->SetYComponent(2,image_data1->GetValue(x,y,z,28));
			  dirs->SetZComponent(2,image_data1->GetValue(x,y,z,29));
			  sigma_theta_vox.sigma_theta=(float *) malloc(3*sizeof(float));
			  sigma_theta_vox.sigma_theta[0]=image_data1->GetValue(x,y,z,16);
			  sigma_theta_vox.sigma_theta[1]=image_data1->GetValue(x,y,z,21);
			  sigma_theta_vox.sigma_theta[2]=image_data1->GetValue(x,y,z,26);
			  sigma_theta_vox.number=3;
			  sigma_theta_vox.calculated=true;
			  break;
			case fan:
			  //cout << "fan\n";
			  dirs=new Directions(3);
			  dirs->SetXComponent(0,image_data1->GetValue(x,y,z,32));
			  dirs->SetYComponent(0,image_data1->GetValue(x,y,z,33));
			  dirs->SetZComponent(0,image_data1->GetValue(x,y,z,34));
			  dirs->SetXComponent(1,image_data1->GetValue(x,y,z,37));
			  dirs->SetYComponent(1,image_data1->GetValue(x,y,z,38));
			  dirs->SetZComponent(1,image_data1->GetValue(x,y,z,39));
			  dirs->SetXComponent(2,image_data1->GetValue(x,y,z,42));
			  dirs->SetYComponent(2,image_data1->GetValue(x,y,z,43));
			  dirs->SetZComponent(2,image_data1->GetValue(x,y,z,44));
			  sigma_theta_vox.sigma_theta=(float *) malloc(3*sizeof(float));
			  sigma_theta_vox.sigma_theta[0]=image_data1->GetValue(x,y,z,31);
			  sigma_theta_vox.sigma_theta[1]=image_data1->GetValue(x,y,z,36);
			  sigma_theta_vox.sigma_theta[2]=image_data1->GetValue(x,y,z,41);
			  sigma_theta_vox.number=3;
			}
		      

		      tmp->SetODFMaxima(dirs,x,y,z);
		      tmp->SetSigmaTheta(sigma_theta_vox,x,y,z);
		    }
		  
		}
	    }
	}     
      image_data1->Release();
      image_data2->Release();
      image_data3->Release();
      image_data4->Release();
      image_data5->Release();
      image_data6->Release();



      image_data1=tmp;
 
    }//prob



  else if (ODF)
    {
      if (ODF_maxima_filename!=NULL)
	{
	  image_data1->SetODFMaxima(ODF_maxima_filename);  
	  //temp:
	  //for prob., set maxima from most likely config:


	}
      if (sigma_theta_name!=NULL)
	{
	  //then you can't have a const input, so just in case:
	  sigma_theta=0;
	  image_data1->SetSigmaTheta(sigma_theta_name); 

	  //temp:
	  //for prob., set sigma-theta from most likely config:


	}

    }

  else //will be e1
    {
      image_data2=new ImageData(file2);
      image_data3=new ImageData(file3);
    }
  
  if (scale)
    {
      scale_image_data=new ImageData(scale_name);
    }


  vectors->SetNumberOfComponents(3);

 
  float height_scale=0.4; //0.5 to align with vector, 0.4 to stick out

 
  counter=0;
  for (x=0;x<image_data1->GetXSize();x++)
    {
      for (y=0;y<image_data1->GetYSize();y++)
	{
	  for (z=0;z<image_data1->GetZSize();z++)
	    {
	      if (mask_image_data->GetValue(x,y,z)>FLT_EPSILON)
		{
		  point=VoxelToWorld(x,y,z,image_data1->GetImageGeom());

		    if (mostlikely)
			{
			  geom=(Label) mostlikelygeom->GetValue(x,y,z);
			}  
		  
		  if (!ODF)
		    {
		      therewasanonfan=true;//true if using e1 vectors
		      points->InsertPoint(counter,point.x,point.y,point.z);
		      
		      //note that the vector is already in world space, and we will plot it in world space.
		      vector.x=image_data1->GetValue(x,y,z);
		      vector.y=image_data2->GetValue(x,y,z);
		      vector.z=image_data3->GetValue(x,y,z);

		      minstep=float_min3(float_abs(image_data1->GetXStep()), float_abs(image_data1->GetYStep()), float_abs(image_data1->GetZStep()));
		      
		      v[0]=(double) -0.5*minstep*vector.x;
		      v[1]=(double) -0.5*minstep*vector.y;
		      v[2]=(double) -0.5*minstep*vector.z;

		      
		      if (scale)
			{
			  v[0]=(double) scale_image_data->GetValue(x,y,z)*v[0];
			  v[1]=(double) scale_image_data->GetValue(x,y,z)*v[1];
			  v[2]=(double) scale_image_data->GetValue(x,y,z)*v[2];
			}
		      vectors->InsertTuple(counter,v);
		      scalars->InsertNextValue(ScalarForRGB(v,lut));
		      counter++;
		      
		      points->InsertPoint(counter,point.x,point.y,point.z);
		      
		      //other way:
		      
		      
		      v[0]=(double) 0.5*minstep*vector.x;
		      v[1]=(double) 0.5*minstep*vector.y;
		      v[2]=(double) 0.5*minstep*vector.z;
		      if (scale)
			{
			  v[0]=(double) scale_image_data->GetValue(x,y,z)*v[0];
			  v[1]=(double) scale_image_data->GetValue(x,y,z)*v[1];
			  v[2]=(double) scale_image_data->GetValue(x,y,z)*v[2];
			}
		      vectors->InsertTuple(counter,v);
		      scalars->InsertNextValue(ScalarForRGB(v,lut));
		      counter++;
		    }
		  else //ODF
		    {
		      //test:
		      //int *p=(int *) malloc(3000*sizeof(p));
		      //free(p);
		      //cout << geom << endl;
		      //memory-related problem here:

		      directions=image_data1->GetODFMaxima(x,y,z);
		      
		      if (geom==fan)
			{
			  //in world space already:
			  vector0=directions->GetDirection(0);
			  vector2=directions->GetDirection(2);

			  //if dirs not in world space:
			  //vector0=VoxelToWorld(directions->GetDirection(0),image_data1->GetImageGeom());
			  //vector2=VoxelToWorld(directions->GetDirection(2),image_data1->GetImageGeom());


			  for (j=0;j<20;j++)//j loop for fan
			  {
			  

			    points->InsertPoint(counter,point.x,point.y,point.z);
			 
			  //-1 because cones face the other way
			  vector.x=-1.0*((20.0-j)/20.0*vector0.x+j/20.0*vector2.x);
			  vector.y=-1.0*((20.0-j)/20.0*vector0.y+j/20.0*vector2.y);
			  vector.z=-1.0*((20.0-j)/20.0*vector0.z+j/20.0*vector2.z);
			  vector=Normalize(vector);
			  

			  minstep=float_min3(float_abs(image_data1->GetXStep()), float_abs(image_data1->GetYStep()), float_abs(image_data1->GetZStep()));

			  v[0]=(double) 0.5*minstep*vector.x;
			  v[1]=(double) 0.5*minstep*vector.y;
			  v[2]=(double) 0.5*minstep*vector.z;
			  if (scale)
			    {
			      v[0]=(double) scale_image_data->GetValue(x,y,z)*v[0];
			      v[1]=(double) scale_image_data->GetValue(x,y,z)*v[1];
			      v[2]=(double) scale_image_data->GetValue(x,y,z)*v[2];
			    }



			  //cone:

			  if (image_data1->SigmaThetaAvailable())
			    {
			      if (image_data1->GetSigmaTheta(x,y,z).calculated)
				{
				  
				  sigma_theta=(20.0-j)/20.0*180/(2*M_PI)*image_data1->GetSigmaTheta(x,y,z).sigma_theta[0]+j/20.0*180/(2*M_PI)*image_data1->GetSigmaTheta(x,y,z).sigma_theta[2];
				  
				  //sigma_theta=180/(2*M_PI)*image_data1->GetSigmaTheta(x,y,z).sigma_theta[i]; 
				}
			      else
				{
				  sigma_theta=0;
				}
			    }

			    

			  if (sigma_theta>FLT_EPSILON)
			    {
			      //cout << sigma_theta << endl;
			      cone_source1=vtkConeSource::New();
			      //cone_source1->SetHeight(float_abs(height_scale*image_data1->GetXStep()));
			      cone_source1->SetHeight(float_abs(height_scale*minstep));
			      cone_source1->SetAngle(sigma_theta);
			      //cone_source1->CappingOff();
			      cone_source1->SetDirection(v);
			      //set height half the vector height:
			      //cone_source1->SetHeight(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/2.0);

			      cone_source1->SetResolution(20);
			      cone_source1->SetCenter(point.x-0.5*cone_source1->GetHeight()*v[0],point.y-0.5*cone_source1->GetHeight()*v[1],point.z-0.5*cone_source1->GetHeight()*v[2]);
			      cone_source1->CappingOn();    
			      mapper=vtkPolyDataMapper::New();
			      mapper->SetInput(cone_source1->GetOutput());

			      c_actor=vtkActor::New();

			      c_actor->SetMapper(mapper);
			      if (!greyscale)
				{
			      c_actor->GetProperty()->SetColor(double_abs(v[0]),double_abs(v[1]),double_abs(v[2]));
				}


			      //display_window->AddObject(mapper);
			      display_window->AddObject(c_actor);
			      mapper->Delete();
			      cone_source1->Delete();
			    }


			  v[0]=-1.0*v[0];
			  v[1]=-1.0*v[1];
			  v[2]=-1.0*v[2];
			  //vectors->InsertTuple(counter,v);
			  //don't show anything:

			  vectors->InsertTuple(counter,nullv);
			  scalars->InsertNextValue(ScalarForRGB(v,lut));

			  counter++;
			  }//j loop
			}//fan
		    
		      if (geom!=fan)
			    {
			      therewasanonfan=true;
		      for (i=0;i<directions->GetNumberOfDirections();i++)
			{

			  points->InsertPoint(counter,point.x,point.y,point.z);
			
			  //in world space already:
			  vector=directions->GetDirection(i);
			  //vector=VoxelToWorld(directions->GetDirection(i),image_data1->GetImageGeom());


			  minstep=float_min3(float_abs(image_data1->GetXStep()), float_abs(image_data1->GetYStep()), float_abs(image_data1->GetZStep()));

			  v[0]=(double) 0.5*minstep*vector.x;
			  v[1]=(double) 0.5*minstep*vector.y;
			  v[2]=(double) 0.5*minstep*vector.z;
			  if (scale)
			    {
			      v[0]=(double) scale_image_data->GetValue(x,y,z)*v[0];
			      v[1]=(double) scale_image_data->GetValue(x,y,z)*v[1];
			      v[2]=(double) scale_image_data->GetValue(x,y,z)*v[2];
			    }



			  //cone:

			  if (image_data1->SigmaThetaAvailable())
			    {
			      if (image_data1->GetSigmaTheta(x,y,z).calculated)
				{
				  
				  sigma_theta=180/(2*M_PI)*image_data1->GetSigmaTheta(x,y,z).sigma_theta[i]; 
				}
			      else
				{
				  sigma_theta=0;
				}
			    }

			    

			  if (sigma_theta>FLT_EPSILON)
			    {
			      cone_source1=vtkConeSource::New();
			      cone_source1->SetHeight(float_abs(height_scale*minstep));
			      cone_source1->SetAngle(sigma_theta);
			      cone_source1->CappingOff();
			      cone_source1->SetDirection(v);
			      //set height half the vector height:
			      //cone_source1->SetHeight(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/2.0);

			      cone_source1->SetResolution(20);
			      cone_source1->SetCenter(point.x-0.5*cone_source1->GetHeight()*v[0],point.y-0.5*cone_source1->GetHeight()*v[1],point.z-0.5*cone_source1->GetHeight()*v[2]);
			      //cone_source1->CappingOn();    
			      mapper=vtkPolyDataMapper::New();
			      mapper->SetInput(cone_source1->GetOutput());
			      c_actor=vtkActor::New();

			      c_actor->SetMapper(mapper);
			      if (!greyscale)
				{
			      c_actor->GetProperty()->SetColor(double_abs(v[0]),double_abs(v[1]),double_abs(v[2]));
				}
			      
			      display_window->AddObject(c_actor);


			      //display_window->AddObject(mapper);
			      mapper->Delete();
			      cone_source1->Delete();
			    }

			  //end cone

			  v[0]=-1.0*v[0];
			  v[1]=-1.0*v[1];
			  v[2]=-1.0*v[2];
			  vectors->InsertTuple(counter,v);
			  
			  scalars->InsertNextValue(ScalarForRGB(v,lut));
			  
			  counter++;
			 

		

			  //other direction:
			 
			  points->InsertPoint(counter,point.x,point.y,point.z);



			  v[0]=(double) -0.5*minstep*vector.x;
			  v[1]=(double) -0.5*minstep*vector.y;
			  v[2]=(double) -0.5*minstep*vector.z;
			  if (scale)
			    {
			      v[0]=(double) scale_image_data->GetValue(x,y,z)*v[0];
			      v[1]=(double) scale_image_data->GetValue(x,y,z)*v[1];
			      v[2]=(double) scale_image_data->GetValue(x,y,z)*v[2];
			    }
			  //vectors->InsertTuple(counter,v);

			  //cone:
			  if (sigma_theta>FLT_EPSILON)
			    {
			      cone_source2=vtkConeSource::New();
			      cone_source2->SetHeight(float_abs(height_scale*minstep));
			      cone_source2->SetAngle(sigma_theta);
			      cone_source2->CappingOff();
			      cone_source2->SetDirection(v);
			      //cone_source2->SetHeight(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])/2.0);
			      cone_source2->SetResolution(20);
			      cone_source2->SetCenter(point.x-0.5*cone_source2->GetHeight()*v[0],point.y-0.5*cone_source2->GetHeight()*v[1],point.z-0.5*cone_source2->GetHeight()*v[2]);
			     
			      
			      mapper=vtkPolyDataMapper::New();
			      mapper->SetInput(cone_source2->GetOutput());
			      c_actor=vtkActor::New();
			      c_actor->SetMapper(mapper);
			      if (!greyscale)
				{
			      c_actor->GetProperty()->SetColor(double_abs(v[0]),double_abs(v[1]),double_abs(v[2]));
				}
			      //display_window->AddObject(mapper);
			      display_window->AddObject(c_actor);



			      mapper->Delete();
			      cone_source2->Delete();
			    }
			  //end cone

			  v[0]=-1.0*v[0];
			  v[1]=-1.0*v[1];
			  v[2]=-1.0*v[2];
			  vectors->InsertTuple(counter,v);

			 
			  scalars->InsertNextValue(ScalarForRGB(v,lut));

			  counter++;
			}//i loop      
		    }//if not fan 
			 
		      //directions->Release();//didn't retain them so can't release

		      

		}
		}
	    }
	}
    }

 
  
  polydata->SetPoints(points);
  polydata->GetPointData()->SetVectors(vectors);

  if (!greyscale)
    {  
      polydata->GetPointData()->SetScalars(scalars);
    }

  hedgehog->SetInput(polydata);
  hedgehog->SetScaleFactor(1);
	 

 

 
  mapper=vtkPolyDataMapper::New();
  mapper->SetInput(hedgehog->GetOutput()); 

  if (!greyscale)
    {
      mapper->ScalarVisibilityOn();
      mapper->SetLookupTable(lut);
      mapper->SetColorModeToMapScalars();
      mapper->SetScalarModeToUsePointData();
    }


  vtkActor *actor=vtkActor::New();
  actor->SetMapper(mapper);
 
  
  if (therewasanonfan)
    {
      display_window->AddObject(actor);
    }
  actor->GetProperty()->SetLineWidth(3);

  if (greyscale) //if (want black vectors)
    {
      actor->GetProperty()->SetColor(0,0,0);
      
    }

  if (!ODF)
    {
      image_data2->Release();
      image_data3->Release();
    }
  //image_data1->Release(); //test: this should be uncommented, but was causing seg fault (from upstream errors, which need to be found...)
  mask_image_data->Release();
  points->Delete();
  vectors->Delete();
  polydata->Delete();
  hedgehog->Delete();

  if (scale)
    {
      scale_image_data->Release();
    }

  scalars->Delete();
  lut->Delete();

  if (mostlikely)
    {
      mostlikelygeom->Release();
    }
}





void Usage()
{

	//IL made some adjustments to the options list
	  //JC added movie option
  printf("Usage: 3Dvis [options]\n\n");
  printf("options:\n\n");
  printf("-3DImage <3D_image.mnc>\n");
  printf("-tracts <tracts_1.vtk>...<tracts_N.vtk>\n");
  printf("-tubes: plot tracts as tubes\n");
    printf("-map_tract_point_scalars: map scalars to opacity at each point on tracts, if scalars exist in tract VTK file (need to be written to tract file by mincFibreTrack)\n");
  printf("-tractcolour <red> <green> <blue>...<redN> <blueN> <greenN>: colour each tract systems with user-chosen colours\n");//still supported
  printf("-colour_by_system: automatically colour each tract system with one different colour: default is to colour by direction\n");
  printf("-scalars <scalars_filename.mnc>: colour tracts using this scalar minc file (takes precendent over other coluring instructions)\n");
  printf("-scalar_range <low high>\n");
  printf("-rois <roi_file_1.mnc> <low_threshold> <high_threshold> ... <roi_file_N.mnc> <low_threshold> <high_threshold> (threshold is inclusive)\n");
  printf("-roicolour <red> <green> <blue>...<redN> <blueN> <greenN>: specify the colour for each roi\n");
  printf("-roiopacity <opacity_1> ... <opacity_n> (as many opacities as rois: must use -rois option first)\n");
  printf("-vector <filename.mnc> [<filename2.mnc> <filename3.mnc>]: plot either e1 if three inputs (e1x, e1y, e1z) are given or maxima if one ODF file given\n");
  printf("-vectormask <mask.mnc>: needs to be used with 'vector' option above\n");
  printf("-vectorscale <scale.mnc>: scale vector lengths s by values in this file\n");
  printf("-greyvectors: plot vectors (and cones of uncertainty) in colour grey; default is colour coded by orientation\n");
  printf("-noresample: assume all inputs are sampled in the same space (default is to resample them)\n");
  printf("-ODF <ODF_filename.mnc>\n");
  printf("-ODFmask <mask.mnc> (required if use the -ODF option)\n");
  printf("-ODFdisplay: Stretch, SubtractMinimum, Physical (default is Physical)\n");
  printf("-ODF_maxima <maxima_filename> (optional: can input these if already calculated, but must use -noresample)\n");
  printf("-sigma_theta_const <sigma_theta>: plots cone subtending twice this angle (in degrees) around vectors (i.e., rotate angle sigma_theta around vectors)\n\n");
 printf("-sigma_theta <sigma_theta.txt>: input sigma_theta for each maximum.  This file must be paired with a -ODF_maxima input.  Plots cone subtending twice this angle (in degrees) around vectors (i.e., rotate angle sigma_theta around vectors)\n\n");
  printf("-ODF_opacity <opacity>\n");
  printf("-ODF_colour <red> <green> <blue>: if nothing given, RGB scheme will be used\n");
  printf("-ProbLabels <problabels>: input ProbLabel format files (for details, see ProbDeconvolve) and plot vectors and uncertainties based on the -label option (below)\n");
  printf("-label <which_label>: <which_label> is one of \"single\",\"double\",\"triple\", \"fan\", or \"mostlikely\".  \n");
  printf("-movie <movie_name>: create .mpg movie_name.mpg.  To do so, hold shift while rotate scene.\n");
  printf("\n\nthere must be only one -tracts, -tractcolour, -rois and -roicolour option\n");
  Error(" ");



}





int main(int argc, char *argv[])
{
  int i,j;

  bool greyvectors=false;

  char *image_filename;
  char *tract_set_filename_array[100]; //could be dynamic but isn't right now...
  int number_of_tract_sets;
  bool plot_tubes=false;


  bool map_tract_point_scalars=false;
  bool interactive_level_set=false;


  char *roi_filename_array[100];
  int number_of_rois;
  float low_thresh[100];
  float high_thresh[100];

  double *roicolour=NULL;
  double *tractcolour=NULL;
  double *roi_opacities; //ilana vtk5 roi_opacities takes double
  bool plot_3Dimage=false;
  bool plot_rois=false;
  bool plot_tracts=false;
  bool plot_odfs=false;
  bool plot_e1=false;
  bool plot_ODF_maxima=false;

  char *e1x_filename;
  char *e1y_filename;
  char *e1z_filename;
  char *ODF_vector_filename;
  char *vector_mask_name;
  char *vector_scale_name=NULL;

  char *ODF_input_name;
  char *ODF_mask_name;
  DISPLAY_MODE ODF_display_mode=Physical;//Stretch; //IRL change the default

  char *moviename=NULL;

  char *ODF_maxima_filename=NULL;
  char *sigma_theta_filename=NULL;

  float ODF_opacity_float;
  double ODF_opacity=1.0; //ilana vtk5 SetOpaciy takes double
  float ODF_colour[3]={-1,0,0}; //this is default and will result in RGB


  char *tempfilename_array[100];

  char templatefilename[200];

  char tempdirname[200];
  char command[300];
  bool templatefilename_set=true;

  bool resample=true;

  bool debian2_6=true;

  bool remove_zero_b=false; //isn't necesary: ODF Plotter checks already


  float sigma_theta=0;

  Label geom=none;
  ImageData *prob_labels;
  char *prob_labels_name=NULL;
  bool mostlikely=true;

  bool colour_by_direction=true;//for tracts; default is true; becomes false if colour_by_system is requested.
  bool colour_by_system=false;

  char *tract_scalars_name=NULL;
  float scalar_high_thresh=FLT_MAX;
  float scalar_low_thresh=FLT_MIN;

  // Check number of input file names
  if (argc < 3)
    {
      Usage();
    }



  // parse options: can be in any order

  for(i = 1; i < argc; i++)
    {
      if(!strcmp(argv[i], "-help"))
	{
	  Usage();
	}
      else if(!strcmp(argv[i], "-tubes"))
	{
	  plot_tubes = true;
	}
      else if(!strcmp(argv[i], "-map_tract_point_scalars"))
	{
	  map_tract_point_scalars = true;
	}
      else if(!strcmp(argv[i], "-noresample"))
	{
	  resample = false;
	}
      else if(!strcmp(argv[i], "-interactive_level_set"))
	{
	  interactive_level_set = true;
	}
      /*else if(!strcmp(argv[i], "-remove_zero_bvalues"))
	{
	  remove_zero_b = true;
	}*/
      else if (!strcmp(argv[i], "-3DImage"))
	{
	  if(++i >= argc) Usage();
	  image_filename = argv[i];
	  plot_3Dimage=true;

	}
      else if (!strcmp(argv[i], "-scalars"))
	{
	  if(++i >= argc) Usage();
	  tract_scalars_name = argv[i];
	}
      else if (!strcmp(argv[i], "-scalar_range"))
	{
	  if(++i >= argc) Usage();
	  sscanf(argv[i++], "%f", &scalar_low_thresh);
	  sscanf(argv[i], "%f", &scalar_high_thresh);
	  
	}

      else if (!strcmp(argv[i], "-movie"))
	{
	  if(++i >= argc) Usage();
	  moviename = argv[i];
	}
      else if (!strcmp(argv[i], "-ProbLabels"))
	{
	  if(++i >= argc) Usage();
	  prob_labels_name = argv[i];
	}
      else if (!strcmp(argv[i], "-label"))
	{
	  if(++i >= argc) Usage();
	  if (!strcmp(argv[i], "single"))
	    {
	      geom=single_curve;
	      mostlikely=false;
	    }
	  else if (!strcmp(argv[i], "double"))
	    {
	      geom=doublecross;
	      mostlikely=false;
	    }
	  else if (!strcmp(argv[i], "triple"))
	    {
	      geom=triplecross;
	      mostlikely=false;
	    }
	  else if (!strcmp(argv[i], "fan"))
	    {
	      geom=fan;
	      mostlikely=false;

	    }
	  else if (!strcmp(argv[i], "mostlikely"))//this is default
	    {
	      mostlikely=true;
	    }
	}
      else if (!strcmp(argv[i], "-tracts"))
	{
	  if(++i >= argc) Usage();
	  j=0;
	  while (argv[i][0] != '-')
	    {
	      tract_set_filename_array[j] = argv[i];
	      i++;
	      j++;
	      if (i>=argc)
		{
		  break;
		}
	    }
	  i--;

	  number_of_tract_sets=j;
	  plot_tracts=true;

	}
      else if (!strcmp(argv[i], "-rois"))
	{
	  if(++i >= argc) Usage();
	  j=0;
	  while (argv[i][0] != '-')
	    {
	      roi_filename_array[j] = argv[i++];
	      sscanf(argv[i++], "%f", &low_thresh[j]);
	      sscanf(argv[i++], "%f", &high_thresh[j]);
	      j++;
	      if (i>=argc)
		{
		  break;
		}
	    }
	  i--;

	  number_of_rois=j;
	  plot_rois=true;
	  //set default opacities:
	  roi_opacities=(double *) malloc(number_of_rois*sizeof(double)); //ilana vtk5 SetOpacity takes double
	  for (j=0;j<number_of_rois;j++)
	    {
	      roi_opacities[j]=0.5;
	    }

	}

      else if (!strcmp(argv[i], "-roicolour"))
	{
	  if(++i >= argc) Usage();
	  j=0;
	  roicolour=(double *) malloc(100*sizeof(double));
	  while (argv[i][0] != '-')
	    {
		  sscanf(argv[i++], "%lf", &roicolour[j]);
	      j++;
	      if (i>=argc)
		  {
		    break;
		  }
	    }
	  i--;

	}
      else if (!strcmp(argv[i], "-tractcolour"))
	{
	  if(++i >= argc) Usage();
	  j=0;
	  tractcolour=(double *) malloc(100*sizeof(double));
	  while (argv[i][0] != '-')
	    {
		  sscanf(argv[i++], "%lf", &tractcolour[j]);
	      j++;
	      if (i>=argc)
		  {
		    break;
		  }
	    }
	  i--;

	}

      else if (!strcmp(argv[i], "-roiopacity"))
	{
	  if(++i >= argc) Usage();

	  for (j=0;j<number_of_rois;j++)
	    {
	      sscanf(argv[i++], "%lf", &roi_opacities[j]); //ilana vtk5 SetOpacity takes double
	    }
	  i--;
	}
      else if (!strcmp(argv[i], "-ODF"))
	{
	  if(++i >= argc) Usage();
	  ODF_input_name=argv[i];
	  plot_odfs=true; //DO: check for ODFmask (should be required...)

	}
      else if (!strcmp(argv[i], "-ODF_opacity"))
	{
	  if(++i >= argc) Usage();
	  sscanf(argv[i], "%f", &ODF_opacity_float); //ilana vtk5 SetOpacity is double
	  ODF_opacity=(double) ODF_opacity_float;

	}
      else if (!strcmp(argv[i], "-ODF_colour"))
	{
	  if(++i >= argc) Usage();
	  sscanf(argv[i++], "%f", &ODF_colour[0]);
	  sscanf(argv[i++], "%f", &ODF_colour[1]);
	  sscanf(argv[i], "%f", &ODF_colour[2]);
	}
      else if (!strcmp(argv[i], "-ODFmask"))
	{
	  if(++i >= argc) Usage();
	  ODF_mask_name = argv[i];
	}
      else if (!strcmp(argv[i], "-sigma_theta_const"))
	{
	  if(++i >= argc) Usage();
	  sscanf(argv[i], "%f", &sigma_theta);
	}
      else if (!strcmp(argv[i], "-sigma_theta"))
	{
	  if(++i >= argc) Usage();
	  sigma_theta_filename = argv[i];
	}
      else if (!strcmp(argv[i], "-ODFdisplay"))
	{
	  if(++i >= argc) Usage();
	  if (!strcmp(argv[i], "SubtractMinimum"))
	    {
	      ODF_display_mode=SubtractMinimum;
	    }
	  if (!strcmp(argv[i], "Physical"))
	    {
	      ODF_display_mode=Physical;
	    }
	}
      else if (!strcmp(argv[i], "-ODF_maxima"))
	{
	  if(++i >= argc) Usage();
	  ODF_maxima_filename = argv[i];

	}
      else if (!strcmp(argv[i], "-vector"))
	{
	  if(++i >= argc) Usage();
	  //check how many args follow to determine whether e1 or maxima:
	  if (i+1>=argc)
	    {
	      ODF_vector_filename=argv[i];
	      plot_ODF_maxima=true;
	    }
	  else if (argv[i+1][0] == '-')
	    {
	      ODF_vector_filename=argv[i];
	      plot_ODF_maxima=true;
	    }
	  else //if there are 3 more (not checking...)
	    {
	      e1x_filename=argv[i];
	      i++;
	      e1y_filename=argv[i];
	      i++;
	      e1z_filename=argv[i];
	      plot_e1=true;
	    }

	  //DO: error check for vectormask

	}

      else if (!strcmp(argv[i], "-vectormask"))
	{
	  if(++i >= argc) Usage();
	  vector_mask_name = argv[i];
	}
      else if (!strcmp(argv[i], "-vectorscale"))
	{
	  if(++i >= argc) Usage();
	  vector_scale_name = argv[i];
	}
      else if (!strcmp(argv[i], "-greyvectors"))
	{
	  greyvectors=true;
	}
     else if (!strcmp(argv[i], "-colour_by_system"))
      	{
    	 colour_by_system=true;
      	}


      //an unknown option:
      else if(argv[i][0] == '-')
	{
	  Usage();
	}


    }







  //for file resampling:


  sprintf(tempdirname,"/tmp/tmp%i",(int) time(NULL));
  sprintf(command,"mkdir %s\n",tempdirname);
  system(command);



  //create DisplayWindow:
  DisplayWindow *display_window = new DisplayWindow();


  if (plot_3Dimage)
    {
      templatefilename_set=false;
      tempfilename_array[0]=(char *) malloc(200*sizeof(char));


      //here, we just give standard dir cos (because vtk imagedata likes it that way)


    if (debian2_6)

	{
	  tempfilename_array[1]=(char *) malloc(200*sizeof(char));

	  ResampleStandardDirCos(image_filename,tempfilename_array[1],tempdirname,true);

	  ResampleNegStep(tempfilename_array[1],tempfilename_array[0],tempdirname);
	}
      else

	{
	  ResampleStandardDirCos(image_filename,tempfilename_array[0],tempdirname,false);
	}



      // IRL want to have a hold on the ThreeOrthogonalPlanesPlotter to be able to delete it later (seg faults if done in the logical place, Plot3DImage)
      ImageData *image_data=new ImageData(tempfilename_array[0]); // TMP
      ThreeOrthogonalPlanesPlotter *three_planes_plotter= new ThreeOrthogonalPlanesPlotter(image_data); //IRL add this here
      //image_data->Release();  // CLAUDE
      Plot3DImage(tempfilename_array[0], display_window, three_planes_plotter );
      //Plot3DImage(tempfilename_array[0], display_window);


     delete three_planes_plotter;   // CLAUDE
     image_data->Release();  // CLAUDE // TMP



      free(tempfilename_array[0]);
      if (debian2_6)
	{
	  free(tempfilename_array[1]);
	}

    }
  if (plot_tracts)
    {

	 
      
	if (colour_by_system) //does this need  || tractcolour!= NULL ?
	    {
	      colour_by_direction=false;
	
	    }
	  
	  PlotTracts(tract_set_filename_array, number_of_tract_sets, display_window, plot_tubes, map_tract_point_scalars, tractcolour, colour_by_direction, tract_scalars_name, scalar_low_thresh, scalar_high_thresh);

    }
  if (plot_rois)
    {

      for (i=0;i<number_of_rois;i++)
	{
	  //templatefilename_set=false;
	  tempfilename_array[i]=(char *) malloc(200*sizeof(char));

	  //Resample(roi_filename_array[i],tempfilename_array[i],tempdirname,0,true,templatefilename,templatefilename_set);
	  //tempfilename_array[number_of_rois+i]=(char *) malloc(200*sizeof(char));

	  if (debian2_6) //appear to need dirs all neg (?) on new Debian systems.  No explanation for this...


	    {
	      tempfilename_array[number_of_rois+i]=(char *) malloc(200*sizeof(char));

	       ResampleStandardDirCos(roi_filename_array[i],tempfilename_array[number_of_rois+i],tempdirname,true);

               //ResampleStandardDirCos(roi_filename_array[i],tempfilename_array[i],tempdirname,true);

	      //sprintf(tempfilename_array[number_of_rois+i],"%s/%s",tempdirname,tempfilename_array[i]);


	      ResampleNegStep(tempfilename_array[number_of_rois+i],tempfilename_array[i],tempdirname);
	      //ResampleNegStep(roi_filename_array[i],tempfilename_array[i],tempdirname);
	    }
	  else

	    {
	      ResampleStandardDirCos(roi_filename_array[i],tempfilename_array[i],tempdirname,false);

	    }

	}



       PlotROIs(tempfilename_array, number_of_rois, display_window, low_thresh, high_thresh,roicolour,roi_opacities, interactive_level_set);
      //PlotROIs(roi_filename_array, number_of_rois, display_window, low_thresh, high_thresh,roicolour,roi_opacities, interactive_level_set);
      if (debian2_6)


	{
	  for (i=0;i<number_of_rois*2;i++)
	    {
	      free(tempfilename_array[i]);
	    }
	}
      else
	{
	  for (i=0;i<number_of_rois;i++)
	    {
	      free(tempfilename_array[i]);
	    }
	}
    }

  if (plot_odfs)
    {
      if (resample)
	{
	  templatefilename_set=false;
	  for (i=0;i<2;i++)
	    {
	      tempfilename_array[i]=(char *) malloc(200*sizeof(char));
	    }
	  Resample(ODF_input_name,tempfilename_array[0],tempdirname,0,true,templatefilename,templatefilename_set,true);
	  //strcpy(templatefilename,tempfilename_array[0]);
	  //templatefilename_set=true;

	  Resample(ODF_mask_name,tempfilename_array[1],tempdirname,0,true,templatefilename,templatefilename_set,true);


	  PlotODFs(tempfilename_array[0],tempfilename_array[1],ODF_display_mode, display_window,ODF_opacity,ODF_colour);
	  for (i=0;i<2;i++)
	    {
	      free(tempfilename_array[i]);
	    }
	}
      else
	{

	  PlotODFs(ODF_input_name,ODF_mask_name,ODF_display_mode, display_window,ODF_opacity,ODF_colour);
	}
    }

  if (plot_e1)
    {
      if (resample)
	{
	  templatefilename_set=false;
	  for (i=0;i<5;i++)
	    {
	      tempfilename_array[i]=(char *) malloc(200*sizeof(char));
	    }
	  Resample(e1x_filename,tempfilename_array[0],tempdirname,0,true,templatefilename,templatefilename_set,false);
	  Resample(e1y_filename,tempfilename_array[1],tempdirname,0,true,templatefilename,templatefilename_set,false);
	  Resample(e1z_filename,tempfilename_array[2],tempdirname,0,true,templatefilename,templatefilename_set,false);
	  Resample(vector_mask_name,tempfilename_array[3],tempdirname,0,true,templatefilename,templatefilename_set,false);
	  if (vector_scale_name!=NULL)
	    {
	      Resample(vector_scale_name,tempfilename_array[4],tempdirname,0,true,templatefilename,templatefilename_set,false);
	      PlotVectors(tempfilename_array[0],tempfilename_array[1],tempfilename_array[2],tempfilename_array[3],tempfilename_array[4],display_window,NULL,sigma_theta,NULL,false,geom,false,greyvectors);

	    }
	  else //Q why two calls?  took one out
	    {
	      PlotVectors(tempfilename_array[0],tempfilename_array[1],tempfilename_array[2],tempfilename_array[3],NULL,display_window,NULL,sigma_theta,NULL,false,geom,false,greyvectors);

	    }

	  for (i=0;i<5;i++)
	    {
	      free(tempfilename_array[i]);
	    }
	}
      else
	{
	  PlotVectors(e1x_filename,e1y_filename,e1z_filename,vector_mask_name,vector_scale_name, display_window,NULL,sigma_theta,NULL,false,geom,false,greyvectors);


	}
    }

  if (prob_labels_name!=NULL)
    {
      if (resample)
	{
	  Error("must call -noresample with this option");
	}


      PlotVectors(prob_labels_name,NULL,NULL,vector_mask_name, vector_scale_name, display_window, ODF_maxima_filename,sigma_theta, sigma_theta_filename,true,geom,mostlikely,greyvectors);

    }

 if (plot_ODF_maxima)
    {

      if (resample)
	{
	  templatefilename_set=false;
	  for (i=0;i<4;i++)
	    {
	      tempfilename_array[i]=(char *) malloc(200*sizeof(char));
	    }
	  Resample(ODF_vector_filename,tempfilename_array[0],tempdirname,0,true,templatefilename,templatefilename_set,true);
	  if (ODF_maxima_filename!=NULL)
	    {
	      Resample(ODF_maxima_filename,tempfilename_array[1],tempdirname,0,true,templatefilename,templatefilename_set,true);
	    }
	  else
	    {
	      tempfilename_array[1]=NULL;
	    }
	  Resample(vector_mask_name,tempfilename_array[2],tempdirname,0,true,templatefilename,templatefilename_set,true);//JC changed to nearest because ODF is so and often plot both - so same voxels plotted
	  if (vector_scale_name!=NULL)
	    {
	      Resample(vector_scale_name,tempfilename_array[3],tempdirname,0,true,templatefilename,templatefilename_set,true);
	      PlotVectors(tempfilename_array[0],NULL,NULL,tempfilename_array[2],tempfilename_array[3] , display_window,tempfilename_array[1],sigma_theta,NULL,false,geom,false,greyvectors);

	    }
	  else
	    {
	      PlotVectors(tempfilename_array[0],NULL,NULL,tempfilename_array[2],NULL, display_window,tempfilename_array[1] ,sigma_theta,NULL,false,geom,false,greyvectors);

	    }

	  for (i=0;i<4;i++)
	    {
	      free(tempfilename_array[i]);
	    }

	}
      else
	{
	  PlotVectors(ODF_vector_filename,NULL,NULL,vector_mask_name, vector_scale_name, display_window, ODF_maxima_filename,sigma_theta, sigma_theta_filename,false,geom,false,greyvectors);


	}

    }


  if (moviename!=NULL)
    {
      display_window->MoviePrep(moviename);
    }


  //plot:
  display_window->Display();




  delete(display_window);
  //delete (three_planes_plotter);   // CLAUDE

  //clean up temp directory:
  sprintf(command,"rm -f %s/tmp*\n",tempdirname);
  //cout << command << endl;


  system(command);

  sprintf(command,"rmdir %s\n",tempdirname);
  //cout << command << endl;

  system(command);

  return(0);



}

