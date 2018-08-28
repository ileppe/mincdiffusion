/************************************************************************

   File: ODFvisualize.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __jc_tools_h
#include "jc_tools.h"
#endif

#ifndef __otherlibs_include_h
#include "otherlibs_include.h"
#endif



void PlotODFs(char *ODF_filename, char *mask_filename, DISPLAY_MODE display_mode)
{

  int counter;
  ODFPlotter **ODFplotter_array;
  int x,y,z;
  gsl_vector_float *voxel_scale_factors;
  float minimum_scale_factor;
  int i;

  //could resample mask to match ODF (with mincresample): currently assume done

  //create ImageData for each input:

  ImageData *ODF_image_data = new ImageData(ODF_filename);
  ImageData *mask_image_data = new ImageData(mask_filename);

  //create Display window:

  DisplayWindow *display_window = new DisplayWindow;


  //allocate max size this array may have to be (wasteful: could use realloc):
  ODFplotter_array=(ODFPlotter **)malloc(ODF_image_data->GetXSize()*ODF_image_data->GetYSize()*ODF_image_data->GetZSize()*sizeof(ODFPlotter *));

  
  ODFPlotter *ODFplotter;
  

  //loop through all voxels of ODF: if mask>0.5; dynamically allocate a new ODFPlotter and prepare it.  keep track of their best scale factors in a gsl_vector_float array.  keep track of how many plotting voxels there are.

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

		  ODFplotter = new ODFPlotter(ODF_image_data);
		  ODFplotter->SetSmoothOn();
		  ODFplotter->SetVoxel(x,y,z);
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
		      
		    }
		  else 
		    {
		      ODFplotter->GetInterpolator()->CalculateLookupTable();
		      
		    }

	
		  ODFplotter->CalculatePlottingPoints();
		
		  if(ODFplotter->GetBestScaleFactorForVoxelFit()<minimum_scale_factor)
		    {
		      minimum_scale_factor=ODFplotter->GetBestScaleFactorForVoxelFit();
		    }
		  
		  ODFplotter_array[counter]=ODFplotter;
		  counter++;
		  
		}
	      
	    }
	  
	}
      
    }
  

  

  for (i=0;i<counter;i++)
    {
      ODFplotter_array[i]->SetScaleFactor(minimum_scale_factor);
      display_window->AddObject(ODFplotter_array[i]->CreateActor());
      
    }
  
  display_window->Display();
      
  //DO: delete all ODF Plotters in array and then:
  free(ODFplotter_array);
  

  
}

    
void Usage()
{
  Error("Usage: ODFvisualize [options] -mask <mask.mnc> <input.mnc>\n\noptions:\n\n-display: Stretch, SubtractMinimum (default is Physical)");
}




int main(int argc, char *argv[])
{
  char *input_name;
  char *mask_name;
  int i;
  DISPLAY_MODE display_mode;
  

  // Check number of input file names 
  if (argc < 4)
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
      else if (!strcmp(argv[i], "-mask")) 
	{
	  if(++i >= argc) Usage();
	  mask_name = argv[i];
	}
      else if (!strcmp(argv[i], "-display")) 
	{
	  if(++i >= argc) Usage();
	  if (!strcmp(argv[i], "SubtractMinimum"))
	    {
	      display_mode=SubtractMinimum;
	    }
	  if (!strcmp(argv[i], "Stretch"))
	    {
	      display_mode=Stretch;
	    }
	}


      //an unknown option:
      else if(argv[i][0] == '-') 
	{
	  Usage();
	}
      
      
      else 
	{
	  input_name=argv[i];
	}
      

    }
  
  //assumes mask file is resampled to match. but could check and resample it here if it doesn't


  PlotODFs(input_name,mask_name,display_mode);
  
  return(0);
  
  
}


