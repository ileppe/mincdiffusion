/************************************************************************

   File: Deconvolve.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2007 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/

#include <float.h>
#include "jc_tools.h"
#include "otherlibs_include.h"


/*

the FORECAST algorithm does not require the eigenvectors: currently, use defaults or input the eigenvalues and a mask.



currently outputs samples of fODF (FOD) at same sample points as input data.  If have super-resolution, should increase number of sample points.



currently always normalizes, to unit volume.

currently checks for nonzero b values (if asked to) on the whole volumes because that's how it's implemented in ImageData.  This loop is time-costly. (e.g., if response were input not as a ImageData, could do just on the one voxel...)




*/


void Deconvolve(char *input_name, char *mask_name, char *output_name, char *response_name, char *response_mask_name, char *history, bool clobber, bool use_e1, char *e1x_name, char *e1y_name, char *e1z_name,int order,char *lambda1_name,char *lambda2_name, bool def)
{

  DeconvolutionAlgorithm algorithm=FORECAST;//hardcoded here for now
  int x,y,z;
  ImageData *deconvolved;
  Deconvolver *deconvolver= new Deconvolver();
  deconvolver->SetOrder(order);
  deconvolver->SetDeconvolutionAlgorithm(algorithm);
  ImageData *image_data=new ImageData(input_name);
  ImageData *nonzero_b_image;
  ImageData *response;
  ImageData *response_nonzero_b;
  ImageData *response_voxel;
  ImageData *mask=new ImageData(mask_name);
  ImageData *e1x;
  ImageData *e1y;
  ImageData *e1z;
  ImageData *lambda1;
  ImageData *lambda2;
  ImageData *tmp_image;
  
  int i;

  /*if (remove_zero_b)
    {*/
      nonzero_b_image=image_data->CreateDWIODF();
      
   /* }
   else 
    {
      nonzero_b_image=image_data; //IRL should not get here, use full diffusion set as input
  
    }*/

      if (nonzero_b_image->GetNSize()<50)
	{
	  tmp_image=Create100DirInput(nonzero_b_image);
	  nonzero_b_image->Release();
	  nonzero_b_image=tmp_image;
	}

  deconvolver->SetInput(nonzero_b_image); 
		  

  if (use_e1)
    {
      /*if (remove_zero_b)
	{*/
	  response_nonzero_b=response->CreateDWIODF();
	/*}
      else 
	{
	  response_nonzero_b=response;
	}*/

      
	
      response=new ImageData(response_name);
      e1x=new ImageData(e1x_name);
      e1y=new ImageData(e1y_name);
      e1z=new ImageData(e1z_name);
      response_voxel=new ImageData(response_mask_name);
      deconvolver->SetResponseFunctionSingleVoxel(response_nonzero_b, response_voxel, e1x, e1y, e1z);
    }
  else if (algorithm==FORECAST)
    {
      //response=new ImageData(response_name);//IRL we don't need this anymore
      //IRL using user input for response_voxel
      if (!def) {
        lambda1=new ImageData(lambda1_name);
        lambda2=new ImageData(lambda2_name);
	//deconvolver->SetResponseFunction(response);//note that for this algorithm, response has b=0 images and they must be at the beginning
	deconvolver->SetResponseFunction(image_data);//IRL use the inputted diffusion set directly
	response_voxel=new ImageData(response_mask_name);
        deconvolver->SetLambdas(lambda1,lambda2,response_voxel);
      
      }else{
	deconvolver->SetResponseFunction(image_data);
	deconvolver->SetLambdasDefault();
	
      }
      
    }
  else //Tournier and use the unfitted DWIs to find direction
    {
      /*if (remove_zero_b)
	{*/
	  response_nonzero_b=response->CreateDWIODF();
	/*}
      else 
	{
	  response_nonzero_b=response;
	}*/
      response_voxel=new ImageData(response_mask_name);
      deconvolver->SetResponseFunctionSingleVoxel(response_nonzero_b, response_voxel);
    }

  //loop through and deconvolve:


  for (x=0;x<nonzero_b_image->GetXSize();x++)
    {
      for (y=0;y<nonzero_b_image->GetYSize();y++)
	{
	  for (z=0;z<nonzero_b_image->GetZSize();z++)
	    {
	      if (mask->GetValue(x,y,z)>FLT_EPSILON)
		{
		  //cout << "deconvolving\n";
		  deconvolver->Deconvolve(x,y,z);
		}
	    }
	}
    }

  deconvolved=deconvolver->GetOutput();


  //normalize, and zero the other voxels: normalize is currently making it harder to see, so not.  it is nowhere near normalizing: haven't looked into why at this point

  if (0) 
    {
  for (x=0;x<nonzero_b_image->GetXSize();x++)
    {
      for (y=0;y<nonzero_b_image->GetYSize();y++)
	{
	  for (z=0;z<nonzero_b_image->GetZSize();z++)
	    {
	      if (mask->GetValue(x,y,z)>FLT_EPSILON)
		{
		  deconvolved->NormalizeODF(x,y,z);
		}
	      else
		{
		  for (i=0;i<deconvolved->GetNSize();i++)
		    {
		      deconvolved->SetValue(x,y,z,i,0.0);
		    }
		}

	    }
	}
    }

    }//normalize and zero or not
  
  ///////////////////////////////////////////////////////////////////////

  //debug: have a look at the results: //mod from  PlotODFs 
  if (0)
    {
      DisplayWindow *display_window = new DisplayWindow();

      
  int counter;
  ODFPlotter **ODFplotter_array;
  int x,y,z;
  gsl_vector_float *voxel_scale_factors;
  float minimum_scale_factor;
  int i;
  float scale_factor;
  DISPLAY_MODE display_mode=Physical;
  double ODF_opacity=1.0;
  float colour[]={1.0,0.0,0.0};

  ODFplotter_array=(ODFPlotter **)malloc(deconvolved->GetXSize()*deconvolved->GetYSize()*deconvolved->GetZSize()*sizeof(ODFPlotter *));
  
  ODFPlotter *ODFplotter;
  

  //loop through all voxels of ODF: if mask>0.5; dynamically allocate a new ODFPlotter and prepare it.  keep track of their best scale factors in a gsl_vector_float array.  keep track of how many plotting voxels there are.




  counter=0;
  minimum_scale_factor=FLT_MAX;

  for (x=0;x<deconvolved->GetXSize();x++)
    {
	for (y=0;y<deconvolved->GetYSize();y++)
	{
	  for (z=0;z<deconvolved->GetZSize();z++)
	    {
	      //test response:
	      //if (x==0 && y==0 && z==0)
	      if (mask->GetValue(x,y,z)>0.5)
		{

		  ODFplotter = new ODFPlotter(deconvolved);
		  ODFplotter->SetSmoothOn();
		  ODFplotter->SetVoxel(x,y,z);
		  ODFplotter->SetOpacity(ODF_opacity);
		  ODFplotter->SetColour(colour[0],colour[1],colour[2]);
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
  



  for (i=0;i<counter;i++)
    {
      ODFplotter_array[i]->SetScaleFactor(minimum_scale_factor);
      display_window->AddObject(ODFplotter_array[i]->CreateActor());
      
    }
  
  //DO: delete all ODFPlotters in array, then:  
  free(ODFplotter_array);

//plot:
  display_window->Display();


  

  
  delete(display_window);

    }


  //debug:no range: should add clamp if leave this way
  deconvolved->WriteMincFileNoRange(output_name,clobber,history);


  image_data->Release();
  deconvolved->Release();
  mask->Release();

  //if (remove_zero_b && algorithm!=FORECAST)
  if (algorithm!=FORECAST)
      {
      response_nonzero_b->Release();
      nonzero_b_image->Release();
    }
  if (use_e1)
    {
      e1x->Release();
      e1y->Release();
      e1z->Release();
      response->Release();
    }
  if (algorithm==FORECAST && !def)
    {
      lambda1->Release();
      lambda2->Release();
    }
  if (!def) /*if not using default reponse*/
    {
     response_voxel->Release();
    }
}



void Usage()
{  
   Error("\nUsage: mincDeconvolve -DWIs <DWIs.mnc> -mask <calc_mask.mnc> [options] -o <out.mnc> \n\n-DWIs:  series of DWIs that must have a single nonzero bvalue.  If response_voxel option is used, there must be at last one b=0 frame as well, and it must be at the *beginning* of the series.\n-mask: voxels in which deconvolution will be performed\n\nAll inputs must be sampled the same way.\n\nUse these options if you would like to specify a subject-specific response [by default, a population-based average is provided]:\n\n-response_voxel: a file in which the voxel or voxels that are considered the single fibre response function are nonzero. \n -lambda1 <lambda1.mnc> -lambda2 <lambda2.mnc>: enter the eigenvalues from a previous tensor fit of this data, to be used in the FORECAST deconvolution algorithm (The FORECAST algorithm is the only algorithm currently available - default).\n\nOther options:\n\n-clobber:overwrite existing <out.mnc>\n");

		  /*-noremove_zero_b: input DWIs have no b=0 images: speeds things up to tell the program so; otherwise (default) the program removes them \n");*/
  //this is only needed for the Tournier approach: which isn't finished:
  //-e1 <e1x.mnc> <e1y.mnc> <e1z.mnc> : use this as the single fibre direction (default finds the minimum of the DWI ODF, which might be more error prone than using the fitted tensor).\n


}


int main(int argc, char *argv[])
{


  int i;
  bool normalize=true;
  bool clobber=false;
  char *output_name;
  char *input_name;
  char *mask_name;
  char *response_name=NULL;
  char *response_mask_name=NULL;
  char *lambda1_name=NULL;
  char *lambda2_name=NULL;
  //bool remove_zero_b=true;
  bool use_e1=false;
  char *e1x_name;
  char *e1y_name;
  char *e1z_name;
  bool def=true;


  int order=8;

  float xsize, ysize, zsize;
  char tempdirname[200];
  ImageData *image_data=NULL;
  ImageData *mask_data=NULL;
  bool templatefilename_set=true;
  char command[1000];
  
  char *tempfilename_array[3];
  int j;

  
  
  
  // Check number of input file names 
  if (argc < 3)
    {
      Usage();
    }
  
   for(i = 1; i < argc; i++) 
    {
      if(!strcmp(argv[i], "-help")) 
	{
	  Usage();
	}
      else if (!strcmp(argv[i], "-DWIs")) 
	{
	  if(++i >= argc) Usage();
	  input_name = argv[i];
	}
      else if (!strcmp(argv[i], "-mask")) 
	{
	  if(++i >= argc) Usage();
	  mask_name = argv[i];
	}
      /*
      else if (!strcmp(argv[i], "-response")) 
	{
	  if(++i >= argc) Usage();
	  response_name = argv[i];
	  def = false; //IRL
	}
      */
      else if (!strcmp(argv[i], "-response_voxel")) 
	{
	  if(++i >= argc) Usage();
	  response_mask_name = argv[i];
	  def = false; //IRL
	}
      else if (!strcmp(argv[i], "-e1")) 
	{
	  if(++i >= argc) Usage();
	  e1x_name = argv[i++];
	  e1y_name = argv[i++];
	  e1z_name = argv[i];
	  use_e1=true;
	  def=false;
	}
      else if (!strcmp(argv[i], "-lambda1")) 
	{
	  if(++i >= argc) Usage();
	  lambda1_name = argv[i];
	  def = false; //IRL
	}
      else if (!strcmp(argv[i], "-lambda2")) 
	{
	  if(++i >= argc) Usage();
	  lambda2_name = argv[i];
	  def = false; //IRL
	}
     /* else if(!strcmp(argv[i], "-noremove_zero_b")) 
	{
	  remove_zero_b = false;
	}*/ /*IRL remove b=0 images every time*/
      else if (!strcmp(argv[i], "-o")) 
	{
	  if(++i >= argc) Usage();
	  output_name = argv[i];
	}
      else if(!strcmp(argv[i], "-clobber")) 
	{
	  clobber = true;
	}
      //an unknown option:
      else if(argv[i][0] == '-') 
	{
	  Usage();
	}
      
    }

   
   if (order>8)
     {
       Error("implemented only up to order 8 so far...");
     }

   char *history=CreateHistoryString(argc, argv);    
   

    //for file resampling:
   sprintf(tempdirname,"/tmp/tmp%i",(int) time(NULL));
   sprintf(command,"mkdir %s\n",tempdirname);
   system(command);
  
   image_data=new ImageData(input_name); /*probably not the most efficient way but easier to get dim sizes*/
   mask_data=new ImageData(mask_name);
   
   /*check if data needs to be resampled*/
   for (j=0;j<3;j++)
   {
     tempfilename_array[j]=(char *) malloc(200*sizeof(char));
   }
   
      
  //Check if mask has to be resampled
  if(mask_data->GetXSize()!=image_data->GetXSize() || mask_data->GetYSize()!=image_data->GetYSize() || mask_data->GetZSize()!=image_data->GetZSize())
  {
     Resample(mask_name,tempfilename_array[0],tempdirname,0,false,input_name,templatefilename_set,true);
  }else
  {
     tempfilename_array[0]=mask_name;
  }
  if(!def)
  {
/*
     ImageData *response_data = new ImageData(input_name);
     
     
     //assume the e1 and lambda files are resampled correctly since they were probably created from the diffusion set, should probably check this
     
     if(response_data->GetXSize()!=image_data->GetXSize() || response_data->GetYSize()!=image_data->GetYSize() || response_data->GetZSize()!=image_data->GetZSize())
      {
       Resample(response_name,tempfilename_array[1],tempdirname,0,false,input_name,templatefilename_set,true);
      }else
      {
	 tempfilename_array[1]=response_name;
      }
*/
     ImageData *response_mask_data = new ImageData(response_mask_name);
     if(response_mask_data->GetXSize()!=image_data->GetXSize() || response_mask_data->GetYSize()!=image_data->GetYSize() || response_mask_data->GetZSize()!=image_data->GetZSize())
      {
       Resample(response_mask_name,tempfilename_array[2],tempdirname,0,false,input_name,templatefilename_set,true);

      }else
      {
      tempfilename_array[2]=response_mask_name;
      }
  }
  
/*   Deconvolve(input_name, mask_name, output_name, response_name, response_mask_name, history, clobber, use_e1, e1x_name, e1y_name, e1z_name, order, lambda1_name, lambda2_name, def);*/

  
  
  Deconvolve(input_name,  tempfilename_array[0], output_name,  tempfilename_array[1],  tempfilename_array[2], history, clobber, use_e1, e1x_name, e1y_name, e1z_name, order, lambda1_name, lambda2_name, def);

      for (j=0;j<3;j++)
      {
	 //free(tempfilename_array[j]); //IRL this is causing glibc errors?
      } 


   return(0);
}


