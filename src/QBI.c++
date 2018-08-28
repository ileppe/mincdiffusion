/************************************************************************

   File: QBI.c++

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

/*
directions text file must be of format dirx diry dirz\n..., and must be unit length


could add to this: computation and output of mean diffusivity, GFA, maxima (these may be computed in DiffusionODFCalculator, but probably elsewhere would be better.


fixing seg fault when alloc ImageData

 */

void CalculateQBIODFs(char *input_filename, char *mask_filename, char *output_filename, char *directions_filename, char *history, bool clobber, bool normalize, bool SH)
{

  int x,y,z;
  int n;
  Directions *directions=new Directions(directions_filename);



  ImageData *DWI_image_data = new ImageData(input_filename);




  ImageData *mask_image_data = new ImageData(mask_filename);  
 
  //create no zero b volume here: DiffusionODFCalculator no longer does (in the SH case)
  ImageData *no_zero_b=DWI_image_data->CreateDWIODF();


  //DiffusionODFCalculator *diffusion_ODF_calculator=new DiffusionODFCalculator(DWI_image_data, directions);
  DiffusionODFCalculator *diffusion_ODF_calculator=new DiffusionODFCalculator(no_zero_b, directions);

  if (SH)
    {
      diffusion_ODF_calculator->SetCalculationMode(QBISH);
    }
  else
    {
      diffusion_ODF_calculator->SetCalculationMode(QBI);
    }
  //cout << "initialize QBI volume to zero\n";
  //this takes too long: if you don't care, comment it out

  //diffusion_ODF_calculator->InitializeToZero();



  for (x=0;x<DWI_image_data->GetXSize();x++)
    {
	for (y=0;y<DWI_image_data->GetYSize();y++)
	{
	  for (z=0;z<DWI_image_data->GetZSize();z++)
	  
	    
	    {
	      
	      if (mask_image_data->GetValue(x,y,z)>0.5)
		{
		  
		  diffusion_ODF_calculator->CalculateODF(x,y,z);
		  //not normalizing because it is failing (although NormalizeODF works) 
		  //cout << diffusion_ODF_calculator->GetODFImageData()->GetValue(x,y,z,25) << endl;
		  if (normalize)//again having a problem here: hangs.  
		    {
		      
		      //diffusion_ODF_calculator->NormalizeODF(volume,x,y,z);
		    }


		  //test
		  //zero here:
		  //cout << diffusion_ODF_calculator->GetODFImageData()->GetValue(x,y,z,25) << endl;

		}
	      
	      else
		{
		  //set background to 0 here:
		  for (n=0;n<directions->GetNumberOfDirections();n++)
		    {

		      diffusion_ODF_calculator->GetODFImageData()->SetValue(x,y,z,n,0.0);

		    }
		  
		}
	      
	     
	      
	      
	      
	    }
  
	}
    }
  


  //ImageData *QBI_image_data= diffusion_ODF_calculator->GetODFImageData();
  //QBI_image_data->WriteMincFile(output_filename, clobber, history, false); // IL changed this to inculde bbox option
  //diffusion_ODF_calculator->GetODFImageData()->WriteMincFile(output_filename, clobber, history,false);
  //JC the calculation of min and max and input of these corrupts ODF files(??)



  diffusion_ODF_calculator->GetODFImageData()->WriteMincFileNoRange(output_filename, clobber, history);

  //free(history);
  //directions->Release();
  //DWI_image_data->Release();
  //mask_image_data->Release();
  //delete(diffusion_ODF_calculator);
  //QBI_image_data->Release();
  no_zero_b->Release();
  
  
}

void Usage()
{
  Error("Usage: QBI [options] -mask <mask.mnc> -directions <directions.txt> <input.mnc> <output.mnc>\n\noptions:\n\n-clobber: overwrite output file name\n-nonormalize don't normalize ODFs to unit volume\n-SH use spherical harmonic analytic formula (default: numeric QBI reconstruction)\n\n\n inputs must be resampled to same space");
  
}



int main(int argc, char *argv[])
{
  bool clobber=false;
  char *output_name;
  char *input_name;
  char *mask_name;
  char *directions_name;
  int i;
  bool got_input_name=false;
  bool normalize=true;
  bool SH=false;

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
      else if(!strcmp(argv[i], "-clobber")) 
	{
	  clobber = true;
	}
      else if(!strcmp(argv[i], "-SH")) 
	{
	  SH = true;
	}
      else if(!strcmp(argv[i], "-nonormalize")) 
	{
	  normalize = false;
	}
      else if (!strcmp(argv[i], "-mask")) 
	{
	  if(++i >= argc) Usage();
	  mask_name = argv[i];
	}
      else if (!strcmp(argv[i], "-directions")) 
	{
	  if(++i >= argc) Usage();
	  directions_name = argv[i];
	}
      //an unknown option:
      else if(argv[i][0] == '-') 
	{
	  Usage();
	}
      
      //an input and output file names
      else if (!got_input_name)
	{
	  input_name=argv[i];
	  got_input_name=true;
	}
      else if (got_input_name)
	{
	  output_name=argv[i];
	}
      
      
      
    }
  
  
  char *history=CreateHistoryString(argc, argv);    

  //assumes mask file is resampled to match. but could check and resample it here if it doesn't

  if (normalize)
    {
      cout << "normalize feature is currently broken: output will not be normalized\n";
      
    }

  CalculateQBIODFs(input_name, mask_name, output_name, directions_name, history, clobber, normalize, SH);

  return(0);
}




