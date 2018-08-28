/************************************************************************

   File: transform_diff_dirs.c++

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

#include <minc_1_rw.h>
/*
add history to the minc file (?)

note: this overwrites the input file's directions: currently no option to create a new file (could add this option by copying the input file before write to header and having an output filename.  

this probbaly isn't handling stretching/scrunching properly: check it


 */

void Usage()
{
  Error("Usage: transform_diff_dirs -xfm <transform.xfm> <input.mnc>");

}


void transform_diffusion_directions(char *input_name, char *transform_name, char *history)
{
  ImageData *image_data=new ImageData(input_name);
  
  Directions *new_directions;
  
  switch (image_data->GetImageType())
    {
    case DWIdata:
      new_directions=image_data->GetDiffusionParameters()->GetDirections()->Transform(transform_name);
      break;
    case ODFdata:
      new_directions=image_data->GetDirections()->Transform(transform_name);
      break;
      
    }
  
  minc::minc_1_reader rdr;
  rdr.open(input_name,false,false,true);
  new_directions->WriteToMincFile(rdr);//true

  //DO: write history string to minc file (with minc_modify_header)

}


int main(int argc, char *argv[])
{
  char *input_name;
  bool got_input_name=false;
  char *transform_name;
  int i;
  

  if (argc < 4)
    {
      Usage();
    }


  for(i = 1; i < argc; i++) 
    {
      if(!strcmp(argv[i], "-help")) 
	{
	  Usage();
	}
      else if(!strcmp(argv[i], "-xfm")) 
	{
	  if(++i >= argc) Usage();
	  transform_name = argv[i];
	}
      else if(argv[i][0] == '-') 
	{
	  Usage();
	}
      else if (!got_input_name)//if(argv[i][0] != '-')
	{
	  input_name=argv[i];
	  got_input_name=true;
	}
     
      
    }
  
  char *history=CreateHistoryString(argc, argv);

  transform_diffusion_directions(input_name,transform_name,history);
  
  free(history);
  
  return(0);

}
