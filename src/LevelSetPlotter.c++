/************************************************************************

   File: LevelSetPlotter.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __LevelSetPlotter_h
#include "LevelSetPlotter.h" 
#endif

#ifndef __otherlibs_include_h
#include "otherlibs_include.h"
#endif

#ifndef __misc_h
#include "misc.h"
#endif

#ifndef __ImageData_h
#include "ImageData.h"
#endif


LevelSetPlotter::LevelSetPlotter(ImageData *image_data, float low_thresh, float high_thresh)
{
  _image_data=image_data;
  switch (_image_data->GetImageType())
    {
    case Image3D_vanilla:
      break;
    default:
      Error("Image type unsupported."); //IRL Error for gcc4
    }

  _opacity=0.5;
  _red=0.8;
  _green=0.8;
  _blue=0.8;
  _low_thresh=low_thresh;
  _high_thresh=high_thresh;
  _smoothing_iterations=100;
  _actor_init=false;
  

}

LevelSetPlotter::LevelSetPlotter(char *image_filename, float low_thresh, float high_thresh)
{


  _image_data=new ImageData(image_filename);
  switch (_image_data->GetImageType())
    {
    case Image3D_vanilla:
      break;
    default:
      Error("Image type unsupported.");//IRL Error for gcc4
    }

  _opacity=0.5;
  _red=0.8;
  _green=0.8;
  _blue=0.8;
  _low_thresh=low_thresh;
  _high_thresh=high_thresh;
  _smoothing_iterations=100;
  _actor_init=false;
  

}

/* this isn't done and isn't necessary
LevelSetPlotter::LevelSetPlotter(char *image_filename, float low_thresh, float high_thresh, char *tempdirname)
{

  //resample file so that there are guaranteed to be zeros around the level set:
  //this assumes temp dir cleaned up later by calling program.
  //this ring around edge *could* be taken care of in the assignment to the vtk image data below.

  int seed;
  char *tempfilename=(char *) malloc(200*sizeof(char));
  ImageData *image_data;
  int xnelements,ynelements,znelements;
  float xstart, ystart, zstart;

  seed = (int) time (NULL);
  srand48(seed); 

  sprintf(tempfilename,"%s/tmp%f.mnc",tempdirname, (double) rand());

  image_data=new ImageData(image_filename);

  //DO: I think this should be changed to just check that it is 3D:
  switch (image_data->GetImageType())
    {
    case Image3D_vanilla:
      break;
    default:
      Error("Image type unsupported.");
    }

  xnelements=(int) image_data->GetXSize()+4;
  ynelements=(int) image_data->GetYSize()+4;
  znelements=(int) image_data->GetZSize()+4;


  if (image_data->GetXStep()


  _image_data=new ImageData(tempfilename);



  _opacity=0.5;
  _red=0.8;
  _green=0.8;
  _blue=0.8;
  _low_thresh=low_thresh;
  _high_thresh=high_thresh;
  _smoothing_iterations=100;
  _actor_init=false;
  

}
*/


LevelSetPlotter::LevelSetPlotter(ImageData *image_data)
{
  _image_data=image_data;
  switch (_image_data->GetImageType())
    {
    case Image3D_vanilla:
      break;
    default:
      Error("Image type unsupported.");//IRL Error for gcc4
    }

  _opacity=0.5;
  _red=0.8;
  _green=0.8;
  _blue=0.8;
  _low_thresh=0.5;
  _high_thresh=0.5;
  _smoothing_iterations=100;
  _actor_init=false;

}



LevelSetPlotter::~LevelSetPlotter()
{
  _actor->Delete();
  
}


void LevelSetPlotter::SetOpacity(double opacity) //ilana vtk5 SetOpacity takes double
{
  _opacity=opacity;
  
}


void LevelSetPlotter::SetColour(double red, double green, double blue)
{
  _red=red;
  _green=green;
  _blue=blue;
  
}

vtkActor *LevelSetPlotter::CreateActor()
{

  int i,j,k;
  float value;
  int offset;
  int x,y,z;
  

  float xstep,ystep,zstep;
  float xstart, ystart, zstart;
  bool xneg,yneg,zneg;
  

  xstart=_image_data->GetXStart();
  xstep=_image_data->GetXStep();
  ystart=_image_data->GetYStart();
  ystep=_image_data->GetYStep();
  zstart=_image_data->GetZStart();
  zstep=_image_data->GetZStep();

  if (xstep<0) {
    xstart=xstart+xstep*(_image_data->GetXSize()-1);
    xstep=-xstep;
    xneg=true;
    
  }
 
  if (ystep<0) {
    ystart=ystart+ystep*(_image_data->GetYSize()-1);
    ystep=-ystep;
    yneg=true;
  }
 
  if (zstep<0) {
    zstart=zstart+zstep*(_image_data->GetZSize()-1);
    zstep=-zstep;
    zneg=true;
  }
    
    

  vtkStructuredPoints *volume = vtkStructuredPoints::New();

  //a buffer to ensure that the level set is surrounded by zeros is included here

  volume->SetDimensions(_image_data->GetXSize()+4, _image_data->GetYSize()+4, _image_data->GetZSize()+4);

  volume->SetSpacing(xstep,ystep,zstep);
    
  //volume->SetOrigin( xstart, ystart, zstart);

  //offset the origin:

  volume->SetOrigin( xstart-2*zstep, ystart-2*ystep, zstart-2*zstep);
  
  volume->SetScalarTypeToFloat(); //Ilana add this here too, vtk5 seems to need this

  vtkFloatArray *scalars=vtkFloatArray::New();


  scalars->Allocate((_image_data->GetXSize()+4)*(_image_data->GetYSize()+4)*(_image_data->GetZSize()+4));  

  for (k=0;k<_image_data->GetZSize()+4;k++)
    {
      for (j=0;j<_image_data->GetYSize()+4;j++)
	{
	  for (i=0;i<_image_data->GetXSize()+4;i++)
	    {
	      if (i<2 || i>_image_data->GetXSize()+1 || j<2 || j>_image_data->GetYSize()+1 || k<2 || k>_image_data->GetZSize()+1)
		{
		  value=0;
		}
	      else 
		{
		  if (xneg)
		    {
		      x=_image_data->GetXSize()-1-(i-2);
		    }
		  else
		    {
		      x=i-2;		  
		    }
		  if (yneg)
		    {
		      y=_image_data->GetYSize()-1-(j-2);
		    }
		  else
		    {
		      y=j-2;		  
		    }
		  if (zneg)
		    {
		      z=_image_data->GetZSize()-1-(k-2);
		    }
		  else
		    {
		      z=k-2;		  
		    }
		 

		  value=_image_data->GetValue(x,y,z);
		}
	      
	      offset=i+j*(_image_data->GetXSize()+4)+k*(_image_data->GetXSize()+4)*(_image_data->GetYSize()+4);	      
	      if (value>=_low_thresh && value<=_high_thresh)
		{
		  scalars->InsertValue(offset,-1.0);
		}
	      else
		{
		  scalars->InsertValue(offset,0.0);
		}
	      
	      
	    }
	}
    }

  volume->GetPointData()->SetScalars(scalars);
  scalars->Delete(); //ilana do some clean up
  
  // Run ContourFilter
  vtkContourFilter *aContour= vtkContourFilter::New();
  aContour->SetInput(volume);
  aContour->ComputeScalarsOff();
  aContour->ComputeGradientsOff();
  aContour->ComputeNormalsOff();
  aContour->UseScalarTreeOn();
  aContour->SetValue(1, -1.0);

  //Smooth the data
  vtkSmoothPolyDataFilter *smoother;
  smoother = vtkSmoothPolyDataFilter::New();
  smoother->SetInput(aContour->GetOutput());
  smoother->SetNumberOfIterations(_smoothing_iterations);

 
  // Compute Normals
  vtkPolyDataNormals *normals;
  normals = vtkPolyDataNormals::New();
  normals->SetInput(smoother->GetOutput());

   
  // Create triangle strips
  vtkStripper *stripper;
  stripper = vtkStripper::New();
  stripper->SetInput(smoother->GetOutput());



  //  map to graphics library
  vtkPolyDataMapper *level_set_mapper=vtkPolyDataMapper::New();
  level_set_mapper->SetInput(stripper->GetOutput());


  // actor coordinates geometry, properties, transformation
  vtkActor *level_set_actor=vtkActor::New();
  level_set_actor->SetMapper(level_set_mapper);
  level_set_actor->GetProperty()->SetColor(_red, _green, _blue);
  level_set_actor->GetProperty()->SetOpacity(_opacity);
  

  if (_actor_init)
    {
      _actor->Delete();
    }
  
  _actor=level_set_actor;
  _actor_init=true;
  
  return level_set_actor;
  


}

vtkActor *LevelSetPlotter::GetActor()
{
  if (!_actor_init)
    {
      return this->CreateActor();
      
    }
  else
    {
      return _actor;
      
    }
  
}

void LevelSetPlotter::SetHighThreshold(float high_thresh)
{
  _high_thresh=high_thresh;
  
}

float LevelSetPlotter::GetHighThreshold()
{
  return _high_thresh;
  
}



void LevelSetPlotter::SetLowThreshold(float low_thresh)
{
  _low_thresh=low_thresh;
  
}

float LevelSetPlotter::GetLowThreshold()
{
  return _low_thresh;
  
}

