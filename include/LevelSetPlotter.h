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
#define __LevelSetPlotter_h 

#ifndef __otherlibs_include_h
#include "otherlibs_include.h"
#endif

#ifndef __misc_h
#include "misc.h"
#endif

#ifndef __ImageData_h
#include "ImageData.h"
#endif

/*
only for Image3D_vanilla type of ImageData

may require positive steps!?! in vtkStructuredPoints

could add setting threshold after initialization...

 */

class LevelSetPlotter
{
 public:
  LevelSetPlotter(ImageData *image_data, float low_thresh, float high_thresh);
  LevelSetPlotter(char *image_filename, float low_thresh, float high_thresh);
  LevelSetPlotter(ImageData *image_data);
  ~LevelSetPlotter();
  void SetOpacity(double opacity); //ilana VTK5 SetOpacity takes double
  void SetColour(double red, double green, double blue); /*set this to double for vtk5 ilana*/
  vtkActor *CreateActor();
  vtkActor *GetActor();
  void SetHighThreshold(float high_thresh);
  float GetHighThreshold();
  void SetLowThreshold(float low_thresh);
  float GetLowThreshold();
  

 private:
  double _opacity; //ilana VTK5 SetOpacity takes double
  double _red; /*set colors to double for vtk5 ilana*/
  double _green;
  double _blue;
  float _low_thresh; 
  float _high_thresh;
  ImageData *_image_data;
  int _smoothing_iterations;
  vtkActor *_actor;
  bool _actor_init;
  

}
;

#endif
