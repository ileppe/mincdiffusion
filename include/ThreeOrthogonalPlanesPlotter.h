/************************************************************************

   File: ThreeOrthogonalPlanesPlotter.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __ThreeOrthogonalPlanesPlotter_h
#define __ThreeOrthogonalPlanesPlotter_h

#ifndef __otherlibs_include_h
#include "otherlibs_include.h"
#endif

#ifndef __ImageData_h
#include "ImageData.h"
#endif

#ifndef __misc_h
#include "misc.h"
#endif



/*
note: you can do a lot more with this: see vtk example ImagePlaneWidget.py: make GUi buttons in renderer (would do in display window) that change the lookup table, etc for the widgets.  could also add to this class methods for changing the LUT, etc, before give it to the DisplayWindow renderer

may want a method to get the picker so that can have the same picker for any other datasets: not sure

DO: look into problem where if zoom leaves no outside the image space, can't zoom, rotate, etc anymore

DO: see what happens if plot two datasets: if aligned originally, should be able to move together (?).  also have a LUT specification here
 */


class ThreeOrthogonalPlanesPlotter
{
 public:
  ThreeOrthogonalPlanesPlotter(ImageData *image_data);
  ~ThreeOrthogonalPlanesPlotter();
  vtkImagePlaneWidget **GetWidgets();
  

 private:
  vtkImageData *_image_data_vtk;
  //  vtkImagePlaneWidget **_image_plane_widgets;   // old
  vtkImagePlaneWidget * _image_plane_widgets[3];    // CLAUDE
  
};



#endif
