/************************************************************************

   File: ODFPlotter.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __ODFPlotter_h
#define __ODFPlotter_h

#ifndef __otherlibs_include_h
#include "otherlibs_include.h"
#endif

#ifndef __ImageData_h
#include "ImageData.h"
#endif

#ifndef __Directions_h
#include "Directions.h"
#endif


#ifndef __misc_h
#include "misc.h"
#endif

#ifndef __S2Interpolator_h
#include "S2Interpolator.h"
#endif

/*

this nly works for symmetric ODFs right now.


display arbitrary function on S2 given ODFdata ImageData object. plots
in minc space


an object to plot a single ODF: invoke
multiple times to do many.  return a pointer to an actor.  add that actor 
to the DisplayWindow

a colour option could be to colour code the surface according to direction: RGB proportional to xyz.  meaning what exactly?  by mag or does vector have to have length 1?  could this work instead of stretching?

currently, _interpolation_kernel_window_radius=_interpolation_kernel_FWHM=X*max expected separation (empirical)
Sugg X=1.5, but play with this
DO: add methods to control smoothing via interpolation kernel, or other means

user gives degree inputs but internally, use radians


this assumes the directions in the _ODF_image_data are isotropically spaced!

stretch factor is hardcoded as 4: could be input


DO: make CalculateInterpolationLookupTable be a separate method. maybe even make  InterpolationLookupTable be an object with ->Calculate().  CalculatePlottingPoints() could call it if not set already, but have optin to change it without calling CalculatePlottingPoints()

RGB is not implemented yet.

DO: figure out how to smooth the surface after plotting it (Gourand shading??)

confusing: will be setting the same interpolator for many ODFPlotters.  when do I destroy it?  when destroy the interpolator.  but the interpolator is not known to the end user?  perhaps it is.

S2Interpolator should probably be deep copied and destroyed with the ODFPlotter: included in the copy constructor would be:

gsl_matrix_float_memcpy(new_interpolation_lookup_table,_interpolation_lookup_table);



currently, must get the interpolator from one of them and destroy it in calling program

angle_step is now fixed and set when initialized.

DO:SetInterpolator should check that the directions match, or at least that the number of directions match

Gourand shading doesn't appear to make a difference

ODF opacity?

 */


typedef enum colour_scheme {RGB, SingleColour, Radial} COLOUR_SCHEME;
typedef enum display_mode {Stretch,SubtractMinimum,Physical} DISPLAY_MODE;



class ODFPlotter
{
 public:
  ODFPlotter(ImageData *ODF_image_data);
  ODFPlotter(ImageData *ODF_image_data, ImageData *polarity, ImageData *full_pdf);
  ODFPlotter(ImageData *ODF_image_data, float angle_step);
  ~ODFPlotter();
  void SetVoxel(int voxel_x, int voxel_y, int voxel_z);
  void CalculatePlottingPoints();
  float GetBestScaleFactorForVoxelFit();
  void SetScaleFactor(float scale_factor);
  vtkActor *CreateActor();
  void SetColour(float red, float green, float blue);
  void SetColourScheme(COLOUR_SCHEME colour_scheme);
  void SetDisplayMode(DISPLAY_MODE display_mode);
  S2Interpolator *GetInterpolator();
  S2Interpolator *GetInterpolator2(); //for two input files
  void SetInterpolator(S2Interpolator *interpolator);
  void SetInterpolator2(S2Interpolator *interpolator_2); //for two input files
  void SetSmoothOn();
  void SetOpacity(double opacity); //ilana vtk5 SetOpacity requires double
  //void SetPlotFan(int plot_fan); //should input to constructor along with the filenames

 private:
  ImageData *m_ODF_image_data;
  float m_angle_step;
  COLOUR_SCHEME m_colour_scheme;
  float m_single_colour[3]; 
  DISPLAY_MODE m_display_mode;
  int m_voxel[3];
  gsl_vector_float *m_rvalues;
  bool m_rvalues_are_allocated;
  float m_scale_factor;
  S2Interpolator *m_interpolator;
  S2Interpolator *m_interpolator_2;
  int m_smoother_iterations;
  float m_specular_fraction;
  float m_diffuse_fraction;
  float m_ambient_fraction;
  float m_specular_power;
  float m_specular_color_red;
  float m_specular_color_green;
  float m_specular_color_blue;
  bool m_smooth_shapes;
  double m_opacity; //ilana in VTK5 SetOpacity requires double
  ImageData *m_full_pdf;
  ImageData *m_polarity;
  int m_plot_fan;
  
};



#endif
