/************************************************************************

   File: S2Interpolator.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __S2Interpolator_h
#define __S2Interpolator_h

#ifndef __otherlibs_include_h
#include "otherlibs_include.h"
#endif

//#ifndef __ImageData_h
//#include "ImageData.h"
//#endif

#ifndef __Directions_h
#include "Directions.h"
#endif


#ifndef __misc_h
#include "misc.h"
#endif



/*

if LUT is not set, will calculate the weight.  if it is set, will use it.  
DO: request for an interpolated value without using the LUT even whenthe LUT is set

user must call CalculateLUT if wants to benefit from speed thereof

FWHM is not used for linear: only radius

the directions (at which one might want to calculate interpolated values) are fixed for this interpolator.


this is assumes the directions of the input data are only in one hemisphere (and function symmetric)

 */

class ImageData;

typedef enum interpolation_type {Gaussian, Linear} INTERPOLATION_TYPE;

class S2Interpolator
{

 public:
  //S2Interpolator(ImageData *image_data); //this gives it the directions and the data
 
  S2Interpolator(ImageData *image_data, Directions *directions, float kernel_window_FWHM, float kernel_window_radius, INTERPOLATION_TYPE interpolation_type); //this gives it the directions and the data and the places at which we will want interpolated values and the type of interpolation, etc
  

  ~S2Interpolator();
  //gsl_matrix_float *GetLookupTable();
  //void SetLookupTable(gsl_matrix_float *interpolation_lookup_table);
  float GetInterpolatedValue(int counter, int voxel_x, int voxel_y, int voxel_z);
  void CalculateLookupTable();
  void SetInterpolationType(INTERPOLATION_TYPE interpolation_type);
  void SetFWHM(float kernel_FWHM);
  void SetWindowRadius(float kernel_window_radius);
  Directions *GetDirections();
  void SetAssymetric();

 private:

  ImageData *m_image_data;
  gsl_matrix_float *m_lookup_table; 
  bool m_lookup_table_is_current;
  bool m_lookup_table_is_allocated;
  INTERPOLATION_TYPE m_interpolation_type;
  float m_kernel_window_radius;
  float m_kernel_FWHM;
  void CalculateWeights(int counter);
  Directions *m_directions;
  bool m_assymetric;


};



#endif
