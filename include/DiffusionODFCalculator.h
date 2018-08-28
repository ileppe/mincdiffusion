/************************************************************************

   File: DiffusionODFCalculator.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __DiffusionODFCalculator_h
#define __DiffusionODFCalculator_h

#ifndef __otherlibs_include_h
#include "otherlibs_include.h"
#endif

#ifndef __ImageData_h
#include "ImageData.h"
#endif

#ifndef __misc_h
#include "misc.h"
#endif

#ifndef __S2Interpolator_h
#include "S2Interpolator.h"
#endif

#ifndef __Directions_h
#include "Directions.h"
#endif

#ifndef __ProcessedDTIData_h
#include "ProcessedDTIData.h"
#endif




/*
currently: give directions at which you will want the ODF calculated to constructor: calculates all at once, for a requested voxel, and writes to _ODF_image_data.  Can Get() this ODF ImageData at any time.  from calling program, can control which voxels with a mask

options: QBI, qspace with integration, DTI: analytic


normalize to different things: trace, mean, mean square: think about and see notes and email to Peter.  

can SetTraceFile, etc, for these things

could also mult by GFA as Tuch does / use either my or Tuch's GFA

problem: the directions may come from someshere else.  freeing the imagedata would free the directions: check on this: do reference counting, or do deep copies

//angle step for great circle will be 0.5*sqrt(4pi/(2*[number of input directions]))

DO: implement min-max normalization and scaling by GFA as in Tuch MRM 2004

DO: add normalization by division by a 3D ImageData (e.g., trace.....)


could add thresholding (currently do afterwards with mincmath)

 */

typedef enum eDiffusionODFCalculationMode {QBI, QBISH, DTI} DiffusionODFCalculationMode;

typedef enum eNormalizationMode {mean,volume} NormalizationMode;
//could add trace: not sure it makes sense: use volume.


class DiffusionODFCalculator
{
 public:
  DiffusionODFCalculator(ImageData *input_DWI_data, Directions *directions_to_calculate);
  DiffusionODFCalculator(ProcessedDTIData *input_DTI_data, Directions *directions_to_calculate);
  ~DiffusionODFCalculator();
  void InitializeToZero();
  ImageData *GetODFImageData();
  void CalculateODF(int voxel_x, int voxel_y, int voxel_z);
  void NormalizeODF(NormalizationMode normalization_mode, int voxel_x, int voxel_y, int voxel_z);
  void SetCalculationMode(DiffusionODFCalculationMode mode);


  
  //this is currently determined by the constructor:
  //void DiffusionODFCalculator::SetCalculationMode(DIFFUSION_ODF_CALCULATION_MODE diffusion_ODF_calculation_mode);
  //this is currently just private:
  //S2Interpolator **DiffusionODFCalculator::GetInterpolatorArray();
  

 private: 
  ImageData *m_ODF_image_data;
  S2Interpolator **m_interpolator_array;
  ImageData *m_input_DWI_data;
  ProcessedDTIData *m_input_DTI_data;
  DiffusionODFCalculationMode m_diffusion_ODF_calculation_mode;
  bool m_interpolator_array_is_set;
  float m_angle_step_for_great_circle;
  ImageData *m_no_zero_b;
  gsl_matrix *m_basis;
  int m_rank;
  SphericalHarmonics *m_sphere_harm;
  bool m_sphere_harm_constructed;
  gsl_matrix *m_P;
  gsl_matrix *m_basis2;

  void CalculateQBIODF(int voxel_x, int voxel_y, int voxel_z);
  void CalculateQBISHODF(int voxel_x, int voxel_y, int voxel_z);
  void CalculateDTIODF(int voxel_x, int voxel_y, int voxel_z);
  

};

#endif
