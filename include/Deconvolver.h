/************************************************************************

   File: Deconvolver.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2007 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __Deconvolver_h
#define __Deconvolver_h


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
could use the same SphericalHarmonics object for in and response.  Haven't here. 

order should vary with N...

currently does only spherical deconvolution

SetInput must be called before SetResponseFunction and both must be called just once, as is.

output directions are input; could change to be something else

assumes input has no b=0: single nonzero b value ODF


*/

typedef enum eDeconvolutionMode {Spherical} DeconvolutionMode;
typedef enum eDeconvolutionAlgorithm {Tournier, FORECAST} DeconvolutionAlgorithm;//note: Tournier doesn't work right now

class Deconvolver
{
 public:
  Deconvolver();
  ~Deconvolver();
  void SetDeconvolutionMode(DeconvolutionMode deconvolution_mode);
  void SetInput(ImageData *ODF_input);
  void SetResponseFunctionSingleVoxel(ImageData *response_function, ImageData *mask_response_voxel);
  void SetResponseFunctionSingleVoxel(ImageData *response_function, int x, int y, int z);
  void SetResponseFunctionSingleVoxel(ImageData *response_function, ImageData *mask_response_voxel, ImageData *e1x, ImageData *e1y, ImageData *e1z);
  void SetResponseFunction(ImageData *response_function);
  void Deconvolve(int x,int y,int z);
  ImageData *GetOutput();
  void SetOrder(int order);
  void SetDeconvolutionAlgorithm(DeconvolutionAlgorithm algorithm);
  void SetLambdas(ImageData *lambda1,ImageData *lambda2,ImageData *mask_response_voxel);
  void SetLambdasDefault(); /*IRL*/
  void GetEigenvaluesAndSo(float &lambda1, float &lambda2, float &So);
  void SetLambdas(float lambda1, float lambda2, float So);
  
 private:
  void Deconv2Sphere(int x,int y,int z);
  void RotateResponse(VECTOR3D vector);
  float DeltaCoeff(int l);
  float A(int l);
  int thisl(int n);
  float Scale(int l);
  ImageData *m_ODF_input;
  DeconvolutionMode m_deconvolution_mode;
  ImageData *m_deconvolved;
  ImageData *m_response_function;
  int m_order;
  ImageData *m_in_SH;
  ImageData *m_response_SH;
  ImageData *m_deconvolved_SH;
  SphericalHarmonics *m_in_sphere_harm;
  SphericalHarmonics *m_response_sphere_harm;
  bool m_in_sphere_harm_constructed;
  bool m_response_sphere_harm_constructed;
  bool m_response_SH_calculated;
  float m_rotation[3][3];
  DeconvolutionAlgorithm m_algorithm;
  float m_lambda1;
  float m_lambda2;
  float m_b;
  float m_response_So;
};

#endif
