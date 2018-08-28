/************************************************************************

   File: SphericalHarmonics.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2008 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/

#ifndef __SphericalHarmonics_h
#define __SphericalHarmonics_h

#include <gsl/gsl_matrix.h>
#include <complex>

#include "otherlibs_include.h"
#include "ImageData.h"
#include "misc.h"

#ifndef __misc_h
#include "misc.h"
#endif


/*
there is a lot of common stuff in CalculateSHCoeffs CalculateSHFit CoeffsToImage.  this should be put in private common fncs

 */
//need full def. of what's needed here...will need to be able to output an array (not an ImageData) as well probably...           

class SphericalHarmonics
{
 public:
  SphericalHarmonics();
  ~SphericalHarmonics();
  void CalculateSHFit(ImageData *input_data_sphere,ImageData *output_data_sphere, int x,int y, int z);
  //coeffs as gsl matrix for one voxel:
  gsl_matrix *CalculateSHCoeffs(ImageData *input_data_sphere, int x,int y, int z);

  //coeffs as ImageData: think about format:
  void CalculateSHCoeffs(ImageData *input_data_sphere,ImageData *output_coeffs, int x,int y, int z);
  //for doing the whole ImageData: output coeffs:
  ImageData *CalculateSHCoeffs(ImageData *input_data_sphere);
  //this assumes the ImageDatas involved are "correct" for *this* SphericalHarmonics
  void CoeffsToImage(ImageData *output_data_sphere, ImageData *coeffs, int x,int y, int z);
  void SetOrder(int order);
  int GetRank();
  void GetLeverageTerms(float *lev, int x, int y, int z, int N);

private:
  std::complex<float> GetSH(int _l, int _m, float theta, float phi);
  bool m_sh_basis_is_computed;
  gsl_matrix *m_BasisFunction;
  gsl_matrix *m_ls;
  int m_order;
  int m_rank;
};

#endif

