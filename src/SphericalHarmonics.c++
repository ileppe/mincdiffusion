/************************************************************************

   File: SphericalHarmonics.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2008 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#include <float.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

#include "SphericalHarmonics.h"
#include "otherlibs_include.h"



SphericalHarmonics::SphericalHarmonics()
{
  m_sh_basis_is_computed=false;
  //default:
  m_order=8;
  m_rank=(m_order+1)*(m_order+2)/2;
}

SphericalHarmonics::~SphericalHarmonics()
{
  if (m_sh_basis_is_computed)
    {
      gsl_matrix_free(m_BasisFunction);
      gsl_matrix_free(m_ls);
    }
}



//this will output fitted values sampled at same as input sampling points:output_data_sphere must be allocated already
//input must have one value only if a DWI file

void SphericalHarmonics::CalculateSHFit(ImageData *input_data_sphere,ImageData *output_data_sphere, int x,int y, int z)
{
  int n;

  int n_s = input_data_sphere->GetDirections()->GetNumberOfDirections();

  gsl_matrix *fitted=gsl_matrix_alloc(n_s,1);

  gsl_matrix *C=CalculateSHCoeffs(input_data_sphere,x,y,z);

  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m_BasisFunction,C,0.0,fitted);

  //cout << "out:\n";
  for(n=0;n<input_data_sphere->GetDirections()->GetNumberOfDirections();n++)
    {
      //output_data_sphere->SetValue(x,y,z,n,(float) gsl_matrix_get(fitted,0,n));        
      output_data_sphere->SetValue(x,y,z,n,(float) gsl_matrix_get(fitted,n,0));        
      //cout << (float) gsl_matrix_get(fitted,n,0) << endl;
    }  

  
  output_data_sphere->SetDirections(input_data_sphere->GetDirections());//Retain()s them

  gsl_matrix_free(fitted);


  //print coeffs:
  /*
  cout << "coeffs ";
  for (n=0;n<m_rank;n++)
    {
      cout << gsl_matrix_get(C, n, 0) << " ";
    }
  cout << endl << endl;
  */
  gsl_matrix_free(C);
  
}


gsl_matrix *SphericalHarmonics::CalculateSHCoeffs(ImageData *input_data_sphere, int x,int y, int z)
{

  //NOTE: regularization is commented out here; have done it with regularization

  float theta;
  int n,r;
  int n_s = input_data_sphere->GetDirections()->GetNumberOfDirections();
  int j;
  float lambda=0.006; //hardcoded here: can make it be an input...0.006 from Max's paper

    

  /*
    data is of size n_s x 1.  
    C is of size rank x 1
    L is of size rank x rank   (Laplace-Beltrami matrix) used if regularizing
    P is of size rank x rank  (ODF transformation matrix) not done here - see DiffusionODFCalculator
    m_ls (least-square transform) is of size rank x n_s
    m_BasisFunction is of size rank x n_s 
  */

  gsl_matrix *data = gsl_matrix_alloc(n_s,1);
  
  gsl_matrix *inter;
  //gsl_matrix_complex *inter2;
  gsl_matrix *inter2;
  gsl_matrix *L;
  //gsl_matrix *P;

  gsl_matrix *C = gsl_matrix_alloc(m_rank,1);

  int p;

  if (!m_sh_basis_is_computed) /* compute SH basis on the first pass, as well as L and P and formula for coeffs, excepting data*/
    {      
      m_sh_basis_is_computed=true;
      m_BasisFunction = ComputeSHMatrix(input_data_sphere->GetDirections(),m_order); 
      
      L = gsl_matrix_alloc(m_rank,m_rank);//(smoothing)
      

      //P = gsl_matrix_alloc(rank,rank);//not used here
      j = 0;

      //note L is sparse: do I need zeros?
      
      /*uncomment if regularizing
      //wrong l: need to go to m_order
      for(int l = 0; l < m_order+1; l+=2) {
	for(int m = -l; m <= l; m++) {
	  
	  //uncomment if regularizing:
	  //gsl_matrix_set(L,j,j,lambda*l*l*(l+1)*(l+1));//this is lamba*L


	  //gsl_matrix_set(P,j,j,2*M_PI*Legendre0(l));
	  j++;
	  
	}
      }	
      */
      
      inter = gsl_matrix_alloc(m_rank,m_rank);
      //inter2 = gsl_matrix_complex_alloc(rank,rank);
      inter2 = gsl_matrix_alloc(m_rank,m_rank);

      m_ls = gsl_matrix_alloc(m_rank,n_s);


      gsl_permutation *permutation=gsl_permutation_alloc(m_rank);

      // Compute LS transform : 

      
      gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,m_BasisFunction,m_BasisFunction,0.0,inter);

      //uncomment if regularizing:
      //gsl_matrix_add(inter,L);

      gsl_linalg_LU_decomp (inter, permutation, &p);

      
      gsl_linalg_LU_invert(inter, permutation, inter2);

      //should inter2 be real? if so, use GSL_REAL(inter2's indiv components in loop); else make new complex matrices for the others; maybe make all complex in the beginning. 
    
      

      gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,inter2,m_BasisFunction,0.0,m_ls);
      
      
      //uncomment if regularizing:
      //gsl_matrix_free(L);
      gsl_matrix_free(inter);
      gsl_matrix_free(inter2);
      gsl_permutation_free(permutation);
    }

      
      
  /* Set data matrix with signal */
  //CImg<float> data(1, ns);

  //cout << "data\n";
  for(n=0;n<n_s;n++) {
    gsl_matrix_set(data,n,0, (double) input_data_sphere->GetValue(x,y,z,n));

    //cout << gsl_matrix_get(data,n,0) << endl; 
    //data(0,n) = m_input_DWI_data->GetValue(x,y,z,n);
  }
  



  /* Compute SH estimation with least-squares (m_ls) transform */

  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m_ls,data,0.0,C);

  gsl_matrix_free(data);

  return C;



}



/*
ImageData *SphericalHarmonics::CalculateSHCoeffs(ImageData *input_data_sphere)
{

	
}
*/

void SphericalHarmonics::SetOrder(int order)
{
  m_order=order;
  m_rank=(m_order+1)*(m_order+2)/2;
}

int SphericalHarmonics::GetRank()
{
  return m_rank;
}

void SphericalHarmonics::CalculateSHCoeffs(ImageData *input_data_sphere,ImageData *output_coeffs, int x,int y, int z)
{


  //call the gsl one
  //this assumes same voxels for each of the imagedatas and that they were init to the right size

  int n;

  gsl_matrix *out_gsl=CalculateSHCoeffs(input_data_sphere,x,y,z);

  //assign ImageData values:
  for (n=0;n<m_rank;n++)
    {
      output_coeffs->SetValue(x,y,z,n,(float) gsl_matrix_get(out_gsl,n,0));
    }

  //clean up:

  gsl_matrix_free(out_gsl);

}


//need to check this
//this outputs only the positive image values
void SphericalHarmonics::CoeffsToImage(ImageData *output_data_sphere, ImageData *coeffs, int x,int y, int z)
{

  /* Project SH representation on the sphere
   */
  int n;
  int n_s=output_data_sphere->GetNSize();

  //if we haven't done anything else with this instance, which is unlikely
  if (!m_sh_basis_is_computed) 
    {      
      m_sh_basis_is_computed=true;
      m_BasisFunction = ComputeSHMatrix(output_data_sphere->GetDirections(),m_order); 
    }

  gsl_matrix *fitted=gsl_matrix_alloc(n_s,1);

  //put coeffs in gsl matrix C:
  gsl_matrix *C=gsl_matrix_alloc(m_rank,1);

  for(n=0;n<m_rank;n++)
    {
     gsl_matrix_set(C,n,0, (double) coeffs->GetValue(x,y,z,n));
     
    }

  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m_BasisFunction,C,0.0,fitted);


  
for(n=0;n<n_s;n++)
    {
      if (gsl_matrix_get(fitted,n,0)>0)
     	{
	 output_data_sphere->SetValue(x,y,z,n,(float) gsl_matrix_get(fitted,n,0));
	}
      else
	{
	  output_data_sphere->SetValue(x,y,z,n,0);
	}
    }  
  

  gsl_matrix_free(fitted);
  gsl_matrix_free(C);

}

//run only if already have computed the coeffs and therefore have m_BasisFunction and m_ls (doesn't check) ; you give it N
void SphericalHarmonics::GetLeverageTerms(float *lev, int x, int y, int z, int N)
	{

	  //H=X(XTX)-1XT; 
	  //(XTX)-1XT is m_ls MxN (m_ls)
	  //X maps coeffs to fitted signal (m_BasisFunction) NxM
	  int n;
	 
	 
	  
	  gsl_matrix *H=gsl_matrix_alloc(N,N);
	
	  
	 
	  

	  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, m_BasisFunction, m_ls, 0.0, H);
	
	  for (n=0;n<N;n++)
	    {
	      //cout << "H value " << gsl_matrix_get(H,n,n) << endl;
	      lev[n]=(float) gsl_matrix_get(H,n,n);
	      
	     }
	  

	}
	
