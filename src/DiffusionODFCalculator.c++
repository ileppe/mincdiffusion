/************************************************************************

   File: DiffusionODFCalculator.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#include <float.h>
#include <gsl/gsl_blas.h>

#include "DiffusionODFCalculator.h"
#include "otherlibs_include.h"
#include "ImageData.h"
#include "misc.h"
#include "S2Interpolator.h"
#include "Directions.h"
#include "ProcessedDTIData.h"



DiffusionODFCalculator::DiffusionODFCalculator(ImageData *input_DWI_data, Directions *directions_to_calculate)
{


  m_diffusion_ODF_calculation_mode=QBI;
  m_input_DWI_data=input_DWI_data;

 
  m_ODF_image_data=new ImageData(ODFdata,m_input_DWI_data->GetXSize(),m_input_DWI_data->GetYSize(),m_input_DWI_data->GetZSize(),directions_to_calculate->GetNumberOfDirections(),m_input_DWI_data->GetXStart(),m_input_DWI_data->GetYStart(),m_input_DWI_data->GetZStart(),0,m_input_DWI_data->GetXStep(),m_input_DWI_data->GetYStep(),m_input_DWI_data->GetZStep(),1,m_input_DWI_data->GetXDirectionCosine(),m_input_DWI_data->GetYDirectionCosine(),m_input_DWI_data->GetZDirectionCosine());

 
  m_ODF_image_data->SetDirections(directions_to_calculate);
  

  
  m_interpolator_array_is_set=false;

  
  m_angle_step_for_great_circle=0.5*sqrt(2*M_PI/m_input_DWI_data->GetNSize());
  
  m_sphere_harm_constructed=false;

 
}

DiffusionODFCalculator::DiffusionODFCalculator(ProcessedDTIData *input_DTI_data, Directions *directions_to_calculate)
{


  m_diffusion_ODF_calculation_mode=DTI;
  m_input_DTI_data=input_DTI_data;




  m_ODF_image_data=new ImageData(ODFdata,m_input_DTI_data->GetXSize(),m_input_DTI_data->GetYSize(),m_input_DTI_data->GetZSize(),directions_to_calculate->GetNumberOfDirections(),m_input_DTI_data->GetXStart(),m_input_DTI_data->GetYStart(),m_input_DTI_data->GetZStart(),0,m_input_DTI_data->GetXStep(),m_input_DTI_data->GetYStep(),m_input_DTI_data->GetZStep(),1,m_input_DTI_data->GetXDirectionCosine(),m_input_DTI_data->GetYDirectionCosine(),m_input_DTI_data->GetZDirectionCosine());


  m_ODF_image_data->SetDirections(directions_to_calculate);
    
  m_interpolator_array_is_set=false;

    

}




DiffusionODFCalculator::~DiffusionODFCalculator()
{
  int n;
  
  if (m_interpolator_array_is_set)
    {
      for(n=0;n<m_ODF_image_data->GetDirections()->GetNumberOfDirections();n++)
	{
	  delete(m_interpolator_array[n]);
	}
      free(m_interpolator_array);
      
    }

  if (m_sphere_harm_constructed)
    {
      m_no_zero_b->Release();
      gsl_matrix_free(m_basis);
      gsl_matrix_free(m_basis2);
      gsl_matrix_free(m_P);
      delete(m_sphere_harm);

    }

  m_ODF_image_data->Release();
  
}

void DiffusionODFCalculator::SetCalculationMode(DiffusionODFCalculationMode mode)
{
  m_diffusion_ODF_calculation_mode=mode;
}

void DiffusionODFCalculator::InitializeToZero()
{
  
  //preset the _ODF_image_data to be all zero: DO: make this more efficient.  ideally, make this part of ImageData constructor
  //probably need to do this directly through minc: this is way too slow otherwise

  int x,y,z,n;
  

  for(x=0;x<m_ODF_image_data->GetXSize();x++) 
    {
      for(y=0;y<m_ODF_image_data->GetYSize();y++) 
	{
	  for(z=0;z<m_ODF_image_data->GetZSize();z++) 
	    {
	      for(n=0;n<m_ODF_image_data->GetNSize();n++)
		{
		  m_ODF_image_data->SetValue(x,y,z,n,0);
		}
	    }
	}
    }
}


ImageData *DiffusionODFCalculator::GetODFImageData()
{
  return m_ODF_image_data;
  
}

void DiffusionODFCalculator::CalculateODF(int voxel_x, int voxel_y, int voxel_z)
{
  switch (m_diffusion_ODF_calculation_mode)
    {
    case QBI:
      this->CalculateQBIODF(voxel_x, voxel_y, voxel_z);
      break;
    case DTI:
      this->CalculateDTIODF(voxel_x, voxel_y, voxel_z);
      break;
    case QBISH:
      this->CalculateQBISHODF(voxel_x, voxel_y, voxel_z);
      break;
    }
  
}

void DiffusionODFCalculator::CalculateQBISHODF(int voxel_x, int voxel_y, int voxel_z)
{

  //I have no idea whether we can switch dirs and hence basis functions in the middle here

  int order=8;
  
  int j,l,m,n;

  if (!m_sphere_harm_constructed)
    {
      m_sphere_harm=new SphericalHarmonics();
      m_sphere_harm->SetOrder(order);
      m_sphere_harm_constructed=true;


      m_rank=m_sphere_harm->GetRank();

      //DO: put this back or insist no 0 b in input!
      // (or check)
      
      //m_no_zero_b=m_input_DWI_data->CreateDWIODF();//this step is quite slow.
      cout << "warning: assume input has no b0 (QBISH)\n";
      m_no_zero_b=m_input_DWI_data;
      m_no_zero_b->Retain();
  
      m_basis=ComputeSHMatrix(m_no_zero_b->GetDirections(),order);

      //make the basis function for the directions to calculate: DO: check whether can do this
      m_basis2 = ComputeSHMatrix(m_ODF_image_data->GetDirections(),order);

      m_P = gsl_matrix_alloc(m_rank,m_rank);

      j = 0;
      for(l = 0; l <= order; l+=2) {
	for(m = -l; m <= l; m++) {
	  gsl_matrix_set(m_P,j,j,2*M_PI*Legendre0(l));
	  
	  j++;
	}
      }	

      //cout << j << endl;
    }

  //DO: some of these could be allocated just once, but leaving as is for now
  gsl_matrix *coeffs = m_sphere_harm->CalculateSHCoeffs(m_no_zero_b,voxel_x,voxel_y,voxel_z);
  gsl_matrix *coeffs_prime = gsl_matrix_alloc(m_rank,1);


 
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m_P,coeffs,0.0,coeffs_prime);


  gsl_matrix *ODF=gsl_matrix_alloc(m_ODF_image_data->GetDirections()->GetNumberOfDirections(),1);
  

  
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,m_basis2,coeffs_prime,0.0,ODF);
  
  
  for(n=0;n<m_ODF_image_data->GetDirections()->GetNumberOfDirections();n++)
    {
      
      m_ODF_image_data->SetValue(voxel_x,voxel_y,voxel_z,n,(float) gsl_matrix_get(ODF,n,0));
      //cout << m_ODF_image_data->GetValue(voxel_x,voxel_y,voxel_z,n) <<endl;
    }
  
  gsl_matrix_free(coeffs_prime);
  gsl_matrix_free(coeffs);
}


void DiffusionODFCalculator::CalculateQBIODF(int voxel_x, int voxel_y, int voxel_z)
{
  //the first time this is called, the interpolator array is set and the lookup tables are calculated: so for future x,y,z, should be faster
  float theta;
  int n;
  S2Interpolator *interpolator;
  int points_per_great_circle=(int)ceil(2*M_PI/m_angle_step_for_great_circle);
  int i;
  Directions *directions;  //for great circle
  float ODFvalue;
  
  float *basis1 = (float *) malloc(3*sizeof(float));
  float *basis2 = (float *) malloc(3*sizeof(float));
  
  float vector[3];

  if (!m_interpolator_array_is_set)
    {
      
      m_interpolator_array=(S2Interpolator **)malloc(m_ODF_image_data->GetDirections()->GetNumberOfDirections()*sizeof(S2Interpolator *));
    }
  


  for(n=0;n<m_ODF_image_data->GetDirections()->GetNumberOfDirections();n++)
    {
      if (!m_interpolator_array_is_set)
	{
	  //make a basis:
	  vector[0]=m_ODF_image_data->GetDirections()->GetXComponent(n);
	  vector[1]=m_ODF_image_data->GetDirections()->GetYComponent(n);
	  vector[2]=m_ODF_image_data->GetDirections()->GetZComponent(n);
	  
	  MakeTwoOrthogonalVectors(vector, basis1, basis2);

	  
	  directions=new Directions(points_per_great_circle);
	  i=0;
	  
	  for (theta=0;theta<2*M_PI;theta+=m_angle_step_for_great_circle)
	    {
	      directions->SetXComponent(i,basis1[0]*cos(theta)+basis2[0]*sin(theta));
	      directions->SetYComponent(i,basis1[1]*cos(theta)+basis2[1]*sin(theta));
	      directions->SetZComponent(i,basis1[2]*cos(theta)+basis2[2]*sin(theta));
	      
		      
	      i++;
	      
	    }
	  
	  //create interpolator:
	  interpolator=new S2Interpolator(m_input_DWI_data,directions,1.5*sqrt(2*M_PI/m_input_DWI_data->GetNSize()), 1.5*sqrt(2*M_PI/m_input_DWI_data->GetNSize()), Gaussian);
	  m_interpolator_array[n]=interpolator;	  
	  m_interpolator_array[n]->CalculateLookupTable();
	  
	}
      
      
      //interpolator now set, at least for this n: calculate the ODF value in this direction and assign to the ImageData:  must loop through  points_per_great_circle
      ODFvalue=0;
      for (i=0;i<points_per_great_circle;i++)
	{
	  
	  ODFvalue+=m_interpolator_array[n]->GetInterpolatedValue(i,voxel_x,voxel_y,voxel_z);
	}
      m_ODF_image_data->SetValue(voxel_x,voxel_y,voxel_z,n,ODFvalue);
      
    }
  
  
  
  
  m_interpolator_array_is_set=true;
  
  free(basis1);
  free(basis2);
  


}


void DiffusionODFCalculator::CalculateDTIODF(int voxel_x, int voxel_y, int voxel_z)
{
  //note: ODF value set to zero if all eigenvalues are zero
  int n;
  float direction_dot_e1;
  float direction_dot_e2;
  float direction_dot_e3;
  float ODFvalue;
  float one,two,three;
  float lambda1,lambda2,lambda3;
  float e1x,e1y,e1z;
  float e2x,e2y,e2z;
  float e3x,e3y,e3z;
  
  e1x=m_input_DTI_data->Gete1x()->GetValue(voxel_x, voxel_y, voxel_z);
  e1y=m_input_DTI_data->Gete1y()->GetValue(voxel_x, voxel_y, voxel_z);
  e1z=m_input_DTI_data->Gete1z()->GetValue(voxel_x, voxel_y, voxel_z);

  e2x=m_input_DTI_data->Gete2x()->GetValue(voxel_x, voxel_y, voxel_z);
  e2y=m_input_DTI_data->Gete2y()->GetValue(voxel_x, voxel_y, voxel_z);
  e2z=m_input_DTI_data->Gete2z()->GetValue(voxel_x, voxel_y, voxel_z);

  e3x=m_input_DTI_data->Gete3x()->GetValue(voxel_x, voxel_y, voxel_z);
  e3y=m_input_DTI_data->Gete3y()->GetValue(voxel_x, voxel_y, voxel_z);
  e3z=m_input_DTI_data->Gete3z()->GetValue(voxel_x, voxel_y, voxel_z);

  lambda1=m_input_DTI_data->Getlambda1()->GetValue(voxel_x, voxel_y, voxel_z);
  lambda2=m_input_DTI_data->Getlambda2()->GetValue(voxel_x, voxel_y, voxel_z);
  lambda3=m_input_DTI_data->Getlambda3()->GetValue(voxel_x, voxel_y, voxel_z);


  for(n=0;n<m_ODF_image_data->GetDirections()->GetNumberOfDirections();n++)
    {
      
      direction_dot_e1=m_ODF_image_data->GetDirections()->GetXComponent(n)*e1x+m_ODF_image_data->GetDirections()->GetYComponent(n)*e1y+m_ODF_image_data->GetDirections()->GetZComponent(n)*e1z;
      direction_dot_e2=m_ODF_image_data->GetDirections()->GetXComponent(n)*e2x+m_ODF_image_data->GetDirections()->GetYComponent(n)*e2y+m_ODF_image_data->GetDirections()->GetZComponent(n)*e2z;
      direction_dot_e3=m_ODF_image_data->GetDirections()->GetXComponent(n)*e3x+m_ODF_image_data->GetDirections()->GetYComponent(n)*e3y+m_ODF_image_data->GetDirections()->GetZComponent(n)*e3z;
      
    
	//DO: normalize direction_dot_e...

      if (lambda1>FLT_EPSILON)
	{
	  one=pow(direction_dot_e1,2)/lambda1;
         
	}
      else
	{
	  one=FLT_MAX;
	}
      if (lambda2>FLT_EPSILON)
	{
	  two=pow(direction_dot_e2,2)/lambda2;
         
	}
      else
	{
	  two=FLT_MAX;
	}
      if (lambda3>FLT_EPSILON)
	{
	  three=pow(direction_dot_e3,2)/lambda3;
         
	}
      else
	{
	  three=FLT_MAX;
	}
      

      if (lambda1<FLT_EPSILON && lambda2<FLT_EPSILON && lambda3<FLT_EPSILON)
	{
	  ODFvalue=0;
	}
      else if (one+two+three<FLT_EPSILON)
	{
	  ODFvalue=FLT_MAX;
	}
      
      else
	{
	  ODFvalue=sqrt(1/(one+two+three));
	}
      
      
      m_ODF_image_data->SetValue(voxel_x,voxel_y,voxel_z,n,ODFvalue);

      
    }
  
  
}


void DiffusionODFCalculator::NormalizeODF(NormalizationMode normalization_mode,int voxel_x, int voxel_y, int voxel_z)
{
  float norm_factor;
  float scale;
  float value;
  
  
  int n;
  

  switch (normalization_mode)
    {
    case volume:
      //normalize so that all ODFs from this acquisition have unit volume:

      //need to bring radius down to (eg, 1) first, or this blows up:

      switch (m_diffusion_ODF_calculation_mode)
	{
	case QBI:

	  norm_factor=0;
	  scale =0;
	  
	  //get the value of the first direction and set a scale factor (onlt done is the first direction had a nonzero value - alt is to fins the max value and use it...not implemented):
	  scale=m_ODF_image_data->GetValue(voxel_x,voxel_y,voxel_z,0);
	  if (scale>FLT_EPSILON)
	    {
	      for(n=0;n<m_ODF_image_data->GetDirections()->GetNumberOfDirections();n++)
		{ 
		  m_ODF_image_data->SetValue(voxel_x,voxel_y,voxel_z,n,m_ODF_image_data->GetValue(voxel_x,voxel_y,voxel_z,n)/scale);
	       }
	    }
	  
	  for(n=0;n<m_ODF_image_data->GetDirections()->GetNumberOfDirections();n++)
	    {

	      norm_factor+=(4*M_PI*pow(m_ODF_image_data->GetValue(voxel_x,voxel_y,voxel_z,n),3))/(3*m_ODF_image_data->GetDirections()->GetNumberOfDirections());
	    }
	  
	  
	  
	  break;
	  case QBISH:

	  norm_factor=0;
	  scale =0;
	  
	  //get the value of the first direction and set a scale factor (onlt done is the first direction had a nonzero value - alt is to fins the max value and use it...not implemented):
	  scale=m_ODF_image_data->GetValue(voxel_x,voxel_y,voxel_z,0);
	  if (scale>FLT_EPSILON)
	    {
	      for(n=0;n<m_ODF_image_data->GetDirections()->GetNumberOfDirections();n++)
		{ 
		  m_ODF_image_data->SetValue(voxel_x,voxel_y,voxel_z,n,m_ODF_image_data->GetValue(voxel_x,voxel_y,voxel_z,n)/scale);
	       }
	    }
	  
	  for(n=0;n<m_ODF_image_data->GetDirections()->GetNumberOfDirections();n++)
	    {

	      norm_factor+=(4*M_PI*pow(m_ODF_image_data->GetValue(voxel_x,voxel_y,voxel_z,n),3))/(3*m_ODF_image_data->GetDirections()->GetNumberOfDirections());
	    }
	  
	  
	  
	  break;
	case DTI:
	  norm_factor=(4.0/3.0)*M_PI*sqrt(m_input_DTI_data->Getlambda1()->GetValue(voxel_x,voxel_y,voxel_z))*sqrt(m_input_DTI_data->Getlambda2()->GetValue(voxel_x,voxel_y,voxel_z))*sqrt(m_input_DTI_data->Getlambda3()->GetValue(voxel_x,voxel_y,voxel_z));
	  
	    break;
	}
      
      norm_factor=pow(norm_factor, ((float) 1.0)/((float) 3.0));
      break;



    case mean:
      norm_factor=0;
      for(n=0;n<m_ODF_image_data->GetDirections()->GetNumberOfDirections();n++)
	{
	  norm_factor+=m_ODF_image_data->GetValue(voxel_x,voxel_y,voxel_z,n);
	}
      

	
      break;

    }
  
  for(n=0;n<m_ODF_image_data->GetDirections()->GetNumberOfDirections();n++)
    {
      if (norm_factor>FLT_EPSILON)
	{
	  
	  m_ODF_image_data->SetValue(voxel_x,voxel_y,voxel_z,n,m_ODF_image_data->GetValue(voxel_x,voxel_y,voxel_z,n)/norm_factor);
	}
      else 
	{
	  m_ODF_image_data->SetValue(voxel_x,voxel_y,voxel_z,n,0);
	  
	}
    }
  


  
}



