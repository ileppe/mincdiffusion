/************************************************************************

   File: S2Interpolator.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#include <float.h>
#include "S2Interpolator.h"
#include "otherlibs_include.h"
#include "ImageData.h"
#include "Directions.h"
#include "misc.h"


S2Interpolator::S2Interpolator(ImageData *image_data, Directions *directions, float kernel_FWHM, float kernel_window_radius, INTERPOLATION_TYPE interpolation_type)
{
  m_image_data=image_data;
  m_kernel_FWHM=kernel_FWHM;
  m_kernel_window_radius=kernel_window_radius;
  m_interpolation_type=interpolation_type;
  

  m_lookup_table_is_current=false;
  m_lookup_table_is_allocated=false;

  m_directions=directions;
  m_directions->Retain();
  
  m_assymetric=false;


}



S2Interpolator::~S2Interpolator()
{
  if (m_lookup_table_is_allocated)
    {
      gsl_matrix_float_free(m_lookup_table);    
    }
  m_directions->Release();
}


float S2Interpolator::GetInterpolatedValue(int counter, int voxel_x, int voxel_y, int voxel_z)
{
  float value;
  int i;
  
  
  if (!m_lookup_table_is_current)
    {
      this->CalculateWeights(counter);
    }

  value=0;
 
  for (i=0;i<m_image_data->GetNSize();i++)
    {
      
      value+=m_image_data->GetValue(voxel_x,voxel_y,voxel_z,i)*gsl_matrix_float_get(m_lookup_table,counter,i);	
      
      
    }
  //cout << value << endl; 
  return value;
  
  
}



void S2Interpolator::SetInterpolationType(INTERPOLATION_TYPE interpolation_type)
{
  m_interpolation_type=interpolation_type;
 
  m_lookup_table_is_current=false;


}

void S2Interpolator::SetFWHM(float kernel_FWHM)
{
  m_kernel_FWHM=kernel_FWHM;
  
  m_lookup_table_is_current=false;
  

}

void S2Interpolator::SetWindowRadius(float kernel_window_radius)
{
  m_kernel_window_radius=kernel_window_radius;

  m_lookup_table_is_current=false;

}

void S2Interpolator::CalculateLookupTable()
{
  int i;
  
  if (!m_lookup_table_is_current)
    {
      if (!m_lookup_table_is_allocated)
	{
	  m_lookup_table=gsl_matrix_float_alloc(m_directions->GetNumberOfDirections(),m_image_data->GetNSize());
	  m_lookup_table_is_allocated=true;
	  
	}
      
      for (i=0;i<m_directions->GetNumberOfDirections();i++)
	{
	  
	      this->CalculateWeights(i);
	      
	}
      
      
      
    }
  m_lookup_table_is_current=true;
  
}

//this calculates weights for one direction and stores them in m_lookup_table
void S2Interpolator::CalculateWeights(int counter)
{
  float *weight;
  float sumweights;
  float direction[3];  
  int i,j;
  bool use_direction=true;
  float distance;
  float factor;
  int start;  

  if (m_assymetric)
    {
      start=1;
    }
  else
    {
      start=-1;
    }

  if (!m_lookup_table_is_allocated)
    {
      m_lookup_table=gsl_matrix_float_alloc(m_directions->GetNumberOfDirections(),m_image_data->GetNSize());
      m_lookup_table_is_allocated=true;
      
    }
  
 
  weight=(float *) malloc(m_image_data->GetNSize()*sizeof(float));

  sumweights=0;
  for (i=0;i<m_image_data->GetNSize();i++) 
    {
      weight[i]=0;
      for (j=start;j<2;j+=2)//i.e., + an - all directions if the input is symmetric (sampled only on half the sphere and only + if not
	{
	 switch (m_image_data->GetImageType())
	   {
	   case DWIdata:
	     //if direction in file was 0 0 0, we don't want to use it:
	     
	     //DO: require length of direction to be ~1
	     if ((pow(m_image_data->GetDiffusionParameters()->GetDirections()->GetXComponent(i),2)+pow(m_image_data->GetDiffusionParameters()->GetDirections()->GetYComponent(i),2)+pow(m_image_data->GetDiffusionParameters()->GetDirections()->GetZComponent(i),2))>FLT_EPSILON)
	       {
		 
		 distance=acos(j*m_directions->GetXComponent(counter)*m_image_data->GetDiffusionParameters()->GetDirections()->GetXComponent(i)+j*m_directions->GetYComponent(counter)*m_image_data->GetDiffusionParameters()->GetDirections()->GetYComponent(i)+j*m_directions->GetZComponent(counter)*m_image_data->GetDiffusionParameters()->GetDirections()->GetZComponent(i));
		 use_direction=true;
		 
	       }
	     else use_direction=false;
	     break;
	   case ODFdata:
	     if ((pow(m_image_data->GetDirections()->GetXComponent(i),2)+pow(m_image_data->GetDirections()->GetYComponent(i),2)+pow(m_image_data->GetDirections()->GetZComponent(i),2))>FLT_EPSILON)
	       {
		 
		 distance=acos(j*m_directions->GetXComponent(counter)*m_image_data->GetDirections()->GetXComponent(i)+j*m_directions->GetYComponent(counter)*m_image_data->GetDirections()->GetYComponent(i)+j*m_directions->GetZComponent(counter)*m_image_data->GetDirections()->GetZComponent(i));
		 use_direction=true;
	       }
	     else use_direction=false;
	     break;
	     
	   }
	 if (distance<-1*FLT_EPSILON && use_direction) 
	   {
	     Error("negative distance!\n"); 
	   }
		
     	if (distance<m_kernel_window_radius && use_direction)
	  {
	    switch (m_interpolation_type)
	      {
	      case Linear:
		factor=m_kernel_window_radius-distance;
		break;
	      case Gaussian:
		factor=exp(-1*pow(distance*M_LN2/m_kernel_FWHM,2));
		break;
	      }
	    
	    weight[i]+=factor;
	  }
	}
		  
		  
      sumweights+=weight[i];
		   
	 
      
    }
  

  for (i=0;i<m_image_data->GetNSize();i++)
    {
      if (sumweights>FLT_EPSILON)
	{
	  
	  gsl_matrix_float_set(m_lookup_table,counter,i,weight[i]/sumweights);
	}
      else
	{
	  cout << "window too small";
	  gsl_matrix_float_set(m_lookup_table,counter,i,0);
	  
	  
	}
      
    }
	      
	      

  free(weight);
  

}

void S2Interpolator::SetAssymetric()
{
  m_assymetric=true;
}
  

Directions *S2Interpolator::GetDirections()
{
  return m_directions;
}

