/************************************************************************

   File: Deconvolver.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2007 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/

#include <float.h>
#include <gsl/gsl_sf_erf.h>

#include "Deconvolver.h"
#include "otherlibs_include.h"
#include "ImageData.h"
#include "Directions.h"
#include "misc.h"
#include "S2Interpolator.h"



Deconvolver::Deconvolver()
{
  //default is spherical
  this->SetDeconvolutionMode(Spherical);


  //hardcode order: this is default
  m_order=8;

  m_algorithm=Tournier;
  m_response_sphere_harm_constructed=false;
  m_in_sphere_harm_constructed=false;
  m_response_SH_calculated=false;
  m_response_function=NULL;
  m_deconvolved=NULL;
  m_in_SH=NULL;
  m_response_SH=NULL;
  m_deconvolved_SH=NULL;

}

void Deconvolver::SetOrder(int order)
{
  m_order=order;
}

void Deconvolver::SetDeconvolutionMode(DeconvolutionMode deconvolution_mode)
{
  m_deconvolution_mode=deconvolution_mode;
}


Deconvolver::~Deconvolver()
{
  if (m_deconvolved !=NULL)
    {
      m_deconvolved->Release();
    }
  if (m_in_SH !=NULL)
    {
      m_in_SH->Release();
    }
  if (m_response_SH!=NULL)
    {
      m_response_SH->Release();
    }
  if (m_deconvolved_SH!=NULL)
    {
      m_deconvolved_SH->Release();
    }
  if (m_response_function!=NULL)
    {
      m_response_function->Release();
    }
  if (m_in_sphere_harm_constructed)
    {
      delete m_in_sphere_harm;
    }
  if (m_response_sphere_harm_constructed)
    {
      delete m_response_sphere_harm;
    }
  

  //need to clean up the e1s if they exist, etc.
}


void Deconvolver::SetInput(ImageData *ODF_input)
{

  m_ODF_input=ODF_input;


  //create m_deconvolved ImageData
  
  m_deconvolved=new ImageData(ODFdata,m_ODF_input->GetImageGeom());
  
  m_deconvolved->SetDirections(m_ODF_input->GetDirections());

  //not nec. because SetDirections retains them:
  //m_ODF_input->GetDirections()->Retain();


  //allocate the SH domain ImageDatas: these will have NSize=rank

  int rank=(m_order+1)*(m_order+2)/2;

  m_in_SH= new ImageData(ODFdata,m_ODF_input->GetXSize(),m_ODF_input->GetYSize(),m_ODF_input->GetZSize(),rank,m_ODF_input->GetXStart(),m_ODF_input->GetYStart(),m_ODF_input->GetZStart(),0,m_ODF_input->GetXStep(),m_ODF_input->GetYStep(),m_ODF_input->GetZStep(),1,m_ODF_input->GetXDirectionCosine(),m_ODF_input->GetYDirectionCosine(),m_ODF_input->GetZDirectionCosine());


  //these starts and steps aren't right - not used.
  m_response_SH = new ImageData(ODFdata,1,1,1,rank,m_ODF_input->GetXStart(),m_ODF_input->GetYStart(),m_ODF_input->GetZStart(),0,m_ODF_input->GetXStep(),m_ODF_input->GetYStep(),m_ODF_input->GetZStep(),1,m_ODF_input->GetXDirectionCosine(),m_ODF_input->GetYDirectionCosine(),m_ODF_input->GetZDirectionCosine());

  m_deconvolved_SH = new ImageData(ODFdata,m_ODF_input->GetXSize(),m_ODF_input->GetYSize(),m_ODF_input->GetZSize(),rank,m_ODF_input->GetXStart(),m_ODF_input->GetYStart(),m_ODF_input->GetZStart(),0,m_ODF_input->GetXStep(),m_ODF_input->GetYStep(),m_ODF_input->GetZStep(),1,m_ODF_input->GetXDirectionCosine(),m_ODF_input->GetYDirectionCosine(),m_ODF_input->GetZDirectionCosine());

  //DO: redo this with a find max function: currently ASSUMES last b value is highest (and that all const ex. zero!!!)

  
  if (0)
    {
      m_b=1000;
    }
  else if (1) //assume always constant b parameter:
    {
      m_b=m_ODF_input->GetDiffusionParameters()->GetConstb();
    }
  else //never
    {
      m_b=m_ODF_input->GetDiffusionParameters()->GetBValue(m_ODF_input->GetNSize()-1);
    }

}


//the response_function ImageData and mask must match input sampling
//this uses the response function minimum as the fibre dir, but DO: should do SH fit first to make minimum detection more robust. prob. use order 2 so really smooth: expect simple shape. (?)
//I recommend using e1 as the reponse min instead of using this method: see below
void Deconvolver::SetResponseFunctionSingleVoxel(ImageData *response_function, ImageData *mask_response_voxel)
{
  //make new ImageData for response: one voxel (which could be average of many input voxels in future);starts NOT adjusted!!!
  m_response_function=new ImageData(ODFdata,1,1,1,response_function->GetNSize(),0,0,0,0,1,1,1,1,response_function->GetXDirectionCosine(),response_function->GetYDirectionCosine(),response_function->GetZDirectionCosine());

  m_response_function->SetDirections(response_function->GetDirections());


  int x,y,z,n;
  int xf,yf,zf;
  int cont=true;
  int index;
  int i;
  float alpha, beta;

  //find the nonzero voxel:

  for (x=0;x<response_function->GetXSize()&&cont;x++)
    {
      for (y=0;y<response_function->GetYSize()&&cont;y++)
	{
	  for (z=0;z<response_function->GetZSize()&&cont;z++)
	    {
	      if (mask_response_voxel->GetValue(x,y,z)>FLT_EPSILON)
		{
		  xf=x;
		  yf=y;
		  zf=z;
		  cont=false;
		}
	    }
	}
    }

  //now set the data:

  for (n=0;n<response_function->GetNSize();n++)
    {
      m_response_function->SetValue(0,0,0,n,response_function->GetValue(xf,yf,zf,n));
    }

  //DO: blur first before getting minimum: might make it more accurate
  //ImageData *blurred_reponse_function;

  //find the minimum:

  m_response_function->GetMinimumNonZeroValue(xf,yf,zf,index);

  //now rotate directions so that the minimum one lies along 0,0,1: can use same ImageData and just swap its Directions.  note that response is assumed to be axially symmetric: fibre along 0,0,1 is all that is assumed

  float rotation[3][3];



  VECTOR3D vector=m_response_function->GetDirections()->GetDirection(index);


  this->RotateResponse(vector);
 
  //debug:
  //m_response_function->WriteMincFileNoRange("test.mnc",true,"debug");



}

void Deconvolver::SetResponseFunctionSingleVoxel(ImageData *response_function, ImageData *mask_response_voxel, ImageData *e1x, ImageData *e1y, ImageData *e1z)
{


  m_response_function=new ImageData(ODFdata,1,1,1,response_function->GetNSize(),0,0,0,0,1,1,1,1,response_function->GetXDirectionCosine(),response_function->GetYDirectionCosine(),response_function->GetZDirectionCosine());

  m_response_function->SetDirections(response_function->GetDirections());

  int n;
  int x,y,z;
  int xf,yf,zf;
  int cont=true;
  int index;
  int i;



  //find the nonzero voxel:

  for (x=0;x<response_function->GetXSize()&&cont;x++)
    {
      for (y=0;y<response_function->GetYSize()&&cont;y++)
	{
	  for (z=0;z<response_function->GetZSize()&&cont;z++)
	    {
	      if (mask_response_voxel->GetValue(x,y,z)>FLT_EPSILON)
		{
		  xf=x;
		  yf=y;
		  zf=z;
		  cont=false;
		}
	    }
	}
    }


  //now set the data:

  for (n=0;n<response_function->GetNSize();n++)
    {
      m_response_function->SetValue(0,0,0,n,response_function->GetValue(xf,yf,zf,n));
    }


  //now rotate directions so that the PDD one lies along ???: can use same ImageData and just swap its Directions.  note that response is assumed to be axially symmetric: fibre along 0,0,1 is all that is assumed

  //could do instead...: just save the rotation and apply it to the result at the end - inverse of this:


  VECTOR3D vector;
  vector.x=e1x->GetValue(xf,yf,zf);
  vector.y=e1y->GetValue(xf,yf,zf);
  vector.z=e1z->GetValue(xf,yf,zf);


  this->RotateResponse(vector);

}

void Deconvolver::SetResponseFunction(ImageData *response_function)
{
  m_response_function=response_function;
}

void Deconvolver::Deconvolve(int x,int y,int z)
{

  Deconv2Sphere(x,y,z);



}


void Deconvolver::Deconv2Sphere(int x,int y,int z)
{
  
  

  int n,l;
  int rank=(m_order+1)*(m_order+2)/2;
  float quotient;

  if (!m_in_sphere_harm_constructed)
    {
      m_in_sphere_harm=new SphericalHarmonics();
      m_in_sphere_harm->SetOrder(m_order);
      m_in_sphere_harm_constructed=true;
    }

  m_in_sphere_harm->CalculateSHCoeffs(m_ODF_input,m_in_SH,x,y,z);

  if (!m_response_SH_calculated && m_algorithm==Tournier)
    {
      m_response_sphere_harm=new SphericalHarmonics();
      m_response_sphere_harm->SetOrder(m_order);
      m_response_sphere_harm_constructed=true;
      m_response_sphere_harm->CalculateSHCoeffs(m_response_function,m_response_SH,0,0,0);
      /*
      for (n=0;n<m_response_function->GetNSize();n++)
	{
	  cout << m_response_function->GetValue(0,0,0,n) << endl;
	}
      for (n=0;n<m_sphere_harm->GetRank();n++)
	{
	  
	  cout << m_response_SH->GetValue(0,0,0,n) << endl;
	}
      */
      m_response_SH_calculated=true;
    }
  //now divide the values to get the SH rep of the output (can we do this when they had different sets of directions???)
  
  //test: if the rotational harmonic formulation is correct and Tournier's simplication of it is correct, we want to set all m!=0 indices for the response to zero.  All m=0 indices of fODF coeffs will be set to the quotient of the SH coeff for the signal and the SH coeff for the response.  this means we use the n(l,m)=(l^2+l+2)/2 indices.  what about the zero values?  should fODF coeff go to inf?

  bool m0=false;



  for(n=0;n<rank;n++)
    {
     
       switch(m_algorithm)
	 {
	 case FORECAST:
      
	
	   quotient=Scale(thisl(n))*exp(m_b*m_lambda2)/(2*M_PI*m_response_So*A(thisl(n))*m_b*(m_lambda1-m_lambda2))*m_in_SH->GetValue(x,y,z,n);
	   break;
      
      
	 case Tournier:
	   m0=false;
	   for(l=0;l<=m_order&&!m0;l+=2)
	     {
	  
	       if ((n+1)==((l*l+l+2)/2))
		 {
		   m0=true;
	      
		 }
	     }
      
	   if (m0 && float_abs(m_response_SH->GetValue(0,0,0,n))>FLT_EPSILON)
	     {
	       //includes conversion from SH(l) to RH(l).

	       quotient=Scale(sqrt((2*thisl(n)+1)/4*M_PI)*thisl(n))*this->DeltaCoeff(thisl(n))*m_in_SH->GetValue(x,y,z,n)/m_response_SH->GetValue(0,0,0,n);

	  



	     }
	   else if (float_abs(m_response_SH->GetValue(0,0,0,n))<=FLT_EPSILON)
	     {
	       quotient=0; //these are not inf: just undefined (i.e. can't gain back frequencies lost to convolution??)
	  
	     }
	   else 
	     {
	       quotient=0;
	     }
	   break;//Tournier
	 }

      m_deconvolved_SH->SetValue(x,y,z,n,quotient);   

    } 


  m_in_sphere_harm->CoeffsToImage(m_deconvolved,m_deconvolved_SH,x,y,z);




}


ImageData *Deconvolver::GetOutput()
{

  if (0)
    {
  //not used anymore.  this is a hack and hints that the original rotation of response is incorrect.  need to rotate the output dirs to align image


  //R90z
  m_rotation[0][0]=0;
  m_rotation[0][1]=1;
  m_rotation[0][2]=0;
  m_rotation[1][0]=-1;
  m_rotation[1][1]=0;
  m_rotation[1][2]=0;
  m_rotation[2][0]=0;
  m_rotation[2][1]=0;
  m_rotation[2][2]=1;
  

  //R90y
  /*
      m_rotation[0][0]=0;
      m_rotation[0][1]=0;
      m_rotation[0][2]=-1;
      m_rotation[1][0]=0;
      m_rotation[1][1]=1;
      m_rotation[1][2]=0;
      m_rotation[2][0]=1;
      m_rotation[2][1]=0;
      m_rotation[2][2]=0;
  */

  Directions *dirs;
  Directions *trans_dirs;
  dirs=m_deconvolved->GetDirections();


  trans_dirs=dirs->Transform(m_rotation);


  m_deconvolved->GetDirections()->Release();//I think should be done
  m_deconvolved->SetDirections(trans_dirs);
    }
    
  m_deconvolved->Retain();//retain it here so user can delete this...
  return m_deconvolved;
  //test:
  //return m_response_function; 
}

void Deconvolver::RotateResponse(VECTOR3D vector)
{
  float alpha, beta;
 

  if (0)//test no rotation
    {
	  m_rotation[0][0]=1;
	  m_rotation[0][1]=0;
	  m_rotation[0][2]=0;
	  m_rotation[1][0]=0;
	  m_rotation[1][1]=1;
	  m_rotation[1][2]=0;
	  m_rotation[2][0]=0;
	  m_rotation[2][1]=0;
	  m_rotation[2][2]=1;
    }




  //note: takes e1 to origin as expected.  this is ill-posed though - needs revision for already along z case.

  if (1)
    {



      alpha=atan(vector.x/vector.z);
      beta=asin(vector.y/sqrt(vector.x*vector.x+vector.y*vector.y+vector.z*vector.z));

      m_rotation[0][0]=cos(alpha);
      m_rotation[0][1]=0.0;
      m_rotation[0][2]=sin(alpha);
      m_rotation[1][0]=-sin(alpha)*sin(beta);
      m_rotation[1][1]=cos(beta);
      m_rotation[1][2]=cos(alpha)*sin(beta);
      m_rotation[2][0]=-sin(alpha)*cos(beta);
      m_rotation[2][1]=-sin(beta);
      m_rotation[2][2]=cos(alpha)*cos(beta);


      //another rotation around z? shouldn't be nec.


    }

  //now apply rotation to the response function:

  Directions *dirs;
  Directions *trans_dirs;
  dirs=m_response_function->GetDirections();

  //cout << "orig: " << dirs->GetDirection(0).x << endl;

  trans_dirs=dirs->Transform(m_rotation);

  //cout << "after: " << trans_dirs->GetDirection(0).x << endl;

  m_response_function->GetDirections()->Release();//I think should be done
  m_response_function->SetDirections(trans_dirs);

  //cout << "after: " << m_response_function->GetDirections()->GetDirection(0).x << endl;

}

float Deconvolver::DeltaCoeff(int l)
{
  
  return sqrt((2*l+1)/4*M_PI)*Legendre0(l);

}

float Deconvolver::A(int l)
{
  //setting a to 1: not sure what it should be???
  double a=1.0;
  switch(l)
    {
    case 0:
      return (float) 1.0/pow(a,0.5)*sqrt(M_PI)*gsl_sf_erf(sqrt(a));
      break;
    case 2:
      return (float) -1.0/(4*pow(a,1.5))*(6*sqrt(a)*exp(-a)+(-3+2*a)*sqrt(M_PI)*gsl_sf_erf(sqrt(a)));
      break;
    case 4:
      return (float) 1.0/(32*pow(a,2.5))*(-10*sqrt(a)*(21+2*a)*exp(-a)+3*(35-20*a+4*pow(a,2))*sqrt(M_PI)*gsl_sf_erf(sqrt(a)));
      break;
    case 6:
      return (float) -1.0/(128*pow(a,3.5))*(42*sqrt(a)*(165+20*a+4*pow(a,2))*exp(-a)+5*(-693+378*a-84*pow(a,2)+8*pow(a,3))*sqrt(M_PI)*gsl_sf_erf(sqrt(a)));
      break;
    case 8:
      return (float) 1.0/(2048*pow(a,4.5))*(-6*sqrt(a)*(225225+30030*a+7700*pow(a,2)+248*pow(a,3))*exp(-a)+35*(19305-10296*a+2376*pow(a,2)-288*pow(a,3)+16*pow(a,4))*sqrt(M_PI)*gsl_sf_erf(sqrt(a)));
      break;
    }
}

int Deconvolver::thisl(int n)
{

  if (n==0)
    {
      return 0;
    }
  else if (0<n && n<6)
    {
      return 2;
    }
  else if (6<=n && n<15)
    {
      return 4;
    }
  else if (15<=n && n<28)
    {
      return 6;
    }
  else if (28<=n && n<45)
    {
      return 8;
    }

}


float Deconvolver::Scale(int l)
{
  //filter for coeffs:
  //these values are empirical and wil need adjusting given the dataset (and algorithm - for FORECAST right now, not the same as Tournier's {1,1,1,0.8,0.1}
  float scale=1.0;


  switch (l)
    {
    case 4:
      scale=0.2;
      break;
    case 6:
      scale=0;
      break;
    case 8:
      scale=0;
      break;
    }
  return scale;
}

void Deconvolver::SetDeconvolutionAlgorithm(DeconvolutionAlgorithm algorithm)
{
  m_algorithm=algorithm;
}

void Deconvolver::GetEigenvaluesAndSo(float &lambda1, float &lambda2, float &So){
  lambda1=m_lambda1;
  lambda2=m_lambda2;
  So=m_response_So;
}

void Deconvolver::SetLambdas(float lambda1, float lambda2, float So)
{
  m_lambda1=lambda1;
  m_lambda2=lambda2;
  m_response_So=So;
}

void Deconvolver::SetLambdas(ImageData *lambda1,ImageData *lambda2,ImageData *mask_response_voxel)
{

  int x,y,z;
  int count=0;

  m_lambda1=0;
  m_lambda2=0;
  m_response_So=0;
  //find the nonzero voxel and average the lambdas:

  for (x=0;x<lambda1->GetXSize();x++)
    {
      for (y=0;y<lambda1->GetYSize();y++)
	{
	  for (z=0;z<lambda1->GetZSize();z++)
	    {
	      if (mask_response_voxel->GetValue(x,y,z)>FLT_EPSILON)
		{
		  m_lambda1+=lambda1->GetValue(x,y,z);
		  m_lambda2+=lambda2->GetValue(x,y,z);
		  //DO: fix this: currently ASSUMES first frame is b=0
		  m_response_So+=m_response_function->GetValue(x,y,z,0);
		  count+=1;
		}
	    }
	}
    }
    
  m_lambda1=m_lambda1/count;
  m_lambda2=m_lambda2/count;
  m_response_So=m_response_So/count;

}
void Deconvolver::SetLambdasDefault()
{
  /*Set default values for lambda and response_So, these are based on an average of 25 voxels in a normal subject,
  central corpus callosum, on the TIM Trio with b=1000*/
  m_lambda1=0.001814152052;
  m_lambda2=0.0002350908912;
  m_response_So=156.5972435;			

}
