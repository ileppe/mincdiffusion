/************************************************************************

   File: ProcessedDTIData.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __ProcessedDTIData_h
#include "ProcessedDTIData.h"
#endif

#ifndef __otherlibs_include_h
#include "otherlibs_include.h"
#endif

#ifndef __ImageData_h 
#include "ImageData.h"
#endif

ProcessedDTIData::ProcessedDTIData(const char *e1x_filename, const char *e1y_filename, const char *e1z_filename)
{
  _e1x=new ImageData(e1x_filename);
  _e1x_flag=true;

  _e1y=new ImageData(e1y_filename);
  _e1y_flag=true;
  
  _e1z=new ImageData(e1z_filename);
  _e1z_flag=true;
}


ProcessedDTIData::ProcessedDTIData(const char *e1x_filename,const char *e1y_filename,const char *e1z_filename, 
                                   const char *e2x_filename,const char *e2y_filename,const char *e2z_filename, 
                                   const char *e3x_filename,const char *e3y_filename,const char *e3z_filename, 
                                   const char *lambda1_filename,const char *lambda2_filename, const char *lambda3_filename)
{
  _e1x=new ImageData(e1x_filename);
  _e1x_flag=true;

  _e1y=new ImageData(e1y_filename);
  _e1y_flag=true;
  
  _e1z=new ImageData(e1z_filename);
  _e1z_flag=true;

  _e2x=new ImageData(e2x_filename);
  _e2x_flag=true;

  _e2y=new ImageData(e2y_filename);
  _e2y_flag=true;
  
  _e2z=new ImageData(e2z_filename);
  _e2z_flag=true;

  _e3x=new ImageData(e3x_filename);
  _e3x_flag=true;

  _e3y=new ImageData(e3y_filename);
  _e3y_flag=true;
  
  _e3z=new ImageData(e3z_filename);
  _e3z_flag=true;

  _lambda1=new ImageData(lambda1_filename);
  _lambda1_flag=true;
  
  _lambda2=new ImageData(lambda2_filename);
  _lambda2_flag=true;
  
  _lambda3=new ImageData(lambda3_filename);
  _lambda3_flag=true;
  
}

ProcessedDTIData::~ProcessedDTIData()
{
  //release all image datas, given that they have been retained!
}

//DO: checks that these exist before returning...

ImageData *ProcessedDTIData::Gete1x()
{
  return _e1x;
  
}

ImageData *ProcessedDTIData::Gete1y()
{
  return _e1y;
  
}

ImageData *ProcessedDTIData::Gete1z()
{
  return _e1z;
  
}

ImageData *ProcessedDTIData::Gete2x()
{
  return _e2x;
  
}

ImageData *ProcessedDTIData::Gete2y()
{
  return _e2y;
  
}

ImageData *ProcessedDTIData::Gete2z()
{
  return _e2z;
  
}

ImageData *ProcessedDTIData::Gete3x()
{
  return _e3x;
  
}

ImageData *ProcessedDTIData::Gete3y()
{
  return _e3y;
  
}

ImageData *ProcessedDTIData::Gete3z()
{
  return _e3z;
  
}

ImageData *ProcessedDTIData::Getlambda1()
{
  return _lambda1;
  
}

ImageData *ProcessedDTIData::Getlambda2()
{
  return _lambda2;
  
}

ImageData *ProcessedDTIData::Getlambda3()
{
return _lambda3;
}

ImageData *ProcessedDTIData::Getchisquare()
{
  return _chisquare;
  
}

ImageData *ProcessedDTIData::GetFA()
{
  return _FA;
  

}

ImageData *ProcessedDTIData::Gettrace()
{
  return _trace;
  
}

int ProcessedDTIData::GetDimensions()
{
  return _e1x->GetDimensions();
  
}

int ProcessedDTIData::GetXSize()
{
  return _e1x->GetXSize();
  
}

int ProcessedDTIData::GetYSize()
{
  return _e1x->GetYSize();
}

int ProcessedDTIData::GetZSize()
{
  return _e1x->GetZSize();
}

float ProcessedDTIData::GetXStep()
{
  return _e1x->GetXStep();
}

float ProcessedDTIData::GetYStep()
{
  return _e1x->GetYStep();
}

float ProcessedDTIData::GetZStep()
{
  return _e1x->GetZStep();
}

float ProcessedDTIData::GetXStart()
{
  return _e1x->GetXStart();
  
}

float ProcessedDTIData::GetYStart()
{
  return _e1x->GetYStart();
}

float ProcessedDTIData::GetZStart()
{
  return _e1x->GetZStart();
}

VECTOR3D ProcessedDTIData::GetXDirectionCosine()
{
  return _e1x->GetXDirectionCosine();
  
}

VECTOR3D ProcessedDTIData::GetYDirectionCosine()
{
  return _e1x->GetYDirectionCosine();
}

VECTOR3D ProcessedDTIData::GetZDirectionCosine()
{
  return _e1x->GetZDirectionCosine();
}

