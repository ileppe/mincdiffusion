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
#define __ProcessedDTIData_h

#ifndef __otherlibs_include_h
#include "otherlibs_include.h"
#endif

#ifndef __ImageData_h 
#include "ImageData.h"
#endif


/*

combine all processed files into one object: multiple ImageData objects.
 all eigenvectors, eigenvalues, AI, chisquare, trace.

flags: will need for everything, private, used in destructor: free them all

constructors would have to be general: could init with e1 and AI, or the whole kit, or any subsets.

currently, several subsets are supported, e1 is required

DO: (in constructors): checking to see if image sizes match...

add public methods to see if , e.g trace exists, and use these in other code eg Diffusion oDF calculcatior


DOL retain and release image datas!!

*/

class ProcessedDTIData 
{
 public:
  ProcessedDTIData(const char *e1x_filename, const char *e1y_filename, const char *e1z_filename);
  
  ProcessedDTIData(const char *e1x_filename, const char *e1y_filename, const char *e1z_filename, const char *FA_filename);
  
  ProcessedDTIData(const char *e1x_filename, const char *e1y_filename, const char *e1z_filename, 
                   const char *e2x_filename, const char *e2y_filename, const char *e2z_filename, 
                   const char *e3x_filename, const char *e3y_filename, const char *e3z_filename,
                   const char *lambda1_filename, const char *lambda2_filename, const char *lambda3_filename);  
                   
  ProcessedDTIData(const char *e1x_filename, const char *e1y_filename, const char *e1z_filename, 
                   const char *e2x_filename, const char *e2y_filename, const char *e2z_filename, 
                   const char *e3x_filename, const char *e3y_filename, const char *e3z_filename, 
                   const char *lambda1_filename, const char *lambda2_filename, const char *lambda3_filename, 
                   const char *FA_filename, const char *trace_filename);   //not implemented
  
  ~ProcessedDTIData();
  ImageData *Gete1x();
  ImageData *Gete1y();
  ImageData *Gete1z();
  ImageData *Gete2x();
  ImageData *Gete2y();
  ImageData *Gete2z();
  ImageData *Gete3x();
  ImageData *Gete3y();
  ImageData *Gete3z();
  ImageData *Getlambda1();
  ImageData *Getlambda2();
  ImageData *Getlambda3();
  ImageData *Getchisquare();
  ImageData *GetFA();
  ImageData *Gettrace();
  int GetDimensions();
  int GetXSize();
  int GetYSize();
  int GetZSize();
  float GetXStep();
  float GetYStep();
  float GetZStep();
  int GetXStepSign();
  int GetYStepSign();
  int GetZStepSign();
  float GetXStart();
  float GetYStart();
  float GetZStart();
  VECTOR3D GetXDirectionCosine();
  VECTOR3D GetYDirectionCosine();
  VECTOR3D GetZDirectionCosine();

  //DO: add Set methods

 private:
  ImageData *_e1x;
  ImageData *_e1y;
  ImageData *_e1z;
  ImageData *_e2x;
  ImageData *_e2y;
  ImageData *_e2z;
  ImageData *_e3x;
  ImageData *_e3y;
  ImageData *_e3z;
  ImageData *_lambda1;
  ImageData *_lambda2;
  ImageData *_lambda3;
  ImageData *_chisquare;
  ImageData *_FA;
  ImageData *_trace;
  bool _e1x_flag;
  bool _e1y_flag;
  bool _e1z_flag;
  bool _e2x_flag;
  bool _e2y_flag;
  bool _e2z_flag;
  bool _e3x_flag;
  bool _e3y_flag;
  bool _e3z_flag;
  bool _lambda1_flag;
  bool _lambda2_flag;
  bool _lambda3_flag;
  bool _chisquare_flag;
  bool _AI_flag;
  bool _trace_flag;
};


#endif
