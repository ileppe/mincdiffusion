/************************************************************************

   File: DiffusionParameters.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __DiffusionParameters_h
#define __DiffusionParameters_h

#ifndef __otherlibs_include_h
#include "otherlibs_include.h"
#endif

#ifndef __Directions_h
#include "Directions.h"
#endif


/*

currently DiffusionParameters is in minc space and a DiffusionParameters in an ImageData is in the correct corresponding minc space

includes the bvalues

can include bmatrix format as well (calculate it, return it, or convert vice versa)

will enable read-in of many file formats for specifying these things

DO: add checking that minc file contained doubles (currently, this is assumed)

*/

class DiffusionParameters
{
 public: 
  DiffusionParameters(int n);
  DiffusionParameters(const char *minc_filename);
  ~DiffusionParameters();
  Directions *GetDirections();
  double *GetBValues();
  float GetBValue(int n);
  void SetBValue(int n, float value);
  int GetConstb();
  void SetConstb(int b);

  int GetNumberOfDWIs()
  {
    return _number_of_DWIs;
  }
  
  
  void WriteToMincFile(minc::minc_1_writer& wrt);
  void SetDirections(Directions *directions);
  
  void Retain();
  void Release();

 private:
  Directions *_directions;
  double *_bvalues;
  int _number_of_DWIs;
  bool _directions_are_set;
  int _reference_count;
  int m_const_b;
  bool m_const_b_set;

};


#endif

//junk:

/*

  bool DiffusionParameters::BValuesExist();
  bool DiffusionParameters::DirectionsExist();
  bool DiffusionParameters::BValuesAndDirectionsExist();



*/
