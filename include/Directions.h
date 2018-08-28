/************************************************************************

   File: Directions.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __Directions_h
#define __Directions_h

#include "misc.h"

/*

DO: add checking that minc file contained doubles (currently, this is assumed).  else could declare as void* and always typecast with methods 

the minc file must have directions in acquisition:direction_x(y,z)

doesn't assume these are always unit vectors: the function _Normalize() can be used to normalize...

DO (maybe): change this to work as does fibretract: never have unassigned directions: add any directions you like, etc.
could do memory allocation within this whenever try to assign a dir that hasn't been allocated

AddDirection: alloc and add 

*/

//typedef struct Vector3D VECTOR3D;
#include <minc_1_rw.h>



class Directions
{
 public: 
  Directions(int n);
  Directions(const char *filename);
  int GetNumberOfDirections();
  float GetXComponent(int n);
  float GetYComponent(int n);
  float GetZComponent(int n);  
  void SetXComponent(int n, float xcomp);
  void SetYComponent(int n, float ycomp);
  void SetZComponent(int n, float zcomp);
  VECTOR3D GetDirection(int n);
  void SetDirection(int n, VECTOR3D vector);
  void WriteToMincFile(minc::minc_1_base& wrt);
  //void WriteToMincFile(char *minc_filename, bool clobber);
  Directions *Transform(const char *xfm_filename);
  Directions *Transform(float rotation[3][3]);
  void IncreaseNumberOfDirections(int additional_directions);
  vtkActor *GetActor();
  vtkActor *GetActor(float red, float green, float blue);
  void Normalize();
  void Normalize(int n);
  void Retain();
  void Release();

 
  protected:
  ~Directions();

 private:
  double *m_directions[3];
  int m_number_of_directions;
  int m_reference_count;
  

};


#endif


//junk:

/*

bool _directions_exist
  bool Directions::DirectionsExist();
*/
