/************************************************************************

   File: VectorImageData.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __VectorImageData_h
#define __VectorImageData_h 

#ifndef __otherlibs_include_h
#include "otherlibs_include.h"
#endif

#ifndef __misc_h
#include "misc.h"
#endif

/*
this is a special class for images of arbitrary numbers of vectors at each voxel

I don't have a file format for this yet: could be just a list of N:x y z .... for each voxel.  I don't believe minc will handle it.

can be plotted with VectorPlotter: takes either vector components (for fixed N) at each voxel, a VectorImageData 

*****************

must allocate VectorData's vector array

3D only

DO: create 3D array here of  VectorDatas.  have AddVector(x,y,z) GetNumberOfVectors(x,y,z), GetVector(

Add: copy constructor, remove/overwrite vector, etc.

 */

typedef struct VectorData
{
  number_of_vectors;
  VECTOR_DATA *vector_array;
  

} VECTOR_DATA;



class VectorImageData
{
 public:
  VectorImageData::VectorImageData();
  VectorImageData::~VectorImageData();
  
  AddVector(x,y,z);
  
 private:


};


#endif
