/************************************************************************

   File: FibreTract.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __BranchedFibreTract_h
#define __BranchedFibreTract_h 

#ifndef __otherlibs_include_h
#include "otherlibs_include.h"
#endif

#ifndef __misc_h
#include "misc.h"
#endif

#ifndef __ImageData_h
#include "ImageData.h"
#endif

#ifndef __FibreTract_h
#include "FibreTract.h"
#endif



/*
here, the number of points will be the total of all sub tracts, and GetPoint (etc) will be reimplemented.  the fact that the tracts don't actually have other copies should be hidden.  The ference counting should be taken care of properly.  The constructor should either take another tract (this tract would append to it), or nothing.  the structure itself needs to have a pointer to the previous tract.  There will always be a m_previous_tract, but it could be empty.  A BranchedFibreTract is a FibreTract.  BranchedFibreTracts must be declared as such.  I'm not sure how private var.s work.  I think those that are shared should be protected.  ReAllocate_Tract, etc, things that assume the TRACT_POINT * needs to be big, need to be redone.  -or-, the private m_number_of_points could be different from GetNumberOfPoints?  destructor would release the previous tracts, so they can be kept to themselves if desired, or be part of other BranchedFibreTracts

!NOTE: we don't really want to write this plus the previous segments to file: we just want to write the unique part.  hence, FibreTractSet, when it writes to file, should ask what sort of tract this is.This assumes that the other segments of this tract are guaranteed to be in the set.  that's why it shuld be an potion, perhaps the user specifies this to the Write... function: write all or write unique

commented functions need to be reimplemented but aren't urgent


*/

class BranchedFibreTract : public FibreTract
{
 public:
  BranchedFibreTract();
  BranchedFibreTract(BranchedFibreTract *tract);
  
  int GetNumberOfPoints();
  POINT3D GetPoint(int n);
  TRACT_POINT GetPointInfo(int n);
  void SetPointInfo(TRACT_POINT tract_point, int n);


  bool GoesThroughVoxel(ImageData *roi, int x, int y, int z);

  bool GoesThroughVoxel(IMAGE_GEOM image_geom, int x, int y, int z);


  void SetIndex(int x,int y, int z, int n);

  void SetIndex(INDEX3D index, int n);
  INDEX3D GetIndex(int n);
  void Retain();
  //void BranchedFibreTract::PrintPoints();


 protected:
  ~BranchedFibreTract();

 private:
  *m_previous_tract;


};

#endif
