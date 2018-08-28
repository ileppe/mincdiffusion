/************************************************************************

   File: otherlibs_include.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/

#ifndef __otherlibs_include_h
#define __otherlibs_include_h

//this is a hack: we need these, but if they are included after volume_io, we get errors.  But volume_io seems to need to come before the rest (or at something in the rest) of vtk
#include "vtk3DWidget.h"
#include "vtkImagePlaneWidget.h"
#include "vtk.h" 

#endif

//for final version, won't do like this (vtk, gsl) with all the includes: find out what you need and add it as you need it
