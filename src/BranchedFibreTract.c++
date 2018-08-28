/************************************************************************

   File: FibreTract.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __BranchedFibreTract_h
#include "BranchedFibreTract.h" 
#endif

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

BranchedFibreTract::BranchedFibreTract()
{
  m_previous_tract=NULL;
  _number_of_points=0;  
  _tract_connectivity_index=0;
  _point_allocation_step=10;

  int *p=(int *)malloc(1000*sizeof(int));
  free(p);

  _points=(TRACT_POINT *)malloc(_point_allocation_step*sizeof(TRACT_POINT));
  _points_allocated=_point_allocation_step;
  _reference_count=1;


}

BranchedFibreTract::BranchedFibreTract(BranchedFibreTract *tract)
{
  m_previous_tract=tract;
  m_previous_tract->Retain();

  _number_of_points=0;  
  _tract_connectivity_index=0;
  _point_allocation_step=10;

  _points=(TRACT_POINT *)malloc(_point_allocation_step*sizeof(TRACT_POINT));
  
  
  _points_allocated=_point_allocation_step;
  _reference_count=1;
}


int BranchedFibreTract::GetNumberOfPoints()
{
  if (m_previous_tract!=NULL)
    {
      return _number_of_points+m_previous_tract->GetNumberOfPoints();
    }
  else
    {
      return _number_of_points; 
    }
}

POINT3D BranchedFibreTract::GetPoint(int n)
{
  if (m_previous_tract!=NULL)
    {
      if (n<m_previous_tract->GetNumberOfPoints())
	{
	  return m_previous_tract->GetPoint(n);
	}
      
      else
	{
	  return _points[n-m_previous_tract->GetNumberOfPoints()].point;
	}
    }
  else
    {
      return _points[n].point;
    }

}


TRACT_POINT BranchedFibreTract::GetPointInfo(int n)
{
  if (m_previous_tract!=NULL)
    {
      if (n<m_previous_tract->GetNumberOfPoints())
	{
	  return m_previous_tract->GetPointInfo(n);
	}
      else
	{
	  return _points[n-m_previous_tract->GetNumberOfPoints()];
	}
    }
  else
    {
      return _points[n];
    } 
}


void BranchedFibreTract::SetPointInfo(TRACT_POINT tract_point, int n)
{
  if (m_previous_tract!=NULL)
    {
      if (n<m_previous_tract->GetNumberOfPoints())
	{
	  m_previous_tract->SetPointInfo(tract_point,n); 
	}
      else
	{
	  _points[n-m_previous_tract->GetNumberOfPoints()]=tract_point;
	}
    }
  else
    {
      _points[n]=tract_point;
    }
}


bool BranchedFibreTract::GoesThroughVoxel(ImageData *roi, int x, int y, int z)
{
  int i;
  POINT3D point;
  bool goes_through_voxel=false;
  
  for (i=0;i<_number_of_points && goes_through_voxel==false;i++)
    {
      point=this->GetPoint(i);

      
      if ((float_abs((roi->GetXStart()+roi->GetXStep()*x)-point.x)<=float_abs(roi->GetXStep()/2.0)) && (float_abs((roi->GetYStart()+roi->GetYStep()*y)-point.y)<=float_abs(roi->GetYStep()/2.0)) && (float_abs((roi->GetZStart()+roi->GetZStep()*z)-point.z)<=float_abs(roi->GetZStep()/2.0)))
	{
	  goes_through_voxel=true;
	}
      
    }
  if (goes_through_voxel==false && m_previous_tract!=NULL)
    {
      goes_through_voxel=m_previous_tract->GoesThroughVoxel(roi, x, y, z);
    }
  
  return goes_through_voxel;
}

bool BranchedFibreTract::GoesThroughVoxel(IMAGE_GEOM image_geom, int x, int y, int z)
{
  int i;
  POINT3D point;
  bool goes_through_voxel=false;
  
  for (i=0;i<_number_of_points && goes_through_voxel==false;i++)
    {
      point=this->GetPoint(i);

      
      if ((float_abs((image_geom.x_start+image_geom.x_step*x)-point.x)<=float_abs(image_geom.x_step/2.0)) && (float_abs((image_geom.y_start+image_geom.y_step*y)-point.y)<=float_abs(image_geom.y_step/2.0)) && (float_abs((image_geom.z_start+image_geom.z_step*z)-point.z)<=float_abs(image_geom.z_step/2.0)))
	{
	  goes_through_voxel=true;
	}
      
    }

    if (goes_through_voxel==false && m_previous_tract!=NULL)
    {
      goes_through_voxel=m_previous_tract->GoesThroughVoxel(image_geom, x, y, z);
    }
    
  return goes_through_voxel;
}

void BranchedFibreTract::SetIndex(int x,int y, int z, int n)
{
  if (m_previous_tract!=NULL)
    {
      if (n<m_previous_tract->GetNumberOfPoints())
	{
	  m_previous_tract->SetIndex(x,y,z,n);
	}
      else
	{
	  
	  _points[n-m_previous_tract->GetNumberOfPoints()].index.x=x;
	  _points[n-m_previous_tract->GetNumberOfPoints()].index.y=y;
	  _points[n-m_previous_tract->GetNumberOfPoints()].index.z=z;
	}
    }
  else
    {
     _points[n].index.x=x;
     _points[n].index.y=y;
     _points[n].index.z=z; 
    }

}


void BranchedFibreTract::SetIndex(INDEX3D index, int n)
{
  if (m_previous_tract!=NULL)
    {
      if (n<m_previous_tract->GetNumberOfPoints())
	{
	  m_previous_tract->SetIndex(index,n);
	}
      else
	{
	  
	  _points[n-m_previous_tract->GetNumberOfPoints()].index.x=index.x;
	  _points[n-m_previous_tract->GetNumberOfPoints()].index.y=index.y;
	  _points[n-m_previous_tract->GetNumberOfPoints()].index.z=index.z;
	}
    }
  else
    {
     _points[n].index.x=index.x;
     _points[n].index.y=index.y;
     _points[n].index.z=index.z; 
    }
}


INDEX3D BranchedFibreTract::GetIndex(int n)
{
  if (m_previous_tract!=NULL)
    {
      if (n<m_previous_tract->GetNumberOfPoints())
	{
	  return m_previous_tract->GetIndex(n);
	}
      else 
	{
	  return _points[n-m_previous_tract->GetNumberOfPoints()].index;
	}
    }
  else
    {
      return _points[n].index;
    }
}


//void BranchedFibreTract::PrintPoints()

void BranchedFibreTract::Retain()
{
  _reference_count+=1;
  if (m_previous_tract!=NULL)
    {
      m_previous_tract->Retain();
    }
}

BranchedFibreTract::~BranchedFibreTract()
{

  free(_points);
  if (m_previous_tract!=NULL)
    {
      m_previous_tract->Release();
    }
  cout << "delete branched tract\n"; 
}



