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
#ifndef __FibreTract_h
#define __FibreTract_h 

#include <vector>
#ifndef __otherlibs_include_h
#include "otherlibs_include.h"
#endif

#ifndef __misc_h
#include "misc.h"
#endif

#ifndef __ImageData_h
#include "ImageData.h"
#endif



/*
a fibre tract is an array of points accompanied by a tract connectivity index (sometimes I call this psi). _points[i] is the three-tuple

allocate each new point as you go and use realloc on the array of pointers to these.

note vtk could be used here: vtkPoints, ... but I'm not doing that

removing points is not implemented at this time and would require data structure modification (could just use a flag for each point, and rewrite when GetPoints() is called

the tract is allocated 10 points at a time

could add more methods, e.g. smoothing, if required

no tract addition so far: not sure when we would want it.  the tract connectivity indices should not be related for the forward and backward connections

not sure we need GetPoints()

should have a min_tangent_diffusion_ODF_value associated with this as well: "strength" of connection/AI directional/ not sure what to call it.  perhaps tract AI: tract index of diffusion anisotropy; index of anisotropy parallel/tangent to tract;

have both mean and min diffusion ODF value tangent to tract (MDOTT)

TractConnectivityIndex=minimum diffusion ODF value tangent to tract

DO: enable error output if try to access points when they don't exist

DO: add scalar(s) at each point of the tract

maybe add index to the point (so searching not nec)

could add flags for whether certain scalars exist;

maybe add changing of points after the fact

could have an array of indices (INDEX3D) int that are associated with a given point, to speed the minc file writeout

add: computemindott, etc

DOTT=DiffusionODFValueTangentToTract


currently, setDOTT needs to be called for only the last point, because it sets cumulative min dott as well: could use gsl min functions to quickly get the cum min DOTT for any point, etc...

should get rif of ...DOTT and just use OTT and have a type for the ODF: D(iffusion) or F(ibre)

could just call these scalars (most general) and have arbitrarily many of these

as is, tract point may not have an index (or DOTT, etc.) associated with it, so be careful using it later.  could add AddPoint() that adds these at the same time (at least index)

DO: change goesthrough roi / goes through voxel to use indices if available. maybe have a flag saying whether these exist


note: FibreTract(tract1,tract2) is for the purpose of multiple rois.  the indices at each point, including cumulative indices, are relative to the first roi given (the seed)  TractConnectivityIndex is assumed not to be used

 */

//typedef enum addition_mode {beginning_reverse, beginning_forward, end_reverse, end_forward} ADDITION_MODE;

struct tract_point
{
  //not sure how to initialize/allocate.  now, scalars are those listed 
  Point3D point;
  float cumulative_CI;
  float cumulative_min_DOTT;
  float DOTT;
  float OTT;
  float cumulative_min_OTT;
  Index3D index;
  
  tract_point()
  {}
  
  tract_point(const Point3D& p):
    point(p) //VF: zero all the elements just in case?
  {}
};

typedef tract_point TRACT_POINT;

class FibreTract
{
 public:
  FibreTract();
  
  FibreTract(FibreTract *tract1, FibreTract *tract2);
  
  FibreTract(FibreTract *tract);
  
  int GetNumberOfPoints(void) const
  {
    if (m_previous_tract!=NULL)
    {
        return _points.size()+m_previous_tract->GetNumberOfPoints();
    } else  {
        return _points.size(); 
    }
  }
  
  float GetTractLength(void) const /*IL*/
  {
    if (m_previous_tract!=NULL)
    {
        return _tract_length+m_previous_tract->GetTractLength();
    } else {
        return _tract_length; 
    }
  }
  
  const Point3D &GetPoint(int n) const 
  {
    if (m_previous_tract!=NULL)
    {
      if (n<m_previous_tract->GetNumberOfPoints())
        return m_previous_tract->GetPoint(n);
      else
        return _points[n-m_previous_tract->GetNumberOfPoints()].point;
    } else {
        return _points[n].point;
    }
  }
  
  const tract_point& GetPointInfo(int n) const
  {
    if (m_previous_tract!=NULL)
    {
        if (n<m_previous_tract->GetNumberOfPoints())
          return m_previous_tract->GetPointInfo(n);
        else
          return _points[n-m_previous_tract->GetNumberOfPoints()];
    }  else 
        return _points[n];
  }  
  
  void SetPointInfo(tract_point tract_point, int n)
  {
    if (m_previous_tract!=NULL)
    {
        if (n<m_previous_tract->GetNumberOfPoints())
        {
          m_previous_tract->SetPointInfo(tract_point,n); 
        }  else  {
          _points[n-m_previous_tract->GetNumberOfPoints()]=tract_point;
        }
    } else  {
        _points[n]=tract_point;
    }
  }
  
  //float **FibreTract::GetPoints();
  void AddPoint(float x, float y, float z);
  void AddPoint(const POINT3D& point);


  void SetTractConnectivityIndex(float tract_connectivity_index);
  float GetTractConnectivityIndex();
  void SetCumulativeCI(float cumulative_CI, int n);
  float GetCumulativeCI(int n);

  void SetDOTT(float DOTT, int n);
  float GetDOTT(int n);
  void SetMinDOTT(float min_DOTT);
  float GetMinDOTT();
  void SetMeanDOTT(float mean_DOTT);
  float GetMeanDOTT();
  void SetCumulativeMinOTT(float CMOTT, int n);
  float GetCumulativeMinDOTT(int n);

  void SetOTT(float OTT, int n);
  float GetOTT(int n);
  float GetBlurredOTT(int n);

  void SetMinOTT(float min_OTT);
  float GetMinOTT();
  float GetCumulativeMinOTT(int n);

  //void FibreTract::AddFibreTract(FibreTract *fibre_tract_to_add, ADDITION_MODE addition_mode); //for concentenating two tracts 
  bool GoesThroughVoxel(ImageData *roi, int x, int y, int z);
  bool GoesThroughVoxel(IMAGE_GEOM image_geom, int x, int y, int z);
  bool GoesThroughROI(ImageData *roi);

  void SetIndex(int x,int y, int z, int n)
  {
     SetIndex(INDEX3D(x,y,z),n);
  }

  void SetIndex(const INDEX3D& index, int n)
  {
    if (m_previous_tract!=NULL)
    {
      if (n<m_previous_tract->GetNumberOfPoints())
      {
        m_previous_tract->SetIndex(index,n);
      } else {
        _points[n-m_previous_tract->GetNumberOfPoints()].index.x=index.x;
        _points[n-m_previous_tract->GetNumberOfPoints()].index.y=index.y;
        _points[n-m_previous_tract->GetNumberOfPoints()].index.z=index.z;
      }
    } else {
      _points[n].index=index;
    }
  }
  
  const INDEX3D& GetIndex(int n) const
  {
    if (m_previous_tract!=NULL)
    {
      if (n<m_previous_tract->GetNumberOfPoints())
      {
        return m_previous_tract->GetIndex(n);
      } else {
        return _points[n-m_previous_tract->GetNumberOfPoints()].index;
      }
    } else {
#ifdef _DEBUG      
        if(n>=GetNumberOfPoints())
          Error("Trying to read beyond the end of tract points");
#endif //_DEBUG        
        return _points[n].index;
    }
  }
  
  void Retain() //VF: why retain is virtual and Release is not !?!?!
  {
    _reference_count+=1;
    if (m_previous_tract!=NULL)
    {
      m_previous_tract->Retain();
    }
  }
  
  void Release()
  {

    if(_reference_count<2)
      delete(this);
    else
      _reference_count-=1;  
  }
 
  void Release(bool debug)
  {

    if(debug)
    {
      std::cout << "reference count " << _reference_count << endl;
    }
    Release();
  }
  
  bool CheckForDuplicateIndex(bool compare_last_only);
  void PrintPoints();
  int GetReferenceCount();
 

 protected:
  virtual ~FibreTract();
  
 private:
  //int _number_of_points;
  float _tract_length; /*IL*/
  std::vector<tract_point> _points; //VF: we live in 21sth century 
  float _tract_connectivity_index;
  float _mean_DOTT;
  float _min_DOTT;
  float _min_OTT;
  
  //int _points_allocated;
  //int _point_allocation_step;
  
  int _reference_count;
  FibreTract *m_previous_tract;
  //void _ReallocateTract();
};

#endif
