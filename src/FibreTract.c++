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
#include <float.h>
#include "FibreTract.h"
#include "otherlibs_include.h"

FibreTract::FibreTract()
{
  m_previous_tract=NULL;
  //_number_of_points=0;
  _tract_length=0; /*IL*/
  _tract_connectivity_index=0;
  //_point_allocation_step=10;

  //_points=(TRACT_POINT *)malloc(_point_allocation_step*sizeof(tract_point));


  //_points_allocated=_point_allocation_step;
  _reference_count=1;

}

FibreTract::FibreTract(FibreTract *tract1, FibreTract *tract2)
{
  int i;

  //_number_of_points=0;
  _tract_connectivity_index=0;
  //_point_allocation_step=10;

  //_points=(TRACT_POINT *)malloc(_point_allocation_step*sizeof(TRACT_POINT));


  //_points_allocated=_point_allocation_step;
  _reference_count=1;

  for (i=tract1->GetNumberOfPoints()-1;i>-1;i--) {
    AddPoint(tract1->GetPoint(i));
    SetPointInfo(tract1->GetPointInfo(i),GetNumberOfPoints()-1);
  }

  for (i=1;i<tract2->GetNumberOfPoints();i++) { //VF: index is starting at 1 ?
    AddPoint(tract2->GetPoint(i));
    SetPointInfo(tract2->GetPointInfo(i),GetNumberOfPoints()-1);
  }

  m_previous_tract=NULL;
}


FibreTract::FibreTract(FibreTract *tract)
{
  m_previous_tract=tract;
  m_previous_tract->Retain();

  //_number_of_points=0;
  //_tract_connectivity_index=0;
  //_point_allocation_step=10;

  //_points=(TRACT_POINT *)malloc(_point_allocation_step*sizeof(TRACT_POINT));


  //_points_allocated=_point_allocation_step;
  _reference_count=1;
}


FibreTract::~FibreTract()
{
  //free(_points);
  if (m_previous_tract!=NULL) {
    m_previous_tract->Release();
  }

  //  cout << "delete tract\n";
}







void FibreTract::AddPoint(float x, float y, float z)
{
  AddPoint(Point3D(x,y,z));

}

void FibreTract::AddPoint(const POINT3D &point)
{
  /*  if (_number_of_points==_points_allocated)
    {
        _ReallocateTract();

    }
    _number_of_points+=1;

    _points[_number_of_points-1].point=point;*/
  _points.push_back(tract_point(point));
  /*Calculate tract length*/ //IL

  if (_points.size()>1) {
    _tract_length = _tract_length + euclidian_distance(_points[_points.size()-2].point,point);
  }
}


void FibreTract::SetTractConnectivityIndex(float tract_connectivity_index)
{
  _tract_connectivity_index=tract_connectivity_index;

}

float FibreTract::GetTractConnectivityIndex()
{
  return _tract_connectivity_index;
}

void FibreTract::SetCumulativeCI(float cumulative_CI, int n)
{
  _points[n].cumulative_CI=cumulative_CI;

}

float FibreTract::GetCumulativeCI(int n)
{
  return _points[n].cumulative_CI;

}

void FibreTract::SetDOTT(float DOTT, int n)
{
  _points[n].DOTT=DOTT;

  if (_min_DOTT>DOTT || _points.size()==1) {
    _min_DOTT=DOTT;

  }

  _points[n].cumulative_min_DOTT=_min_DOTT;


}

float FibreTract::GetDOTT(int n)
{
  return _points[n].DOTT;

}

void FibreTract::SetMeanDOTT(float mean_DOTT)
{
  _mean_DOTT=mean_DOTT;

}

float FibreTract::GetMeanDOTT()
{
  return _mean_DOTT;

}



void FibreTract::SetMinDOTT(float min_DOTT)
{
  _min_DOTT=min_DOTT;

}

float FibreTract::GetMinDOTT()
{
  return _min_DOTT;

}



float FibreTract::GetCumulativeMinDOTT(int n)
{
  return _points[n].cumulative_min_DOTT;

}




void FibreTract::SetOTT(float OTT, int n)
{
  _points[n].OTT=OTT;

  if (_min_OTT>OTT || _points.size()==1) {
    _min_OTT=OTT;

  }

  _points[n].cumulative_min_OTT=_min_OTT;



}

float FibreTract::GetOTT(int n)
{
  return _points[n].OTT;

}


//note: I am not sure whether this function alone has changed all maps requiring the blurred value in FibreTracking.c++
//note: really need to recalculate the values for the actual point offsets, which aren't voxel dim exactly.
//all of this checking for nan is just debugging: this routine results in non-monotonic CMOTT and I haven't figured out why
float FibreTract::GetBlurredOTT(int n)
{
  //tracts shorter than 5 will not be blurred.  near the seed won't be blurred either

  float toohigh=1000;//this needs to be adjusted by the user for experiments were legitimate high values occur: for debugging here

  //if (_number_of_points<5)
  if (_points.size()<5)	  //IRL 
    {
      return this->GetOTT(n);
    }

  if (n==0)
    {
      if (_points[n].OTT>toohigh || _points[n+1].OTT>toohigh || _points[n+2].OTT>toohigh )
	{
	  return _points[n].OTT;
	}
	   
      //cout << _points[n].OTT << " " << _points[n+1].OTT << " " << _points[n+2].OTT << endl;
      else
	{
	  return (0.0545+0.2442+0.4026)*_points[n].OTT+0.2442*_points[n+1].OTT+0.0545*_points[n+2].OTT;
	}
      
    }
  else if (n==1)
    {
      if (_points[n-1].OTT>toohigh || _points[n].OTT>toohigh || _points[n+1].OTT>toohigh || _points[n+2].OTT>toohigh )
	{
	  return _points[n].OTT;
	}
	   
      //cout << _points[n-1].OTT << " " << _points[n].OTT << " " << _points[n+1].OTT << " " << _points[n+2].OTT << endl;
      else
	{
	  return (0.0545+0.2442)*_points[n-1].OTT+0.4026*_points[n].OTT+0.2442*_points[n+1].OTT+0.0545*_points[n+2].OTT;
	}
    }
  //else if (n==_number_of_points-2)
    else if (n==_points.size()-2) //IRL
    {
      if (_points[n-1].OTT>toohigh || _points[n].OTT>toohigh || _points[n+1].OTT>toohigh || _points[n-2].OTT>toohigh )
	{
	  return _points[n].OTT;
	}
      //cout << _points[n-2].OTT  << " " << _points[n-1].OTT << " " << _points[n].OTT << " " << _points[n+1].OTT << endl;
      else
	{
	  return 0.0545*_points[n-2].OTT+0.2442*_points[n-1].OTT+0.4026*_points[n].OTT+(0.0545+0.2442)*_points[n+1].OTT;
	}
    }
  //else if (n==_number_of_points-1)
  else if(n==_points.size()-1) //IRL	  
    {

      if (_points[n-1].OTT>toohigh || _points[n].OTT>toohigh || _points[n-2].OTT>toohigh )
	{
	  return _points[n].OTT;
	} 
      //cout << _points[n-2].OTT  << " " << _points[n-1].OTT << " " << _points[n].OTT  << endl;
      else
	{
      return 0.0545*_points[n-2].OTT+0.2442*_points[n-1].OTT+(0.0545+0.2442+0.4026)*_points[n].OTT;
	}
    }
  else
    {

      if ( _points[n-2].OTT>toohigh || _points[n-1].OTT>toohigh || _points[n].OTT>toohigh || _points[n+1].OTT>toohigh || _points[n+2].OTT>toohigh )
	{
	  return _points[n].OTT;
	}
      //cout << _points[n-2].OTT  << " " << _points[n-1].OTT << " " << _points[n].OTT << " " << _points[n+1].OTT << " " << _points[n+2].OTT << endl;
      else
	{
      return  0.0545*_points[n-2].OTT+0.2442*_points[n-1].OTT+0.4026*_points[n].OTT+0.2442*_points[n+1].OTT+0.0545*_points[n+2].OTT;
	}
    }
}

void FibreTract::SetMinOTT(float min_OTT)
{
  _min_OTT=min_OTT;

}

float FibreTract::GetMinOTT()
{
  return _min_OTT;

}

void FibreTract::SetCumulativeMinOTT(float CMOTT, int n)
{
  _points[n].cumulative_min_OTT=CMOTT;

}


float FibreTract::GetCumulativeMinOTT(int n)
{
  return _points[n].cumulative_min_OTT;

}






bool FibreTract::GoesThroughVoxel(ImageData *roi, int x, int y, int z)
{
  int i;
  POINT3D point;
  bool goes_through_voxel=false;

  for (i=0;i<_points.size() && !goes_through_voxel;i++) {
    point=GetPoint(i);


    if ((float_abs((roi->GetXStart()+roi->GetXStep()*x)-point.x)<=float_abs(roi->GetXStep()/2.0)) && (float_abs((roi->GetYStart()+roi->GetYStep()*y)-point.y)<=float_abs(roi->GetYStep()/2.0)) && (float_abs((roi->GetZStart()+roi->GetZStep()*z)-point.z)<=float_abs(roi->GetZStep()/2.0))) {
      goes_through_voxel=true;
    }

  }

  if (goes_through_voxel==false && m_previous_tract!=NULL) {
    goes_through_voxel=m_previous_tract->GoesThroughVoxel(roi, x, y, z);
  }

  return goes_through_voxel;
}


bool FibreTract::GoesThroughVoxel(IMAGE_GEOM image_geom, int x, int y, int z)
{
  int i;
  POINT3D point;
  bool goes_through_voxel=false;

  for (i=0;i<_points.size() && !goes_through_voxel;i++) {
    point=GetPoint(i);

    //VF: direction cosines not used ?!
    if ((fabs((image_geom.x_start+image_geom.x_step*x)-point.x)<=fabs(image_geom.x_step/2.0)) && 
        (fabs((image_geom.y_start+image_geom.y_step*y)-point.y)<=fabs(image_geom.y_step/2.0)) && 
        (fabs((image_geom.z_start+image_geom.z_step*z)-point.z)<=fabs(image_geom.z_step/2.0))) 
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



bool FibreTract::GoesThroughROI(ImageData *roi)
{
  int x,y,z;
  bool goes_through_roi=false;

  //going through all ROI voxels is slow, but otherwise need a way of assigning multiple voxel indices to edge points.  Hence, for all ROI voxels that are greater than zero, check if there is a tract point in it or on the edge

  for (x=0;x<roi->GetXSize() && !goes_through_roi;x++) {
    for (y=0;y<roi->GetYSize() && !goes_through_roi;y++) {
      for (z=0;z<roi->GetZSize() && !goes_through_roi;z++) {
        if (roi->GetValue(x,y,z)>FLT_EPSILON) {
          if (GoesThroughVoxel(roi,x,y,z)) {
            goes_through_roi=true;

          }
        }
      }
    }
  }

  return goes_through_roi;
}

//there might be a problem with this when tract length gets long
/*void FibreTract::_ReallocateTract()
{

  _points=(TRACT_POINT *) realloc(_points,(_points_allocated+_point_allocation_step)*sizeof(TRACT_POINT));

  _points_allocated=_points_allocated+_point_allocation_step;
  _point_allocation_step*=2;
}
*/



bool FibreTract::CheckForDuplicateIndex(bool compare_last_only)
{
  if (!compare_last_only) {
    Error("Must compare last only");

  }

  for (int i=0;i<(GetNumberOfPoints()-1);i++) 
    if ( GetIndex(GetNumberOfPoints()-1)==GetIndex(i) ) 
      return true;
    
  return false;
}

void FibreTract::PrintPoints()
{
  int i;

  cout << "Points " << GetNumberOfPoints() << endl;

  for (i=0;i<GetNumberOfPoints();i++) {
    cout << GetPoint(i).x << " " << GetPoint(i).y << " " << GetPoint(i).z << endl;
  }

  cout << endl;
}

int FibreTract::GetReferenceCount()
{
  return _reference_count;
}


