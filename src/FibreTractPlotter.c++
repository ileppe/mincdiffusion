/************************************************************************

   File: FibreTractPlotter.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#include <vtkPolyDataReader.h>
#include <vtkTubeFilter.h>
#include <vtkWindowLevelLookupTable.h>
#include <vtkScalarBarActor.h>

#include "FibreTractPlotter.h"
//#include "otherlibs_include.h"
#include "misc.h"
#include "FibreTractSet.h"



FibreTractPlotter::FibreTractPlotter(char *tracts_filename)
{
  vtkPolyDataReader *tracts_reader=vtkPolyDataReader::New();
  tracts_reader->SetFileName(tracts_filename);
  _tracts=tracts_reader->GetOutput();
  _plot_as_tubes=false;
  _red=2;
  _green=0;
  _blue=0;
  _map_point_scalars=false;
  //_map_point_rgb=true;
  m_colour_by_direction=false;
  m_tracts_filename=tracts_filename;
  m_map_file_scalars=false;
  m_file_low=0.0;
  m_file_high=1.0;
}

FibreTractPlotter::~FibreTractPlotter()
{
}

void FibreTractPlotter::SetTubesOn()
{
  _plot_as_tubes=true;
  _tube_radius=0.2; //default
  
}

void FibreTractPlotter::SetTubesOn(float tube_radius)
{
  _plot_as_tubes=true;
  _tube_radius=tube_radius;  
  
}

void FibreTractPlotter::SetColour(double red, double green, double blue)
{
  _red=red;
  _green=green;
  _blue=blue;
  
}

void FibreTractPlotter::MapPointScalarsOn()
{
  //DO: check that the tracts have point scalars and set LUT to have reasonable range.
  _map_point_scalars=true;
}

/*
void FibreTractPlotter::MapPointRGBOn()
{
  //check that the tracts have point scalars and set LUT to have reasonable range.
  _map_point_rgb=true;
  //_tracts->
}
*/


void FibreTractPlotter::ColourByScalarFile(char* tract_scalars_name, float low, float high)
{

  m_colour_by_direction=false;
  m_map_file_scalars=true;
  
  m_lut= vtkLookupTable::New();
  //m_lut->SetRange(0.0, 1.0);
  if (high<FLT_MAX)
    {
      m_file_low=low;
      m_file_high=high;
      m_lut->SetRange(m_file_low, m_file_high); 
      
    }
  //m_lut->SetRange(m_file_low, m_file_high); 


  m_lut->SetHueRange(1.0,0.0);//blue to red 
  m_lut->SetValueRange(0.6,1.0);
  //m_lut->SetSaturationRange(0.0,1.0);
  
  m_lut->Build();  

  

  
  ImageData *scalarsimage=new ImageData(tract_scalars_name);
  //not nec. - done in GetImageGeom() the first time 
  //scalarsimage->GetImageGeom()->point_inverse_xfm=NULL;
 

  POINT3D vpoint;
  
  vtkFloatArray *scalars=vtkFloatArray::New(); 

  vtkIdList* pointIDlist=vtkIdList::New();

  double point1[3];
  double point2[3];

  VECTOR3D direction;




  //here's the read-it-myself version: I assume a very strict file structure, not general polydata file (i.e., we had to write it).  

  //vtkPolyData *tracts;
  vtkPoints *points=vtkPoints::New();
  FILE *fp;
  fp=fopen(m_tracts_filename,"r");
  int ii,jj;
  char tmp[200];
  int numpoints;
  int numlines;
  
  float pointxf;
  float pointyf;
  float pointzf;
  
  float *pointx;
  float *pointy;
  float *pointz;
  
  bool ascii=false;

  //float test;
  /*
  float pointx;
  float pointy;
  float pointz;
  */

  float tmp2;

  int pointindex1,pointindex2;

 

  vtkIdList *line;

  IMAGE_GEOM image_geom_tmp;

  //preamble text:
  for (ii=0;ii<10;ii++)
    {
      fscanf(fp,"%s",tmp);
      if (ii==6)
	{
	  if (!strcmp(tmp,"ASCII"))
	    {
	      ascii=true;
	      m_colour_by_direction=true;
	    }
	  else
	    {
	      m_colour_by_direction=false;
	      return;
	    }
	}
    }
  //number of points:
  fscanf(fp,"%i",&numpoints);
  //datatype:
  fscanf(fp,"%s",tmp);


  //not at the right place here: need two more places for \n?
  //fscanf(fp,"%s",tmp);
  
  //cout << "numpoints " << numpoints << endl;

  //now put them in the polydata:

  //if not ascii, need offset to get to the points: 
  if (!ascii)
    {
      fread(&tmp2,1,1,fp);
    }

  vtkCellArray *lines=vtkCellArray::New();
  for (ii=0;ii<numpoints;ii++)
    {
      if (ascii)
	{
	  fscanf(fp,"%f",&pointxf);
	  fscanf(fp,"%f",&pointyf);
	  fscanf(fp,"%f",&pointzf);
	}
      else
	{
	  fread(&pointxf,sizeof(float),1,fp);
	  vtkByteSwap::Swap4BE(&pointxf);
	  fread(&pointyf,sizeof(float),1,fp);
	  vtkByteSwap::Swap4BE(&pointyf);
	  fread(&pointzf,sizeof(float),1,fp);
	  vtkByteSwap::Swap4BE(&pointzf);
	}

      //cout << "point " << pointxf << " " << pointyf << " " << pointzf << endl;
      points->InsertPoint(ii,pointxf,pointyf,pointzf);
      //points->InsertPoint(ii,pointx,pointy,pointz);
    }



  //lines:
  fscanf(fp,"%s",tmp);//"LINES"
  //cout << "LINES " << tmp << endl;
  fscanf(fp,"%i",&numlines);
  //cout << "numlines " << numlines << endl;
  fscanf(fp,"%s",tmp);//size 
  //cout << "size " << tmp << endl;
  
  //if I don't do this, numbers are offset below. need 4 bytes
  if (!ascii)
    {
      fread(&tmp2,4,1,fp);
    }

  //as we go through the indices, we set the derivative scalars

  //byte swapping seems not to work from here on


  for (ii=0;ii<numlines;ii++)
  
    {

      if (ascii)
	{
	  fscanf(fp,"%i",&numpoints);
	  //fscanf(fp,"%i",&pointindex1);
	}
      else
	{
	  fread(&numpoints,sizeof(int),1,fp);//numpoints is now number of points per line 
	  
	  fread(&pointindex1,sizeof(int),1,fp);
	}

   
      //cout << "numpoints " << numpoints << endl;
      //cout << "pointindex " << pointindex1 << endl;
      

      line=vtkIdList::New();
      //line->InsertNextId(pointindex1);
      //cout << "new line\n";
 
      image_geom_tmp=scalarsimage->GetImageGeom();

      for (jj=0;jj<numpoints;jj++)
	{
	  //if (jj>0)
	  //{
	  //  pointindex1=pointindex2;
	  // }
	 
	  
	  fscanf(fp,"%i",&pointindex1);
	  line->InsertNextId(pointindex1);
	  points->GetPoint(pointindex1,point1);
	  /*
	  if (ascii)
	    {
	      fscanf(fp,"%i",&pointindex2);
	    }
	  else
	    {
	      fread(&pointindex2,sizeof(int),1,fp);
	    }
	  */
	  //cout << "pointindex " << pointindex2 << endl;

	  //line->InsertNextId(pointindex2);





	  //points->GetPoint(pointindex2,point2);

	  //direction.x=point2[0]-point1[0];
	  //direction.y=point2[1]-point1[1];
	  //direction.z=point2[2]-point1[2];


	  //direction=Normalize(direction);
	  
	  
	   
	  //test:
	  //scalars->InsertNextValue(0.1);
	  vpoint=WorldToVoxel((float) point1[0], (float) point1[1], (float) point1[2],&image_geom_tmp);
	  //cout << scalarsimage->GetValue(vpoint.x,vpoint.y,vpoint.z) << endl;
	  
	  
	  if (scalarsimage->GetValue(vpoint.x,vpoint.y,vpoint.z) < low)
	    {
	      scalars->InsertNextValue(low);
	      //test=low;
	    }
	  else if (scalarsimage->GetValue(vpoint.x,vpoint.y,vpoint.z) > high)
	    {
	      scalars->InsertNextValue(high);
	      //test=high;
	    }
	  else 
	    {
	      scalars->InsertNextValue(scalarsimage->GetValue(vpoint.x,vpoint.y,vpoint.z));
	      //test=scalarsimage->GetValue(vpoint.x,vpoint.y,vpoint.z);
	    }
	  //cout << test << endl;
	    
	}
      
      //vpoint=WorldToVoxel((float) point2[0], (float) point2[1], (float) point2[2],scalarsimage->GetImageGeom());
      //scalars->InsertNextValue(scalarsimage->GetValue(vpoint.x,vpoint.y,vpoint.z));
      

      lines->InsertNextCell(line);
      line->Delete();
    }
  


  //give the scalars to the tracts we have already:

  //aborts here:
  //_tracts->Delete();

  _tracts=vtkPolyData::New();
  _tracts->SetPoints(points);
  _tracts->SetLines(lines);
  _tracts->GetPointData()->SetScalars(scalars);
  
  scalarsimage->Release();
   
 


}

vtkActor *FibreTractPlotter::ColourBar()
{
  vtkScalarBarActor *scalarBar =  vtkScalarBarActor::New();
  scalarBar->SetLookupTable(m_lut);
  //scalarBar->SetTitle("Scalars");
  scalarBar->SetNumberOfLabels(10);
  return (vtkActor*) scalarBar;

}


void FibreTractPlotter::ColourByDirection()
{
  


  m_lut= vtkLookupTable::New();
  m_lut->SetRange(0.0, 1.0);
  m_lut->Build();  



  //here we set scalars for each point on the tracts based on the orientation: these scalars will get mapped to colours 


  //DO: combine opacity and directional coding


  if (1)//full comment
    {
  
  vtkFloatArray *scalars=vtkFloatArray::New(); 

  vtkIdList* pointIDlist=vtkIdList::New();

  double point1[3];
  double point2[3];

  VECTOR3D direction;




  //here's the read-it-myself version: I assume a very strict file structure, not general polydata file (i.e., we had to write it).  

  //vtkPolyData *tracts;
  vtkPoints *points=vtkPoints::New();
  FILE *fp;
  fp=fopen(m_tracts_filename,"r");
  int ii,jj;
  char tmp[200];
  int numpoints;
  int numlines;
  
  float pointxf;
  float pointyf;
  float pointzf;
  
  float *pointx;
  float *pointy;
  float *pointz;
  
  bool ascii=false;


  /*
  float pointx;
  float pointy;
  float pointz;
  */

  float tmp2;

  int pointindex1,pointindex2;

 

  vtkIdList *line;

  //preamble text:
  for (ii=0;ii<10;ii++)
    {
      fscanf(fp,"%s",tmp);
      if (ii==6)
	{
	  if (!strcmp(tmp,"ASCII"))
	    {
	      ascii=true;
	      m_colour_by_direction=true;
	    }
	  else
	    {
	      m_colour_by_direction=false;
	      return;
	    }
	}
    }
  //number of points:
  fscanf(fp,"%i",&numpoints);
  //datatype:
  fscanf(fp,"%s",tmp);


  //not at the right place here: need two more places for \n?
  //fscanf(fp,"%s",tmp);
  
  //cout << "numpoints " << numpoints << endl;

  //now put them in the polydata:

  //if not ascii, need offset to get to the points: 
  if (!ascii)
    {
      fread(&tmp2,1,1,fp);
    }

  vtkCellArray *lines=vtkCellArray::New();
  for (ii=0;ii<numpoints;ii++)
    {
      if (ascii)
	{
	  fscanf(fp,"%f",&pointxf);
	  fscanf(fp,"%f",&pointyf);
	  fscanf(fp,"%f",&pointzf);
	}
      else
	{
	  fread(&pointxf,sizeof(float),1,fp);
	  vtkByteSwap::Swap4BE(&pointxf);
	  fread(&pointyf,sizeof(float),1,fp);
	  vtkByteSwap::Swap4BE(&pointyf);
	  fread(&pointzf,sizeof(float),1,fp);
	  vtkByteSwap::Swap4BE(&pointzf);
	}

      //cout << "point " << pointxf << " " << pointyf << " " << pointzf << endl;
      points->InsertPoint(ii,pointxf,pointyf,pointzf);
      //points->InsertPoint(ii,pointx,pointy,pointz);
    }



  //lines:
  fscanf(fp,"%s",tmp);//"LINES"
  //cout << "LINES " << tmp << endl;
  fscanf(fp,"%i",&numlines);
  //cout << "numlines " << numlines << endl;
  fscanf(fp,"%s",tmp);//size 
  //cout << "size " << tmp << endl;
  
  //if I don't do this, numbers are offset below. need 4 bytes
  if (!ascii)
    {
      fread(&tmp2,4,1,fp);
    }

  //as we go through the indices, we set the derivative scalars

  //byte swapping seems not to work from here on


  for (ii=0;ii<numlines;ii++)
  
    {

      if (ascii)
	{
	  fscanf(fp,"%i",&numpoints);
	  fscanf(fp,"%i",&pointindex1);
	}
      else
	{
	  fread(&numpoints,sizeof(int),1,fp);//numpoints is now number of points per line 
	  
	  fread(&pointindex1,sizeof(int),1,fp);
	}

   
      //cout << "numpoints " << numpoints << endl;
      //cout << "pointindex " << pointindex1 << endl;
      

      line=vtkIdList::New();
      line->InsertNextId(pointindex1);

 
      
      for (jj=0;jj<numpoints-1;jj++)
	{
	  if (jj>0)
	    {
	      pointindex1=pointindex2;
	    }
	 
	  
	  
	 
	  points->GetPoint(pointindex1,point1);

	  if (ascii)
	    {
	      fscanf(fp,"%i",&pointindex2);
	    }
	  else
	    {
	      fread(&pointindex2,sizeof(int),1,fp);
	    }
	  
	  //cout << "pointindex " << pointindex2 << endl;

	  line->InsertNextId(pointindex2);





	  points->GetPoint(pointindex2,point2);

	  direction.x=point2[0]-point1[0];
	  direction.y=point2[1]-point1[1];
	  direction.z=point2[2]-point1[2];


	  direction=Normalize(direction);
	  
	  //test: there is one 0
	  /*
	  if (ScalarForRGB(direction,m_lut)<FLT_EPSILON)
	    {
	      cout << ScalarForRGB(direction,m_lut) << endl;
	      cout << direction.x << " " << direction.y << " " << direction.z << endl; 
	      scalars->InsertNextValue(1.0);
	    }
	  else 
	    {
	  */
	  
	   
	  //test:
	  //scalars->InsertNextValue(0.5);
	  
	  scalars->InsertNextValue(ScalarForRGB(direction,m_lut));
	    
	}
      //last point gets same as second to last:

      scalars->InsertNextValue(ScalarForRGB(direction,m_lut));
      //test:
      //scalars->InsertNextValue(0.5);
      lines->InsertNextCell(line);
      line->Delete();
    }
  

  /*

SetLines(vtkCellArray *l)
vtkpolyline - set of 1D lines 
don't see how to set anything
vtkline (prob could add one by one to polydata)
don't see how to set anything
neither is a cellarray
want vtkCellArray directly
InsertCellPoint(vtkIdType id)
prob. needs to be line by line: add them line by line, because I don't see how to delineate the lines.  see if can do multiple times


//might as well set the scalars as we go through this

  */

  //this is my attempt to do this with vtk's reader, which I failed at:
  /*
  for (i=0;i<_tracts->GetNumberOfCells();i++)
    {
      //first point gets the same derivative as second:
      //cout << "Cell " << i <<endl;
      _tracts->GetCellPoints(i,pointIDlist);
      for (j=0;j<pointIDlist->GetNumberOfIds()-1;j++)
	{
	  _tracts->GetPoint(pointIDlist->GetId(j),point1);
	  _tracts->GetPoint(pointIDlist->GetId(j+1),point2);

	  direction.x=point2[0]-point1[0];
	  direction.y=point2[1]-point1[1];
	  direction.z=point2[2]-point1[2];

	  scalars->SetValue(j,ScalarForRGB(direction,m_lut));
	  cout << ScalarForRGB(direction,m_lut) << endl;
	}
      //last point gets same as second to last:
      scalars->SetValue(j,ScalarForRGB(direction,m_lut));



    }
  */


  //give the scalars to the tracts we have already:

  //aborts here:
  //_tracts->Delete();

  _tracts=vtkPolyData::New();
  _tracts->SetPoints(points);
  _tracts->SetLines(lines);
  _tracts->GetPointData()->SetScalars(scalars);
  
  
    }//end full comment
 


}


vtkActor *FibreTractPlotter::CreateActor()
{
  float range[2];
  
  vtkPolyDataMapper *tracts_mapper = vtkPolyDataMapper::New();
  if (!_plot_as_tubes)
    {
      
      tracts_mapper->SetInput(_tracts);
    }
  else
    {
      vtkTubeFilter *tube_filter=vtkTubeFilter::New();
      tube_filter->SetInput(_tracts);
      
      tube_filter->SetNumberOfSides(3);
      tube_filter->SetRadius(_tube_radius);
      tube_filter->CappingOn();
      //tube_filter->SetVaryRadiusToVaryRadiusByScalar();      
      tracts_mapper->SetInput(tube_filter->GetOutput());

    }


  //if(map_tract_scalars)
  //{
  //tracts_mapper->SetScalarModeToUseCellData();
  //etc. 
  //tracts_mapper->SetScalarRange(_tracts->GetCellData->GetScalars()->GetRange());
  //}

  

  if (_map_point_scalars)
    {
      
      

      //color code
      if (1)
	{
	  //vtkLookupTable *lookup_table = vtkLookupTable::New();
      //this makes no difference: it is using the LUT in the file!?
      //lookup_table->SetRange(0.0, 200.0);
      //lookup_table->SetNumberOfColors(20);
      //lookup_table->Build(); 
      tracts_mapper->ScalarVisibilityOn();
      
      tracts_mapper->SetColorModeToMapScalars();
      tracts_mapper->SetScalarModeToUsePointData();

      //tracts_mapper->SetLookupTable(lookup_table);
      //lookup_table->Delete();
	}


      

      //opacity encode
      if (0)
	{

	  tracts_mapper->ScalarVisibilityOn();
      
	  tracts_mapper->SetColorModeToMapScalars();
	  tracts_mapper->SetScalarModeToUsePointData();

	  vtkWindowLevelLookupTable *window_level_lookup_table = vtkWindowLevelLookupTable::New();


      //0: transparent
      //1: opaque
      window_level_lookup_table->SetMinimumTableValue(_red,_green,_blue,0);
      window_level_lookup_table->SetMaximumTableValue(_red,_green,_blue,0.1);
      



      window_level_lookup_table->Build();          
      tracts_mapper->SetLookupTable(window_level_lookup_table);
      window_level_lookup_table->Delete();
	}



   }
  /*
  else if(_map_point_rgb) //not sure why this is here
   {
      vtkLookupTable *lookup_table = vtkLookupTable::New();
      tracts_mapper->ScalarVisibilityOn();
      tracts_mapper->SetColorModeToMapScalars();
      tracts_mapper->SetScalarModeToUsePointData();
      lookup_table->SetHueRange(1.0,0.0);
      //window_level_lookup_table->SetValueRange(0.9,0.9);
      //window_level_lookup_table->GetWindow();
      //window_level_lookup_table->SetLevel(0.5);
   
      lookup_table->Build();          
      tracts_mapper->SetLookupTable(lookup_table);
   }
  */
  else if (m_colour_by_direction || m_map_file_scalars)
    {
      
      
      tracts_mapper->ScalarVisibilityOn();
      tracts_mapper->SetLookupTable(m_lut);
      tracts_mapper->SetColorModeToMapScalars();
      tracts_mapper->SetScalarModeToUsePointData();

      tracts_mapper->SetScalarRange(m_file_low,m_file_high);
      //m_lut->Delete();
      
    }
  else 
    {
      tracts_mapper->ScalarVisibilityOff();
    }  

  tracts_mapper->ImmediateModeRenderingOn();  
  vtkActor *tracts_actor=vtkActor::New();
  tracts_actor->SetMapper(tracts_mapper);
  tracts_mapper->Delete();

  if (_map_point_scalars )
    {
      //tracts_actor->GetProperty()->SetOpacity(0.99);
 
    }
  else if (m_colour_by_direction)
    {
      //do nothing
 
    }
  else
    {
      tracts_actor->GetProperty()->SetColor(_red,_green,_blue);
      
    }
  
  tracts_actor->GetProperty()->SetLineWidth(3);
  return tracts_actor;
  


}


