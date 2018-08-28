/************************************************************************

   File: ThreeOrthogonalPlanesPlotter.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __ThreeOrthogonalPlanesPlotter_h
#include "ThreeOrthogonalPlanesPlotter.h"
#endif

#ifndef __otherlibs_include_h
#include "otherlibs_include.h"
#endif

#ifndef __ImageData_h
#include "ImageData.h"
#endif

#ifndef __misc_h
#include "misc.h"
#endif


ThreeOrthogonalPlanesPlotter::ThreeOrthogonalPlanesPlotter(ImageData *image_data)
{
  int i;
  int x,y,z;
  int ijk[3];

  //create the vtkImageData: DO: check if need to flip (I'm pretty sure no)
  switch (image_data->GetImageType())
    {
    case Image3D_vanilla:
      _image_data_vtk=vtkImageData::New();
      _image_data_vtk->SetOrigin((double)image_data->GetXStart(),(double)image_data->GetYStart(),(double)image_data->GetZStart());
      _image_data_vtk->SetSpacing((double)image_data->GetXStep(),(double)image_data->GetYStep(),(double)image_data->GetZStep());
      _image_data_vtk->SetDimensions(image_data->GetXSize(),image_data->GetYSize(),image_data->GetZSize());
      _image_data_vtk->SetScalarTypeToFloat();//CLAUDE this fixes seg fault!!!!!!!
    
      // should keep a pointer of this for float_array->Delete();   // CLAUDE
      vtkFloatArray *float_array=vtkFloatArray::New();
      /*float_array->Allocate( image_data->GetXSize() * image_data->GetYSize() *
		             image_data->GetZSize() );   // CLAUDE*/
      

      for(x=0;x<image_data->GetXSize();x++)
	{
	  for(y=0;y<image_data->GetYSize();y++) 
	    { 
	      for(z=0;z<image_data->GetZSize();z++) 
		{
	
		  ijk[0]=x;
		  ijk[1]=y;
		  ijk[2]=z;

		  float_array->InsertValue(_image_data_vtk->ComputePointId(ijk),image_data->GetValue(x,y,z));

		  
		}
	    }
	}
      
      _image_data_vtk->GetPointData()->SetScalars(float_array);
      float_array->Delete(); //ilana clean up
      //this doesn't work: want SetScalars
      //_image_data_vtk->GetPointData()->AddArray(float_array);
      break;
      
    }
  //create the picker:
  vtkCellPicker *picker = vtkCellPicker::New();
  //picker->SetTolerance(0.005);
  

  //there is a default lookup table, so this isn't necessary:
  //vtkWindowLevelLookupTable *lookup_table=vtkWindowLevelLookupTable::New();

  //this makes no difference:
  //float maximum_value=image_data->GetMaximumValue();
  //float minimum_value=image_data->GetMinimumValue();

  //lookup_table->SetWindow((maximum_value-minimum_value)/(4*maximum_value));
  //lookup_table->SetLevel(minimum_value+(maximum_value-minimum_value)/(2*maximum_value));
  

  //do I need to place the widget? no
  //float bounds[6];
  
  //bounds[0]=image_data->GetXStart();
  //bounds[1]=image_data->GetXSize()*image_data->GetXStep();
  //bounds[2]=image_data->GetYStart();
  //bounds[3]=image_data->GetYSize()*image_data->GetYStep();
  //bounds[4]=image_data->GetZStart();
  //bounds[5]=image_data->GetZSize()*image_data->GetZStep();

  for (i=0;i<3;i++) 
    {
      _image_plane_widgets[i]=vtkImagePlaneWidget::New();
      _image_plane_widgets[i]->SetInput(_image_data_vtk);
      _image_plane_widgets[i]->DisplayTextOn();
      _image_plane_widgets[i]->RestrictPlaneToVolumeOn();
      //_image_plane_widgets[i]->UserControlledLookupTableOn();            

      switch (i)
	{
	case 0:
	  
	  _image_plane_widgets[i]->SetPlaneOrientationToXAxes();
	  //_image_plane_widgets[i]->SetKeyPressActivationValue("x");
	  _image_plane_widgets[i]->SetSliceIndex((int)image_data->GetXSize()/2);
	  break;
	case 1:
	  _image_plane_widgets[i]->SetPlaneOrientationToYAxes();
	  //_image_plane_widgets[i]->SetKeyPressActivationValue("y");
	  _image_plane_widgets[i]->SetSliceIndex((int)image_data->GetYSize()/2);
	
	  break;
	case 2: 
	  _image_plane_widgets[i]->SetPlaneOrientationToZAxes();
	  //_image_plane_widgets[i]->SetKeyPressActivationValue("z");
	  _image_plane_widgets[i]->SetSliceIndex((int)image_data->GetZSize()/2);
	  
	  break;
	}
      
      if (i==0)
	{
	  //_image_plane_widgets[i]->SetLookupTable(lookup_table);
	}
      
      else//if (i>0)
	{	 
	 _image_plane_widgets[i]->SetLookupTable(_image_plane_widgets[0]->GetLookupTable());
	}
      _image_plane_widgets[i]->SetPicker(picker);
      //_image_plane_widgets[i]->PlaceWidget(bounds);
      
      
    }
}

ThreeOrthogonalPlanesPlotter::~ThreeOrthogonalPlanesPlotter()
{
  int i;
  for (i=0;i<3;i++)
    {
      //_image_plane_widgets[i]->Delete(); //ilana this causes a seg fault
    }
  _image_data_vtk->Delete();
  
  
  //I don't know whether the float array defined in contructor gets deleted with the image_data_vtk
  
}

vtkImagePlaneWidget **ThreeOrthogonalPlanesPlotter::GetWidgets()
{
  return _image_plane_widgets;
  
}
