/************************************************************************

   File: ODFPlotter.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#include <float.h>
#include "ODFPlotter.h"
#include "otherlibs_include.h"
#include "ImageData.h"
#include "Directions.h"
#include "misc.h"
#include "S2Interpolator.h"



ODFPlotter::ODFPlotter(ImageData *ODF_image_data)
{
  float theta,phi;
  int counter;
  
  m_ODF_image_data=ODF_image_data;
  
  //set defaults:
  //m_colour_scheme=SingleColour;
  m_colour_scheme=RGB;
  m_single_colour[0]=0.5;
  m_single_colour[1]=0.5;
  m_single_colour[2]=0.5;
  

  m_angle_step=0.1745;
  //m_angle_step=0.04;
  //m_angle_step=0.08;//causes seg faults
  m_display_mode=Physical;
  m_scale_factor=1;

  m_smooth_shapes=false;
  m_smoother_iterations=100;
  m_specular_fraction=0.7;
  m_diffuse_fraction=0.4;
  m_ambient_fraction=0.3;
  m_specular_power=7.0;
  m_specular_color_red=0.6;
  m_specular_color_green=0.6;
  m_specular_color_blue=0.6;
  m_opacity=1.0;


  Directions *directions=new Directions((int)(ceil(2*M_PI/m_angle_step)*ceil(M_PI/m_angle_step)));


  m_plot_fan=false;
 //for CI, read in fullpdf to plot it instead of histo_most_likely.  must match.
  if (m_plot_fan)
    {
      cout << "plot fan\n";
      m_full_pdf=new ImageData("full_pdf.mnc");
      m_polarity=new ImageData("polarity.mnc");
    }
 

  counter=0;
  for (phi=0;phi<M_PI;phi+=m_angle_step)
    {
      for (theta=0;theta<2*M_PI;theta+=m_angle_step)
	{


	
	  directions->SetXComponent(counter, cos(theta)*sin(phi));
	  directions->SetYComponent(counter, sin(theta)*sin(phi));
	  directions->SetZComponent(counter, cos(phi));



	  counter++;
	}
    }
  

  m_interpolator=new S2Interpolator(m_ODF_image_data, directions, 1.5*sqrt(2*M_PI/m_ODF_image_data->GetNSize()), 1.5*sqrt(2*M_PI/m_ODF_image_data->GetNSize()), Gaussian);


  
  if (m_plot_fan)
    {
      m_interpolator_2=new S2Interpolator(m_full_pdf, directions, 1.5*sqrt(2*M_PI/m_ODF_image_data->GetNSize()), 1.5*sqrt(2*M_PI/m_ODF_image_data->GetNSize()), Gaussian);
    }
  
  m_rvalues=gsl_vector_float_alloc((int)(ceil(2*M_PI/m_angle_step)*ceil(M_PI/m_angle_step)));

  m_rvalues_are_allocated=true;




}


ODFPlotter::ODFPlotter(ImageData *ODF_image_data, ImageData *polarity, ImageData *full_pdf)
{
  m_plot_fan=false;
  m_polarity=polarity;
  m_full_pdf=full_pdf;

  float theta,phi;
  int counter;
  
  m_ODF_image_data=ODF_image_data;
  
  //set defaults:
  //m_colour_scheme=SingleColour;
  m_colour_scheme=RGB;
  m_single_colour[0]=0.5;
  m_single_colour[1]=0.5;
  m_single_colour[2]=0.5;
  

  m_angle_step=0.1745;
  //m_angle_step=0.04;
  //m_angle_step=0.08;//causes seg faults
  m_display_mode=Physical;
  m_scale_factor=1;

  m_smooth_shapes=false;
  m_smoother_iterations=100;
  m_specular_fraction=0.7;
  m_diffuse_fraction=0.4;
  m_ambient_fraction=0.3;
  m_specular_power=7.0;
  m_specular_color_red=0.6;
  m_specular_color_green=0.6;
  m_specular_color_blue=0.6;
  m_opacity=1.0;


  Directions *directions=new Directions((int)(ceil(2*M_PI/m_angle_step)*ceil(M_PI/m_angle_step)));



 

  counter=0;
  for (phi=0;phi<M_PI;phi+=m_angle_step)
    {
      for (theta=0;theta<2*M_PI;theta+=m_angle_step)
	{


	
	  directions->SetXComponent(counter, cos(theta)*sin(phi));
	  directions->SetYComponent(counter, sin(theta)*sin(phi));
	  directions->SetZComponent(counter, cos(phi));



	  counter++;
	}
    }
  

  m_interpolator=new S2Interpolator(m_ODF_image_data, directions, 1.5*sqrt(2*M_PI/m_ODF_image_data->GetNSize()), 1.5*sqrt(2*M_PI/m_ODF_image_data->GetNSize()), Gaussian);


  
  if (m_plot_fan)
    {
      m_interpolator_2=new S2Interpolator(m_full_pdf, directions, 1.5*sqrt(2*M_PI/m_ODF_image_data->GetNSize()), 1.5*sqrt(2*M_PI/m_ODF_image_data->GetNSize()), Gaussian);
    }
  
  m_rvalues=gsl_vector_float_alloc((int)(ceil(2*M_PI/m_angle_step)*ceil(M_PI/m_angle_step)));

  m_rvalues_are_allocated=true;




}


ODFPlotter::ODFPlotter(ImageData *ODF_image_data, float angle_step)
{


 
  float theta,phi;
  int counter;
  m_ODF_image_data=ODF_image_data;
  
  //set defaults:
  m_colour_scheme=SingleColour;
  m_single_colour[0]=1;
  m_single_colour[1]=0;
  m_single_colour[2]=0;

  m_opacity=1.0;  

  m_angle_step=angle_step;
  m_display_mode=Physical;
  m_scale_factor=1;

  Directions *directions=new Directions((int)(ceil(2*M_PI/m_angle_step)*ceil(M_PI/m_angle_step)));

  counter=0;
  for (phi=0;phi<M_PI;phi+=m_angle_step)
    {
      for (theta=0;theta<2*M_PI;theta+=m_angle_step)
	{


	  directions->SetXComponent(counter, cos(theta)*sin(phi));
	  directions->SetYComponent(counter, sin(theta)*sin(phi));
	  directions->SetZComponent(counter, cos(phi));
	  counter++;
	}
    }
  

  m_interpolator=new S2Interpolator(m_ODF_image_data, directions, 1.5*sqrt(2*M_PI/m_ODF_image_data->GetNSize()), 1.5*sqrt(2*M_PI/m_ODF_image_data->GetNSize()), Gaussian);



  //m_interpolator_2=new S2Interpolator(m_ODF_image_data, directions, 1.5*sqrt(2*M_PI/m_ODF_image_data->GetNSize()), 1.5*sqrt(2*M_PI/m_ODF_image_data->GetNSize()), Gaussian);

  m_rvalues=gsl_vector_float_alloc((int)(ceil(2*M_PI/m_angle_step)*ceil(M_PI/m_angle_step)));

  m_rvalues_are_allocated=true;


}


//DO: if no other pointers to the interpolator, destroy it.  need to do reference counting to do this
ODFPlotter::~ODFPlotter()
{
  if (m_rvalues_are_allocated)
    {
      gsl_vector_float_free(m_rvalues);
    }
  
}

//void ODFPlotter::SetPlotFan(int plot_fan)

void ODFPlotter::SetVoxel(int voxel_x, int voxel_y, int voxel_z)
{
  m_voxel[0]=voxel_x;
  m_voxel[1]=voxel_y;
  m_voxel[2]=voxel_z;
  
}




void ODFPlotter::CalculatePlottingPoints()
{
  float theta,phi;
  int counter;
  int i;
  float value;
  
  if (!m_rvalues_are_allocated) //should never happen 
    {
       m_rvalues=gsl_vector_float_alloc((int)(ceil(2*M_PI/m_angle_step)*ceil(M_PI/m_angle_step)));
       

      m_rvalues_are_allocated=true; 
    }
    
  

  counter=0;
  

  for (phi=0;phi<M_PI;phi+=m_angle_step)
    {
      for (theta=0;theta<2*M_PI;theta+=m_angle_step)
	{
	  
	  
	  //full_pdf not input file if polarity vector aligns:
	  if (m_plot_fan)
	    {
	      if (m_interpolator_2->GetDirections()->GetXComponent(counter)*m_polarity->GetValue(m_voxel[0],m_voxel[1],m_voxel[2],0)+m_interpolator_2->GetDirections()->GetYComponent(counter)*m_polarity->GetValue(m_voxel[0],m_voxel[1],m_voxel[2],1)+m_interpolator_2->GetDirections()->GetZComponent(counter)*m_polarity->GetValue(m_voxel[0],m_voxel[1],m_voxel[2],2)>0)
		{
		  //test:

		  //gsl_vector_float_set(m_rvalues,counter,m_interpolator_2->GetInterpolatedValue(counter,m_voxel[0],m_voxel[1],m_voxel[2]));
		  gsl_vector_float_set(m_rvalues,counter,m_interpolator_2->GetInterpolatedValue(counter,m_voxel[0],m_voxel[1],m_voxel[2])/3.6);
		}
	      else
		{
		  gsl_vector_float_set(m_rvalues,counter,m_interpolator->GetInterpolatedValue(counter,m_voxel[0],m_voxel[1],m_voxel[2]));
		}
	    }
	  else
	    {
	      gsl_vector_float_set(m_rvalues,counter,m_interpolator->GetInterpolatedValue(counter,m_voxel[0],m_voxel[1],m_voxel[2]));
	    }
	  counter+=1;
			     
	}
      
      
    
  
    }





  
  switch (m_display_mode) 
    {
    case Stretch:
      for (i=0;i<(int)(ceil(2*M_PI/m_angle_step)*ceil(M_PI/m_angle_step));i++)
	{
	  value=gsl_vector_float_get(m_rvalues,i);
	  gsl_vector_float_set(m_rvalues,i,pow(value,4));
	  
	}
      
      break;
    case SubtractMinimum:
      gsl_vector_float_add_constant(m_rvalues,-1*gsl_vector_float_min(m_rvalues));
      break;
    case Physical:
      //do nothing
      break;
      
    }
  
  
    }
  
 
float ODFPlotter::GetBestScaleFactorForVoxelFit()
{
  float scale_factor;
  float maximum_rvalue;
  float minimum_voxel_dimension;
  

  maximum_rvalue=(float)gsl_vector_float_max(m_rvalues);
  minimum_voxel_dimension=float_min3(m_ODF_image_data->GetXStep(),m_ODF_image_data->GetYStep(),m_ODF_image_data->GetZStep());
  

  scale_factor=minimum_voxel_dimension/(2*maximum_rvalue);
  
  return scale_factor;
  


}

void ODFPlotter::SetScaleFactor(float scale_factor)
{
  m_scale_factor=scale_factor;
  
}



  
 
vtkActor *ODFPlotter::CreateActor()
{
  float theta,phi;
  int counter;
  VECTOR3D direction;
  
  POINT3D point; 
    
  vtkPoints *points=vtkPoints::New();
  vtkFloatArray *scalars=vtkFloatArray::New(); 

  vtkPolyData *polydata=vtkPolyData::New();

  vtkCellArray *polygons=vtkCellArray::New();
  
  int *polygoncorners=(int *) malloc(4*sizeof(int));

  vtkPolyDataMapper *ODFMapper=vtkPolyDataMapper::New();

  vtkActor *ODFActor = vtkActor::New();

  vtkStripper *stripper = vtkStripper::New();
  vtkSmoothPolyDataFilter *smoother = vtkSmoothPolyDataFilter::New();
  

  counter=0;

  point=VoxelToWorld(m_voxel[0],m_voxel[1],m_voxel[2],m_ODF_image_data->GetImageGeom());



  //create the lut for the Radial or RGB color coding:
  //note that for radial encoding of things like the DWI profile, we really need to clamp this range to a tighter range, else all the same color still

  vtkLookupTable *lut= vtkLookupTable::New();
  //vtkScalarsToColors *lut= vtkScalarsToColors::New();  
  lut->SetRange(0.0, 1.0);
  //lut->SetNumberOfTableValues(1000);//default 256 is fine
  lut->Build();      



  //DO: could rewrite this as loop through directions; keep directions from constructor as _directions  
  for (phi=0;phi<M_PI;phi+=m_angle_step)
    {
      for (theta=0;theta<2*M_PI;theta+=m_angle_step)
	{

	 
	  

	  direction.x=cos(theta)*sin(phi);
	  direction.y=sin(theta)*sin(phi);
	  direction.z=cos(phi);


	 
	  //everything is done in world space
	  //direction=VoxelToWorld(direction,m_ODF_image_data->GetImageGeom());

	  //need to make the vector independant of signed direction:
	  
	  
	  points->InsertPoint(counter,point.x+m_scale_factor*gsl_vector_float_get(m_rvalues,counter)*direction.x, point.y+m_scale_factor*gsl_vector_float_get(m_rvalues,counter)*direction.y, point.z+m_scale_factor*gsl_vector_float_get(m_rvalues,counter)*direction.z);  
	  	  
	  if (m_colour_scheme==Radial)
	    {

	      scalars->InsertNextValue(gsl_vector_float_get(m_rvalues,counter));

	    }
	  
	  else if (m_colour_scheme==RGB)
	    {
	      scalars->InsertNextValue(ScalarForRGB(direction,lut));
	      
	    }

	  if (phi>FLT_EPSILON && theta>FLT_EPSILON) 
	    {
	      polygons->InsertNextCell(4);
	      /*
	      polygoncorners[0]=(int) counter;
	      polygoncorners[1]=(int) counter-1;
	      polygoncorners[2]=(int) counter-(int)ceil(2*M_PI/m_angle_step)-1;
	      polygoncorners[3]=(int) polygoncorners[2]+1;*/
	      
	      polygons->InsertCellPoint((int)counter);
	      polygons->InsertCellPoint((int) counter-1);
	      int tmp=(int) counter-(int)ceil(2*M_PI/m_angle_step)-1;
	      polygons->InsertCellPoint(tmp);
	      polygons->InsertCellPoint(tmp+1);
	      	      
	      //polygons->InsertNextCell(4,polygoncorners); 
	    }
	  
	  
	  counter++;
	  
	}
    }
  
  polydata->SetPoints(points);
  polydata->SetPolys(polygons);

  vtkPolyDataNormals *normals=vtkPolyDataNormals::New();

  if (m_smooth_shapes)
    {
      smoother->SetInput(polydata);
      smoother->SetNumberOfIterations(m_smoother_iterations);
      normals->SetInput(smoother->GetOutput());
      stripper->SetInput(normals->GetOutput());


      ODFMapper->SetInput(stripper->GetOutput());
    }
  else 
    {
      ODFMapper->SetInput(polydata);
    }
  
  
  ODFActor->SetMapper(ODFMapper);

  if (m_colour_scheme==SingleColour)
    {
      ODFActor->GetProperty()->SetColor(m_single_colour[0],m_single_colour[1],m_single_colour[2]);  
    }
  else
    {
    
      polydata->GetPointData()->SetScalars(scalars);
     
      ODFMapper->ScalarVisibilityOn();
      ODFMapper->SetLookupTable(lut);
      ODFMapper->SetColorModeToMapScalars();
      ODFMapper->SetScalarModeToUsePointData();
   
    }
  

  ODFActor->GetProperty()->SetOpacity(m_opacity);


  if (m_smooth_shapes)
    {
      
      ODFActor->GetProperty()->SetSpecular(m_specular_fraction);
      ODFActor->GetProperty()->SetDiffuse(m_diffuse_fraction);   
      ODFActor->GetProperty()->SetAmbient(m_ambient_fraction);          
      ODFActor->GetProperty()->SetSpecularPower(m_specular_power);
      ODFActor->GetProperty()->SetSpecularColor(m_specular_color_red, m_specular_color_green, m_specular_color_blue);  
      ODFActor->GetProperty()->SetInterpolationToGouraud();
    }
  

  free(polygoncorners);

 
  scalars->Delete();
  polygons->Delete();
  polydata->Delete();
  points->Delete();
  ODFMapper->Delete();
  lut->Delete();
  normals->Delete();


  return(ODFActor);
  
	  
}
      


void ODFPlotter::SetColour(float red, float green, float blue)
{
  m_single_colour[0]=red;
  m_single_colour[1]=green;
  m_single_colour[2]=blue;
  m_colour_scheme=SingleColour;
}

void ODFPlotter::SetColourScheme(COLOUR_SCHEME colour_scheme)
{
  m_colour_scheme=colour_scheme;
  
}




void ODFPlotter::SetDisplayMode(DISPLAY_MODE display_mode)
{
  m_display_mode=display_mode;
  
}

S2Interpolator *ODFPlotter::GetInterpolator()
{
  return m_interpolator;
  
}

S2Interpolator *ODFPlotter::GetInterpolator2()
{
  return m_interpolator_2;
  
}


void ODFPlotter::SetInterpolator(S2Interpolator *interpolator)
{
  m_interpolator=interpolator;
  
}

void ODFPlotter::SetInterpolator2(S2Interpolator *interpolator_2)
{
  m_interpolator_2=interpolator_2;
  
}


void ODFPlotter::SetSmoothOn()
{
  m_smooth_shapes=true;
  
}

void ODFPlotter::SetOpacity(double opacity) //ilana vtk5 SetOpacity takes double
{
  m_opacity=opacity;

}
