/************************************************************************

   File: ImageData.c++

   Author: Jennifer Campbell

   Created:

   Revisions:

   Description:

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.

***********************************************************************/
#include <float.h>
#include <gsl/gsl_statistics.h>

#include "ImageData.h"
#include "otherlibs_include.h"
#include "DiffusionParameters.h"
#include "Directions.h"
#include "S2Interpolator.h"
#include "misc.h"
#include <assert.h>
#include <unistd.h>
#include <string>
//namespace volume_io
//{
//#include "volume_io.h"
//}
/*



figure out how to read in a fixed order eg know which minc dims are which .  can check pretty easily and match xspace, yspace, zspace, time to x,y,z,n: but want more general

can do this with minc (not volume_io): ask what dimensions are and then match them; may be able to get numbers for the dimensions (not sure though...)

for now, assumption is t,z,y,x



//DO; create a constructor that sets all values to a given number.  also a method (there is a minc function to do this I think (???))

note: currently, float-double conversion means that direction cosines get, e.g., 1.32640338655939e-26 instead of 0.  this is fine. should perhaps check for < EPSILON and set to 0. 1.3664357297781e-08.  taking rations would be horrendous!, 0.00078389935885533.  should prob. make the directions cosines be real always

DO: copy constructor


If you want to require step sign to be +, run mincreshape +direction

in WriteMincFile: could write a random-numberbased filename if user chooses not to write: could put it in /tmp and tell the user it is there

also in WriteMincFile: user could choose to write any format e.g short int


*/

/*
//not doing at this point
ImageData::ImageData()
{
  //defaults
  //create a new volume.  don't allocate.

}
*/

/*
ImageData::ImageData(ImageData* template_image_data)
{

//create a new volume: based on that of the template

  //  copy_volume_definition or copy_volume_definition_no_alloc

}
  */

//copy constructor requires const:
/*
ImageData::ImageData(const ImageData &image_data)
{
  //deep copy everything
  //copy_volume
}
*/



//not implemented: probably won't use.  this would do a complete copy.
  //can pass volume efficiently because it is actually a pointer to a struct
/*
ImageData::ImageData(Volume image_minc_volume)
{
  //copy_volume


}
*/

/*
ImageData::ImageData(Volume *image_minc_volume,DiffusionParameters *diffusion_parameters)
{
}
*/

 //DO: need to set 
 //m_image_geom.vector_inverse_xfm=NULL;
 //m_image_geom.point_inverse_xfm=NULL;
 //in all constructors.

ImageData::ImageData():
  _direction_cosines_exist(false),
  _image_array_is_current(false),
  _image_array_exists(false),
  _ODF_max_interpolator_set(false),
  m_sphere_harm_constructed(false),
  m_maxima_array_exists(false),
  m_sigma_theta_array_exists(false),
  _reference_count(1),
  _directions_are_set(false),
  _image_type(Image3D_vanilla),
  m_image_geom_init(false)
  //m_image_geom.vector_inverse_xfm(NULL)
{
}

//make an ImageData directly from the minc file name
ImageData::ImageData(const char *minc_filename):
  _direction_cosines_exist(false),
  _image_array_is_current(false),
  _image_array_exists(false),
  _ODF_max_interpolator_set(false),
  m_sphere_harm_constructed(false),
  m_maxima_array_exists(false),
  m_sigma_theta_array_exists(false),
  _reference_count(1),
  _directions_are_set(false),
  _image_type(Image3D_vanilla),
  m_image_geom_init(false)
  //m_image_geom.vector_inverse_xfm(NULL)
{
  
  LoadMincFile(minc_filename);
}

ImageData::ImageData(IMAGE_TYPE image_type, const IMAGE_GEOM& image_geom):
  _direction_cosines_exist(false),
  _image_array_is_current(false),
  _image_array_exists(false),
  _ODF_max_interpolator_set(false),
  m_sphere_harm_constructed(false),
  m_maxima_array_exists(false),
  m_sigma_theta_array_exists(false),
  _reference_count(1),
  _directions_are_set(false),
  _image_type(Image3D_vanilla),
  m_image_geom_init(false)
  //m_image_geom.vector_inverse_xfm(NULL)
{
//could add 3D option: currently 4D only

  this->Construct(image_type, 
                  image_geom.x_size, image_geom.y_size, image_geom.z_size, image_geom.n_size, 
                  image_geom.x_start, image_geom.y_start, image_geom.z_start, image_geom.n_start, 
                  image_geom.x_step, image_geom.y_step, image_geom.z_step, image_geom.n_step, 
                  image_geom.x_direction_cosine, image_geom.y_direction_cosine, image_geom.z_direction_cosine);
}

ImageData::ImageData(IMAGE_TYPE image_type, 
                     int x_size, int y_size, int z_size, int n_size, 
                     float x_start, float y_start, float z_start, float n_start, 
                     float x_step, float y_step, float z_step, float n_step, 
                     const VECTOR3D& x_direction_cosine,  const VECTOR3D& y_direction_cosine, const VECTOR3D& z_direction_cosine):
  _direction_cosines_exist(false),
  _image_array_is_current(false),
  _image_array_exists(false),
  _ODF_max_interpolator_set(false),
  m_sphere_harm_constructed(false),
  m_maxima_array_exists(false),
  m_sigma_theta_array_exists(false),
  _reference_count(1),
  _directions_are_set(false),  
  _image_type(Image3D_vanilla),
  m_image_geom_init(false)
  //m_image_geom.vector_inverse_xfm(NULL)
{

 Construct(image_type, 
           x_size, y_size, z_size, n_size,
           x_start, y_start, z_start, n_start, 
           x_step, y_step, z_step, n_step,
           x_direction_cosine, y_direction_cosine, z_direction_cosine);
}


ImageData::ImageData(IMAGE_TYPE image_type, 
                     int x_size, int y_size, int z_size, 
                     float x_start, float y_start, float z_start, 
                     float x_step, float y_step, float z_step,
                     const VECTOR3D& x_direction_cosine,  const VECTOR3D& y_direction_cosine, const VECTOR3D& z_direction_cosine):
                     
  _direction_cosines_exist(false),
  _image_array_is_current(false),
  _image_array_exists(false),
  _ODF_max_interpolator_set(false),
  m_sphere_harm_constructed(false),
  m_maxima_array_exists(false),
  m_sigma_theta_array_exists(false),
  _reference_count(1),
  _directions_are_set(false),
  _image_type(Image3D_vanilla)
  //m_image_geom.vector_inverse_xfm(NULL)
{
  _reference_count=1;
  _directions_are_set=false;

  Construct(image_type, 
            x_size, y_size, z_size, 1,
            x_start, y_start, z_start, 0, 
            x_step, y_step, z_step, 0,
            x_direction_cosine, y_direction_cosine, z_direction_cosine);
  _image_type=Image3D_vanilla;      
  _direction_cosines_exist=false;
  _image_array_is_current=false;  
  _image_array_exists=false;
  _ODF_max_interpolator_set=false;
  m_sphere_harm_constructed=false;
  m_maxima_array_exists=false;
  m_sigma_theta_array_exists=false;
  m_image_geom_init=false;
}


void ImageData::LoadMincFile(const char *minc_filename)
{
  int dimensions;
  int cdfid;
  int acquisition_varid;
  
  //first open the minc volume to get the number of dimensions. and check if certain attributes exist. assumes it will be either t,z,y,x (with directions and with (DWIdata) or without (ODFdata) bvalues) or z,y,x (Image3D_vanilla)

  
  try
  {

    //cdfid=miopen(minc_filename,NC_NOWRITE);
    //ncinquire(cdfid,&dimensions,NULL,NULL,NULL);
    minc::minc_1_reader rdr;
    
    rdr.open(minc_filename);
    minc::load_4d_volume(rdr,_image_minc_volume);
    dimensions=rdr.dim_no();
    
    switch (dimensions) 
    {
      case 3:
        _image_type=Image3D_vanilla;
        break;
        
      case 4:
        //acquisition_varid=ncvarid(cdfid,"acquisition"); //not used

        //figure out if there is a bvalues attribute:

        if (this->_bvaluesAttributeExists(rdr))
        {
            _image_type=DWIdata;
        } else if(this->_directionsAttributeExists(rdr)) {
            _image_type=ODFdata;
          }
        /*commented because it doesn't appear that this could return true unless file writeout includes an attribute...
        else if(this->_LabelProbAttributeExists(minc_filename))
          {
            _image_type=ProbLabel;
          }
        */
        else {
            _image_type=Maxima3; 
          }
    }
    _direction_cosines_exist=true;
    
    _x_direction_cosine.x=_image_minc_volume.direction_cosines(0)[0];
    _x_direction_cosine.y=_image_minc_volume.direction_cosines(0)[1];
    _x_direction_cosine.z=_image_minc_volume.direction_cosines(0)[2];
    
    _y_direction_cosine.x=_image_minc_volume.direction_cosines(1)[0];
    _y_direction_cosine.y=_image_minc_volume.direction_cosines(1)[1];
    _y_direction_cosine.z=_image_minc_volume.direction_cosines(1)[2];
    
    _z_direction_cosine.x=_image_minc_volume.direction_cosines(2)[0];
    _z_direction_cosine.y=_image_minc_volume.direction_cosines(2)[1];
    _z_direction_cosine.z=_image_minc_volume.direction_cosines(2)[2];
    
    switch(_image_type)
    {
      case DWIdata:
        _diffusion_parameters = new DiffusionParameters(minc_filename);//TODO:fix this
        m_neighbour_LUT_calculated=false;
        break;
        
      case ODFdata:
        _directions = new Directions(minc_filename);
        _directions_are_set=true;
        m_neighbour_LUT_calculated=false;
    }
    
    
  }  catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    //std::cerr << err.msg()<<std::endl;
    Error(err.msg());
  }
}


//not done: might as well write diffusion parameters to minc header before reading in
/*
ImageData::ImageData(char* minc_filename,char* diffusion_parameters_filename)
{
}
*/


#if 0 //unused, not implemented
ImageData::ImageData(IMAGE_TYPE image_type, int x_size, int y_size, int z_size)
{
  /*
starts 0
steps 1
direction cosines z001 y010 x100
maybe all values init to 0? or not initialized.

   */
  cout << "not implemented\n";
  const char *dimension_names[3]={"zspace","yspace","xspace"}; //IRL const for gcc4
  _reference_count=1;
  _directions_are_set=false;
  m_sphere_harm_constructed=false;
  _image_array_is_current=false;
}


ImageData::ImageData(IMAGE_TYPE image_type, int x_size, int y_size, int z_size, int n_size)
{
  /*
starts 0
steps 1
direction cosines z001 y010 x100
maybe all values init to 0? or not initialized.

   */
  const char *dimension_names[4]={"time","zspace","yspace","xspace"};  //IRL const ofr gcc4
  _reference_count=1;
  cout << "not implemented\n";
  _directions_are_set=false;
  m_sphere_harm_constructed=false;
  m_image_geom_init(false);
}



#endif // unused


void ImageData::Construct(IMAGE_TYPE image_type, 
                          int x_size, int y_size, int z_size, int n_size, 
                          float x_start, float y_start, float z_start, float n_start, 
                          float x_step, float y_step, float z_step, float n_step, 
                          const VECTOR3D& x_direction_cosine,  const VECTOR3D& y_direction_cosine, const VECTOR3D& z_direction_cosine)
{

  int i;
  _reference_count=1;
  _directions_are_set=false;

  switch(image_type)
    {
    case DWIdata:

      _image_type=DWIdata;
      //currently no preset values
      //flags:
      //now create DiffusionParameters:

      _diffusion_parameters=new DiffusionParameters(n_size);
      m_neighbour_LUT_calculated=false;
      _directions_are_set=false;
      break;

    case ODFdata:


      _image_type=ODFdata;
      //currently no preset values
      //flags

      _directions=new Directions(n_size);
      _directions_are_set=true;
      m_neighbour_LUT_calculated=false;

      break;

    case ProbLabel:

    	_image_type=ProbLabel;

    	// For Problabel data, it is already known what n_size(length of 4th dimension) is. (Maxima and sigmas associated with each voxel)
    	n_size= 36;//JC changed to 36 for now (non-fanning case).  need to input depending on whether or not there is fanning
        break;
        
    case Maxima3:

      _image_type=Maxima3;

      //currently no preset values

      //flags


      break;
    }
  //
  _image_minc_volume.resize(x_size,y_size,z_size,n_size);//simple
  _image_minc_volume.start()=minc::IDX<double>(x_start,y_start,z_start);
  _image_minc_volume.t_start()=n_start;
  _image_minc_volume.step()=minc::IDX<double>(x_step,y_step,z_step);
  _image_minc_volume.t_step()=n_step;
  
  _image_minc_volume.direction_cosines(0)=minc::IDX<double>(x_direction_cosine.x,x_direction_cosine.y,x_direction_cosine.z);
  _image_minc_volume.direction_cosines(1)=minc::IDX<double>(y_direction_cosine.x,y_direction_cosine.y,y_direction_cosine.z);
  _image_minc_volume.direction_cosines(2)=minc::IDX<double>(z_direction_cosine.x,z_direction_cosine.y,z_direction_cosine.z);

  //_direction_cosines_exist=false;//bullshit
  _direction_cosines_exist=true;
  
  _x_direction_cosine=x_direction_cosine;
  _y_direction_cosine=y_direction_cosine;
  _z_direction_cosine=z_direction_cosine;
  
  _image_array_is_current=false;  
  _image_array_exists=false;
  _ODF_max_interpolator_set=false;
  m_sphere_harm_constructed=false;
  m_maxima_array_exists=false;
  m_sigma_theta_array_exists=false;
}



ImageData::~ImageData()
{

  int x,y,z;
  switch (_image_type)
  {
      
    case DWIdata:
      _diffusion_parameters->Release();
      if (_directions_are_set)
	{
	  _directions->Release();
	}
      break;


    case ODFdata:
      if (_directions_are_set)
	{
	  _directions->Release();
	}
      break;

    case Image3D_vanilla:
      //nothing
      break;

    case Maxima3:
      //nothing
      break;
    case ProbLabel:
      break; 
  }
  
  if (_image_array_exists)
  {
      gsl_vector_float_free(_image_array);
  }

  if (_ODF_max_interpolator_set)
  {
      delete(_ODF_max_interpolator);
  }

  if (m_maxima_array_exists==true) //there is a problem here: m_maxima_array_exists goes from false to an integer (40, 248, 62...) and shouldn't...
  {
      for (x=0;x<this->GetXSize();x++)
	{
	 for (y=0;y<this->GetYSize();y++)
	   {
	     for (z=0;z<this->GetZSize();z++)
	       {
		 if (m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated)
		   {
		     m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima->Release();
		   }
	       }
	   }
	}
      free(m_maxima_array);
  }

  if (m_sigma_theta_array_exists==true)
    {
      for (x=0;x<this->GetXSize();x++)
	{
	 for (y=0;y<this->GetYSize();y++)
	   {
	     for (z=0;z<this->GetZSize();z++)
	       {
		 if (m_sigma_theta_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated==true)
		   {
		     free(m_sigma_theta_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].sigma_theta);
		   }
	       }
	   }
	}
      free(m_sigma_theta_array);
    }
  
  if (_image_type==DWIdata || _image_type==ODFdata)
    {
  NeighbourLUTEntry *new_entry;
  if (m_neighbour_LUT_calculated==true && m_neighbour_LUT!=NULL)
    {
      for (int i=0; i<this->GetNSize(); i++)
	{
	  while (m_neighbour_LUT[i]!=NULL)

	    {
	      new_entry = m_neighbour_LUT[i];
	      m_neighbour_LUT[i] = new_entry->next;
	      free(new_entry);
	    }
	}
      free(m_neighbour_LUT);
    }
    }


  if (m_sphere_harm_constructed)
    {
      delete(m_sphere_harm);
      m_SH_fit->Release();
    }

 

}

/*
bool ImageData::WriteMincFile(char* minc_filename)
{


    write the minc file (could use volume_io)


    if this->DiffusionParametersExist()
    write diffusion parameters to the minc file header

    miattputdbl




}
*/

float ImageData::GetValue(int x,int y,int z)
{
  return   _image_minc_volume.get(x,y,z,0);
}



float ImageData::GetValue(int x,int y,int z,int n)
{
  return   _image_minc_volume.get(x,y,z,n);
}

void ImageData::SetValue(int x,int y,int z,float value)
{
  _image_minc_volume.set(x,y,z,0,value);
  _image_array_is_current=false;
}


void ImageData::SetValue(int x,int y,int z,int n,float value)
{
  _image_minc_volume.set(x,y,z,n,value);
  _image_array_is_current=false;

}

int ImageData::GetDimensions()
{
  return _image_minc_volume.frames()>1?4:3;
}

//the arrays for the number of (sizes,starts,steps,etc) are much larger than the dimensions we expect (currently 10)
//if 3D or 4D, the indices corresponding to zyx change.

int ImageData::GetXSize()
{
  return _image_minc_volume.dim(0);
}

int ImageData::GetYSize()
{
  return _image_minc_volume.dim(1);
}


int ImageData::GetZSize()
{
  return _image_minc_volume.dim(2);
}

int ImageData::GetNSize()
{
  return _image_minc_volume.frames();
}

double ImageData::GetXStep()
{
  return _image_minc_volume.step()[0];
}

double ImageData::GetYStep()
{
  return _image_minc_volume.step()[1];
}

double ImageData::GetZStep()
{
  return _image_minc_volume.step()[2];
}

double ImageData::GetNStep()
{
  return _image_minc_volume.t_step();
}


double ImageData::GetXStart()
{
  return _image_minc_volume.start()[0];
}

double ImageData::GetYStart()
{
  return _image_minc_volume.start()[1];
}

double ImageData::GetZStart()
{
  return _image_minc_volume.start()[2];
}

double ImageData::GetNStart()
{
  return _image_minc_volume.t_start();
}

int ImageData::GetXStepSign() //Hmmmm very thorough
{
  if (this->GetXStep()<0)
    {
      return -1;

    }
  else
    {
      return 1;
    }

}

int ImageData::GetYStepSign()
{
  if (this->GetYStep()<0)
    {
      return -1;

    }
  else
    {
      return 1;
    }

}

int ImageData::GetZStepSign()
{
  if (this->GetZStep()<0)
    {
      return -1;

    }
  else
    {
      return 1;
    }

}

int ImageData::GetNStepSign()//irrelevant period
{
  if (this->GetNStep()<0)
    {
      return -1;

    }
  else
    {
      return 1;
    }

}

VECTOR3D ImageData::GetXDirectionCosine()
{
  if (!_direction_cosines_exist)
    Error("Direction cosines not set!");
  return _x_direction_cosine;
}

VECTOR3D ImageData::GetYDirectionCosine()
{
  if (!_direction_cosines_exist)
    Error("Direction cosines not set!");
  return _y_direction_cosine;
}

VECTOR3D ImageData::GetZDirectionCosine()
{
  if (!_direction_cosines_exist)
    Error("Direction cosines not set!");
  
  return _z_direction_cosine;
}

VECTOR3D ImageData::GetVoxelFromWorldCoordinates(double x_world, double y_world, double z_world)
{
  
  _vector voxel=_image_minc_volume.frame(0).world_to_voxel_c(minc::IDX<double>(x_world,y_world,z_world));
  
  VECTOR3D VOXEL;
  VOXEL.x = voxel[0];
  VOXEL.y = voxel[1];
  VOXEL.z = voxel[2];
  return VOXEL;      
}



DiffusionParameters *ImageData::GetDiffusionParameters()
{
  //check if DWIdata
  return _diffusion_parameters;

}


void ImageData::SetDiffusionParameters(DiffusionParameters *diffusion_parameters)
{
  //DO: check if DWIdata
  //if it is, that means it has diffusion parameters already, we delete the one it has

  _diffusion_parameters->Release();
  _diffusion_parameters=diffusion_parameters;
  _diffusion_parameters->Retain();

}

Directions *ImageData::GetDirections()
{
  switch (_image_type)
    {
    case ODFdata:
      return _directions;
      break;

    case DWIdata:
      return this->GetDiffusionParameters()->GetDirections();
      break;
    }

}

void ImageData::SetDirections(Directions *directions)
{

  int i;
  if (_image_type==ODFdata)
    {

      //DO: only release if had some before!

      if (_directions_are_set)
	{
	  _directions->Release();
	}


      _directions=directions;
      _directions->Retain();

      _directions_are_set=true;
    }

  if (_image_type==DWIdata)
    {

      this->GetDiffusionParameters()->SetDirections(directions);

    }


}
  template<class T> void generate_info_alternative(const minc::simple_4d_volume<T>& vol,minc::minc_info& info)
  {
     bool have_time=vol.frames()>1||vol.t_step()!=0.0; //assume that it is 3D file otherwise
     
     bool is_vector=false;
      
      if(typeid(T)==typeid(minc::fixed_vec<3,float>)) {
        is_vector=true;
      } 
      
      info.resize(3+(is_vector?1:0)+(have_time?1:0));
      
      if(is_vector)
      {
        info[0].dim=minc::dim_info::DIM_VEC;
        info[0].length=3;
        info[0].step=1;
        
      }
      
      if(have_time) 
      {
        info[(is_vector?1:0)].dim=minc::dim_info::DIM_TIME;
        info[(is_vector?1:0)].step=vol.t_step();
        info[(is_vector?1:0)].start=vol.t_start();
        info[(is_vector?1:0)].length=vol.frames();
      }
          
      //All this hanky-panky because Ilana says so
      for(int i=0;i<3;i++)
      {
        info[(2-i)+(is_vector?1:0)+(have_time?1:0)].dim=minc::dim_info::dimensions( minc::dim_info::DIM_X+i);
  
        info[(2-i)+(is_vector?1:0)+(have_time?1:0)].length=vol.dim(i);
        info[(2-i)+(is_vector?1:0)+(have_time?1:0)].step=vol.step()[i];
        info[(2-i)+(is_vector?1:0)+(have_time?1:0)].start=vol.start()[i];
        
        info[(2-i)+(is_vector?1:0)+(have_time?1:0)].have_dir_cos=true;
        for(int j=0;j<3;j++)
          info[(2-i)+(is_vector?1:0)+(have_time?1:0)].dir_cos[j]=vol.direction_cosines(i)[j];//Ilana also says there was a bug here
      }
          
  }


//could return true or false depending whether write or no
//this assumes _diffusion_parameters, directions, etc, exist
//max length of string should not be hardcoded...
void ImageData::WriteMincFile(const char *minc_filename, bool clobber)
{

  bool may_output=true;
  char use_new_filename;
  char new_minc_filename[1000];


  //makes the buffer bigger:
  strcpy(new_minc_filename,minc_filename);


  while (may_output)
  {
      if (clobber || access(new_minc_filename, F_OK))  { 
          minc::minc_info _info;
          generate_info_alternative(_image_minc_volume,_info);
          //generate_info(_image_minc_volume,_info);
          minc::minc_1_writer wrt;
          
          wrt.open(new_minc_filename,_info,2,NC_FLOAT,true);
          
          //output_volume(new_minc_filename,NC_FLOAT,TRUE,(Real)this->GetMinimumValue(),(Real)this->GetMaximumValue(),_image_minc_volume,NULL,NULL);
	  
	  
	  switch(_image_type) 
	  {
	    case DWIdata:      
	      //this->_diffusion_parameters->WriteToMincFile(new_minc_filename);
              
	      break;
	    case ODFdata:
	      //this->_directions->WriteToMincFile(new_minc_filename); //IRL why is this commented?
              this->_directions->WriteToMincFile(wrt); //IRL add this
	      break;
	    case Image3D_vanilla:
	      //nothing more
	      break;
	    case ProbLabel://write the acquisition attribute: we aren't using this variable
	      break;
	    case Maxima3:
	      //nothing more
	      break;

	  }
       minc::save_4d_volume(wrt,_image_minc_volume);
       break;
          
      } else  {//prompt user for new filename
        cout << "WriteMincFile:"<<new_minc_filename<<" exists!"<<std::endl;
        cout << "Do you wish to use an alternate file name? [y] [n]"<<std::endl;
      cin >> use_new_filename;
      if (use_new_filename=='y')
      {
        cout << "new file name:"<<std::endl;
	      cin >> new_minc_filename;
      } else {
        cout << "Overwriting..."<<std::endl;
	clobber = true;
	      
      }
      }
    }
}


void ImageData::WriteMincFile(const char* minc_filename, bool clobber, const char *history_string, bool bbox)
{
  //TODO convert this to new style, right now using simple interface
  WriteMincFile(minc_filename,clobber);


  //write the history to this file:
  char command[1000];
  //cout << "write history string " << history_string << endl;
  sprintf(command,"minc_modify_header -sinsert :history=\"%s\" %s\n",history_string,minc_filename);
  //cout << command << endl;
  system(command);


   /*
  bool may_output=true;
  char use_new_filename;
  char new_minc_filename[1000];

  //makes the buffer bigger:
  strcpy(new_minc_filename,minc_filename);

    // IL want temp space
  char tempdirname[200];
  char tempfilename[200];
  char command[1000];
  char output[50];
  //set_volume_real_range(_image_minc_volume,(Real)this->GetMinimumValue(),(Real)this->GetMaximumValue());
  



  while (may_output)
  {
      
      if (!clobber)
      {
	  clobber=!access(new_minc_filename,F_OK);
      }
  
      if (clobber)
      {
	 if(bbox) //IL to save smallest bounding box to new_minc_filename
	 {
	  //IL create temporary directory and file name
	  sprintf(tempdirname,"/tmp/tmp%i",(int)time(NULL));
	  sprintf(command,"mkdir %s\n",tempdirname);
          system(command);
	  sprintf(tempfilename,"%s/tmp%f.mnc",tempdirname,(double)rand());
	  output_volume(tempfilename,NC_FLOAT,TRUE,(Real)this->GetMinimumValue(),(Real)this->GetMaximumValue(),_image_minc_volume,history_string,NULL);

	  //If volume is empty, save volume of size 0
	  if(this->GetMaximumValue()==0)
	  {
		sprintf(command,"mincresample -clobber -nelements 1 1 1 %s %s\n",tempfilename,new_minc_filename);
		system(command);
	  }
	  else
	  {
		sprintf(command,"mincbbox -mincreshape %s\n",tempfilename); // get bounding box
		FILE *fp;
		char c;
		fp = popen(command, "r");
		fgets (output, 50, fp );

		sprintf(command,"mincreshape -clobber %s %s %s\n",tempfilename,new_minc_filename,output); //use mincreshape to get smallest possible file
		system(command);
	  }
	   // clean up
          sprintf(command,"rm -rf %s\n",tempdirname);
          system(command);
         }
	 else //IL save entire volume to new_minc_filename
	 {	
	   output_volume(new_minc_filename,NC_FLOAT,TRUE,(Real)this->GetMinimumValue(),(Real)this->GetMaximumValue(),_image_minc_volume,history_string,NULL);
         }

	  switch(_image_type)
	    {
	    case DWIdata:
	      this->_diffusion_parameters->WriteToMincFile(new_minc_filename);
	      break;
	    case ODFdata:
	      this->_directions->WriteToMincFile(new_minc_filename);
	      break;
	    case Image3D_vanilla:
	      //nothing more
	      break;
	    case ProbLabel://write the acquisition attribute: we aren't using this variable
	      break;
	    case Maxima3:
	      //nothing more
	      break;
	    }
	  may_output=false;

	}
      else //prompt user for new filename
	{
	  cout << "Do you wish to use an alternate file name? [y] [n]\n";
	  cin >> use_new_filename;
	  if (use_new_filename=='y')
	    {
	      cout << "new file name:\n";
	      cin >> new_minc_filename;
	    }
	  else
	    {
	      may_output=false;

	    }


	}
	
    }*/
 

}

void ImageData::WriteMincFileNoRange(const char* minc_filename, bool clobber, const char *history_string)
{
  WriteMincFile(minc_filename,clobber);
  //TODO: figure out what was needed 
  /*
  bool may_output=true;
  char use_new_filename;
  char new_minc_filename[1000];

  //makes the buffer bigger:
  strcpy(new_minc_filename,minc_filename);

  //this was wrecking the result: ImageData values changed! need to check how storage is happening and what's up
  //set_volume_real_range(_image_minc_volume,(Real)this->GetMinimumValue(),(Real)this->GetMaximumValue());
  //however, the result can't be visualized with Display or register if we don't do this...


   while (may_output)
    {

      if (!clobber)
	{
	  clobber=check_clobber_file(new_minc_filename);
	}

      if (clobber)
	{

	  output_volume(new_minc_filename,NC_FLOAT,TRUE,(Real)this->GetMinimumValue(),(Real)this->GetMaximumValue(),_image_minc_volume,history_string,NULL);


	  switch(_image_type)
	    {
	    case DWIdata:
	      this->_diffusion_parameters->WriteToMincFile(new_minc_filename);
	      break;
	    case ODFdata:
	      this->_directions->WriteToMincFile(new_minc_filename);
	      break;
	    case Image3D_vanilla:
	      //nothing more
	      break;
	    case ProbLabel://write the acquisition attribute: we aren't using this variable
	      break;
	    case Maxima3:
	      //nothing more
	      break;
	    }
	  may_output=false;

	}
      else //prompt user for new filename
	{
	  cout << "Do you wish to use an alternate file name? [y] [n]\n";
	  cin >> use_new_filename;
	  if (use_new_filename=='y')
	    {
	      cout << "new file name:\n";
	      cin >> new_minc_filename;
	    }
	  else
	    {
	      may_output=false;

	    }


	}
    }
  */
  
  //write the history to this file:
  char command[1000];
  //cout << "write history string " << history_string << endl;
  sprintf(command,"minc_modify_header -sinsert :history=\"%s\" %s\n",history_string,minc_filename);
  //cout << command << endl;
  system(command);

}






IMAGE_TYPE ImageData::GetImageType()
{
  return _image_type;

}

//note: volume_io functions for getting max and min [get_volume_real_max(_image_minc_volume)]don't actualy calculate them: just return variables.  only reason this matters is that programs like Display and register design the lookup tables based on these values, and apparently don't reculate (xdisp recalculates)

//using gsl array is really just an exercise at this point: could keep track of min and max and every time we ->SetValue, update (but that might be even slower...)

float ImageData::GetMaximumValue()
{
  if (_image_array_is_current)
    {
      return gsl_vector_float_max(_image_array);

    }
  else
    {
      this->_PutVolumeInArray();
      return gsl_vector_float_max(_image_array);
    }


}



//this finds max less than exclude
float ImageData::GetMaximumValue(float exclude)
{
  float maxval=FLT_MIN;
  int i,x,y,z;
  float val;

  assert(_image_type==Image3D_vanilla ||  _image_type==Image3D_with_diffusion_parameters);

  //for(i=0;i<this->GetNSize();i++)
  //{
      for (x=0;x<this->GetXSize();x++)
	{
	  for (y=0;y<this->GetYSize();y++)
	    {
	      for (z=0;z<this->GetZSize();z++)
		{
		  val=this->GetValue(x,y,z);
		  if (val>maxval && val<exclude)
		    {
		      maxval=val;
		    }
		}
	    }
	}
	  //}


      return maxval;
}


float ImageData::GetMaximumValue(int x,int y, int z)
{
  int i;
  float max=FLT_EPSILON;


  if (this->GetDimensions()!=4)
    {
      Error("This function is for 4D data only");//IRL Error for gcc4

    }


  for(i=0;i<this->GetNSize();i++)
    {
      if (this->GetValue(x,y,z,i)>max)
	{
	  max=this->GetValue(x,y,z,i);

	}

    }

  return max;


}


float ImageData::GetMaximumValue(int x, int y, int z, int &index)
{
  int i;
  float max=-FLT_EPSILON;
  int index_c;

  if (this->GetDimensions()!=4)
    {
      Error("This function is for 4D data only");

    }

  for(i=0;i<this->GetNSize();i++)
    {
      if (this->GetValue(x,y,z,i)>max)
	{
	  max=this->GetValue(x,y,z,i);
	  index_c=i;
	}
    }

  index=index_c;
  return max;
}

float ImageData::GetMaximumValue(int x, int y, int z, int &index, int start, int stop)
{
  int i;
  float max=-FLT_EPSILON;
  int index_c;

  if (this->GetDimensions()!=4)
    {
      Error("This function is for 4D data only");

    }

  for(i=start;i<stop;i++)
    {
      if (this->GetValue(x,y,z,i)>max)
	{
	  max=this->GetValue(x,y,z,i);
	  index_c=i;
	}
    }

  index=index_c;
  return max;
}


float ImageData::GetNumberOfVoxelsAbove(float threshold)
{
  float cnt=0;
  int i,x,y,z;
  
  for(i=0;i<GetNSize();i++)
  {
    for(z=0;z<GetZSize();z++)
      for(y=0;y<GetYSize();y++)
         for(x=0;x<GetXSize();x++)
          if(GetValue(x,y,z,i)>threshold ) cnt+=1;
  }
  return cnt;
}

float ImageData::GetMinimumValue()
{
  if (_image_array_is_current)
    {
      return gsl_vector_float_min(_image_array);

    }
  else
    {
      this->_PutVolumeInArray();
      return gsl_vector_float_min(_image_array);
    }

}

float ImageData::GetMinimumNonZeroValue(int x,int y, int z)
{
  int i;
  float min=FLT_MAX;


  if (this->GetDimensions()!=4)
    {
      Error("This function is for 4D data only"); //IRL Error for gcc4

    }


  for(i=0;i<this->GetNSize();i++)
    {
      if (this->GetValue(x,y,z,i)<min && this->GetValue(x,y,z,i)>FLT_EPSILON)
	{
	  min=this->GetValue(x,y,z,i);

	}

    }

  return min;


}

float ImageData::GetMinimumNonZeroValue(int x, int y, int z, int &index)
{
  int i;
  float min=FLT_MAX;
  int index_c;

  if (this->GetDimensions()!=4)
    {
      Error("This function is for 4D data only");

    }


  for(i=0;i<this->GetNSize();i++)
    {
      if (this->GetValue(x,y,z,i)<min && this->GetValue(x,y,z,i)>FLT_EPSILON)
	{
	  min=this->GetValue(x,y,z,i);
	  index_c=i;
	}

    }

  index=index_c;
  return min;
}


void ImageData::_PutVolumeInArray()
{
  int x,y,z,n;
  switch (this->GetDimensions())
    {
    case 3:
      _image_array=gsl_vector_float_alloc(this->GetXSize()*this->GetYSize()*this->GetZSize());

      int x,y,z;
      for(x=0;x<this->GetXSize();x++)
	{
	  for(y=0;y<this->GetYSize();y++)
	    {
	      for(z=0;z<this->GetZSize();z++)
		{
		  gsl_vector_float_set(_image_array,z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x,this->GetValue(x,y,z));
		}
	    }
	}

	break;

    case 4:
      _image_array=gsl_vector_float_alloc(this->GetXSize()*this->GetYSize()*this->GetZSize()*this->GetNSize());

      for(x=0;x<this->GetXSize();x++)
	{
	  for(y=0;y<this->GetYSize();y++)
	    {
	      for(z=0;z<this->GetZSize();z++)
		{
		  for(n=0;n<this->GetNSize();n++)
		    {
		      gsl_vector_float_set(_image_array,n*this->GetZSize()*this->GetYSize()*this->GetXSize()+z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x,this->GetValue(x,y,z,n));
		    }
		}
	    }
	}

      break;
    }

  _image_array_exists=true;
  _image_array_is_current=true;


}





/*
  Parya has added this method to check for the IMAGE_TYPE in the case where a minc file is read which contatins probabilistic labeling info.
*/

bool ImageData::_LabelProbAttributeExists(char *minc_filename)
{
  char command[300];
  char tempfilename[200];
  int problabel_exist;


  sprintf(tempfilename,"/tmp/tmp%i.txt",(int) time(NULL));


  sprintf(command,"mincheader %s | grep acquisition:problabel | wc -w > %s\n",minc_filename,tempfilename);
  system(command);


  FILE * fd;
  fd = fopen(tempfilename,"r");
  fscanf(fd,"%d\n",&problabel_exist);
  fclose(fd);

  sprintf(command,"rm -f %s\n",tempfilename);

  system(command);


  if (problabel_exist>0)
    {
      return true;
    }
  else
    {
      return false;

    }
}

/*
  Here ends the method added by Parya.
*/

bool ImageData::_bvaluesAttributeExists(minc::minc_1_reader& rdr)
{
  //not very elegant, but easy to type :)
  std::vector<double> bvalue=rdr.att_value_double("acquisition","bvalues");
  return !bvalue.empty();

}


bool ImageData::_directionsAttributeExists(minc::minc_1_reader& rdr)
{
  //not very elegant, but easy to type :)
  std::vector<double> dir_x=rdr.att_value_double("acquisition","direction_x");
  return !dir_x.empty();
}

void ImageData::SetAllVoxels(float value)
{
  int x,y,z;

  switch (this->GetDimensions())
    {
    case 3:
      for (x=0;x<this->GetXSize();x++)
	{
	  for (y=0;y<this->GetYSize();y++)
	    {
	      for (z=0;z<this->GetZSize();z++)
		{
		  this->SetValue(x,y,z,value);

		}
	    }
	}
      break;

    }

}

void ImageData::Retain()
{
  _reference_count+=1;

}

void ImageData::Release()
{
  if (_reference_count<2)
    {
      delete(this);

    }
  else
    {
      _reference_count-=1;
    }
}



//this is faster but maxima are sometimes a bit off (constrained to be original sampling directions)
Directions *ImageData::GetODFMaximaFast(int x,int y,int z)
{

  //for now, the ODF is assumed to be blurred in such a way that the *centre* of a broad maximum (due e.g. to fanning or curvature) is selected.
  //CreateNeighbourLUT determines how big a region is considered.  Has only been tested for QBI data (i.e. certain finite angular resolution expected) and N~100
  //blurring should be implemented for more noisy ODFs.  if called from GetODFMaxima, the SH fit has been done so noise should be taken care of

  //this has only been tested for N~100 because same for CreateNeighbourLUT()

  //tolerance isn't completed so is set to 0.  need to deal with large areas and determine the right tolerance number and fix seg fault...



  //test: is it this method that explodes memory? no
  /*
  Directions *dirs=new Directions(1);
  dirs->SetXComponent(0,0.5);
  dirs->SetYComponent(0,0.5);
  dirs->SetZComponent(0,0.5);

  return dirs;
  */



  m_maximum_cutoff_fraction=0.5;//this is empirical and should become user-input - not implemented yet


  int i;

  int xx,yy,zz;

  Directions *directions;
  float maxmax;

  VECTOR3D maxmaxdir;
  bool maxmax_changed;
  bool maxmax_passed;

  bool pass;


  float this_value;
  bool first_dir;
  NeighbourLUTEntry *this_entry;



  Vector3D *close_dir;
  int close;
  int ii;
  float tolerance=0;
  Vector3D avg_dir;
  if (!(_image_type==DWIdata || _image_type==ODFdata))
    {
      Error("ImageData is not an ODF");
    }


  //set the array of ODF maxima that can be used if this function has been called for this voxel:

  if (!m_maxima_array_exists)
    {
      m_maxima_array=(ODF_MAXIMA *) malloc(this->GetXSize()*this->GetYSize()*this->GetZSize()*sizeof(ODF_MAXIMA));


      for (xx=0;xx<this->GetXSize();xx++)
	{
	  for (yy=0;yy<this->GetYSize();yy++)
	    {
	      for (zz=0;zz<this->GetZSize();zz++)
		{

		  m_maxima_array[zz*this->GetYSize()*this->GetXSize()+yy*this->GetXSize()+xx].calculated=false;
		}
	    }
	}
      m_maxima_array_exists=true;
    }


  if (m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated)
    {
      m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima->Retain();
      return m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima;
    }


  if (!m_neighbour_LUT_calculated)
    {
      this->CreateNeighbourLUT();
    }

  //for each direction, check whether all neighbours are less.  if so, this is a maximum,
  //check for whether the maximum is very close to neighbour; if so, take intermediate.
  //DO: propagate this through many points so can identify a clump or line or similar values and take an intermediate direction.

  //DO: find out what max #neighbours is...  right now, this is just way too big

  close_dir=(Vector3D *) malloc(this->GetNSize()*sizeof(Vector3D));

  directions=new Directions(1);
  first_dir=true;
  maxmax=FLT_EPSILON;
  maxmax_passed=true;
  int *indices_of_maxima=(int *) malloc(this->GetNSize()*sizeof(int));

  for (i=0;i<this->GetNSize();i++)//for all directions
    {
      close=0;
      pass=true;
      this_entry=m_neighbour_LUT[i];

      if (this->GetValue(x,y,z,i)>maxmax)
	{
	  maxmax=this->GetValue(x,y,z,i);
	  maxmaxdir=this->GetDirections()->GetDirection(i);
	  maxmax_changed=true;
	  maxmax_passed=false;
	}
      else
	{
	  maxmax_changed=false;
	}

      if (m_neighbour_LUT[i]==NULL)
	{
	  pass=false;
	}

      while (pass && i<this->GetNSize())//the second arg is because some weird memory issue causes i to be equal and hence to crash the program...DO: fix this...
	{

	  if (this->GetValue(x,y,z,i)<this->GetValue(x,y,z,this_entry->index))

	    {
	      pass=false;
	    }
	  else if (this->GetValue(x,y,z,i)-this->GetValue(x,y,z,this_entry->index)<tolerance && this->GetValue(x,y,z,i)-this->GetValue(x,y,z,this_entry->index)>0)
	    {
	      close_dir[close]=this->GetDirections()->GetDirection(this_entry->index);
	      close+=1;
	    }
	  if (this_entry->next==NULL)
	    {
	      break;
	    }
	  else
	    {
	      this_entry=this_entry->next;

	    }

	}


      if (pass)
	{
	  if (maxmax_changed)
	    {
	      maxmax_passed=true;
	    }
	  if (!first_dir)
	    {
	      directions->IncreaseNumberOfDirections(1);
	    }
	  else//it was the first dir
	    {
	      first_dir=false;
	    }

	  if (close>0)//take intermediate:
	    {
	      cout << "close " << close << endl;

	      avg_dir.x=this->GetDirections()->GetDirection(i).x;
	      avg_dir.y=this->GetDirections()->GetDirection(i).y;
	      avg_dir.z=this->GetDirections()->GetDirection(i).z;

	      for (ii=0;ii<close;ii++)
		{
		  avg_dir.x+=close_dir[close].x;
		  avg_dir.y+=close_dir[close].y;
		  avg_dir.z+=close_dir[close].z;

		}
	      avg_dir.x=avg_dir.x/(close+1.0);
	      avg_dir.y=avg_dir.y/(close+1.0);
	      avg_dir.z=avg_dir.z/(close+1.0);

	      directions->SetDirection(directions->GetNumberOfDirections()-1,avg_dir);
	    }
	  else
	    {
	      directions->SetDirection(directions->GetNumberOfDirections()-1,this->GetDirections()->GetDirection(i));
	      indices_of_maxima[directions->GetNumberOfDirections()-1]=i;
	    }

	}

    }

  //if we didn't get anything, add maxmax.  (I don't think this should ever happen)

  if (first_dir)
    {
      directions->SetDirection(0,maxmaxdir);
      //cout << "maxmax: " << maxmax << endl;
    }
  else if (!maxmax_passed)
    {
      directions->IncreaseNumberOfDirections(1);
      directions->SetDirection(directions->GetNumberOfDirections()-1,maxmaxdir);
    }


  if (1)
    {
      //Now go through the maxima and remove those that have value a small fraction of the largest maximum: this applies to deconvolution where there might be ringing, but is currently always run:



      first_dir=true;
      int count=0;
      Directions *cut_directions=new Directions(1);

      for (i=0;i<directions->GetNumberOfDirections();i++)
	{
	  if (this->GetValue(x,y,z,indices_of_maxima[i])>m_maximum_cutoff_fraction*maxmax)
	    {
	      if (!first_dir)
		{

		  cut_directions->IncreaseNumberOfDirections(1);

		}
	      else//it was the first dir
		{
		  first_dir=false;
		}
	      cut_directions->SetDirection(count,directions->GetDirection(i));
	      count++;
	    }
	}

      directions->Release();
      directions=cut_directions;


    }//cutoff threshold

  m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima=directions;
  //m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima->Retain();
  m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated=true;

  free(close_dir);
  free(indices_of_maxima);

  //just in case:
  directions->Normalize();

  //cout << directions->GetNumberOfDirections();
  return directions;

}



//note; for this method, the calling program MUST retain the directions if they will be modified or deleted later
Directions *ImageData::GetODFMaxima(int x,int y,int z)
{
  //this does SH fit first, then samples it and takes maxima from the sampling directions
  int order=8;
  int xx,yy,zz;

  //  int *test;

  //  test=(int *)malloc(1000*sizeof(int));

  if (!m_maxima_array_exists)
    {

      //test=(int *)malloc(1000*sizeof(int));

      m_maxima_array=(ODF_MAXIMA *) malloc(this->GetXSize()*this->GetYSize()*this->GetZSize()*sizeof(ODF_MAXIMA));


      for (xx=0;xx<this->GetXSize();xx++)
	{
	  for (yy=0;yy<this->GetYSize();yy++)
	    {
	      for (zz=0;zz<this->GetZSize();zz++)
		{

		  m_maxima_array[zz*this->GetYSize()*this->GetXSize()+yy*this->GetXSize()+xx].calculated=false;
		}
	    }
	}
      m_maxima_array_exists=true;
    }


  if (m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated)
    {

      //m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima->Retain();
      return m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima;
    }


  if (!m_sphere_harm_constructed)
    {
      m_sphere_harm=new SphericalHarmonics();
      m_sphere_harm->SetOrder(order);
      m_sphere_harm_constructed=true;
      m_SH_fit=new ImageData(ODFdata,GetImageGeom());
    }

  m_sphere_harm->CalculateSHFit(this,m_SH_fit,x,y,z);



  m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima=m_SH_fit->GetODFMaximaFast(x,y,z);
   //test:need to retain so that can delete m_SH_fit...
  m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima->Retain();
  m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated=true;



  return m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima;

  //so what happens here is that the SH fit ImageData calls its own GetODFMaxima.  should only get called if they haven't been calculated, because otherwise they were assigned to this ImageData's maxima.


}


Directions *ImageData::GetODFMaximaHR(int x,int y,int z)
{
  //this calculates ODF as done in ODFPlotter for maxima consideration


  float angle_step=0.1745; //empirical; can look at changing: may lose info for some N
  int counter;
  float theta,phi;
  int i,j,k;
  int xx,yy,zz;
  //float neighbour_theta,neighbour_phi;
  //gsl_vector_float *resampled_ODF;
  Directions *directions;
  float maxmax;

  VECTOR3D maxmaxdir;
  bool maxmax_changed;
  bool maxmax_passed;

  bool pass;


  float this_value;
  bool first_dir;

  bool last_passed;

  int theta_points,phi_points;



  theta_points=(int) ceil(2*M_PI/angle_step);
  phi_points=(int) ceil(M_PI/angle_step);

  if (!(_image_type==DWIdata || _image_type==ODFdata))
    {
      Error("ImageData is not an ODF");
    }


  //set the array of ODF maxima that can be used if this function has been called for this voxel:

  if (!m_maxima_array_exists)
    {
      m_maxima_array=(ODF_MAXIMA *) malloc(this->GetXSize()*this->GetYSize()*this->GetZSize()*sizeof(ODF_MAXIMA));


      for (xx=0;xx<this->GetXSize();xx++)
	{
	  for (yy=0;yy<this->GetYSize();yy++)
	    {
	      for (zz=0;zz<this->GetZSize();zz++)
		{

		  m_maxima_array[zz*this->GetYSize()*this->GetXSize()+yy*this->GetXSize()+xx].calculated=false;
		}
	    }
	}
      m_maxima_array_exists=true;
    }


  if (m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated)
    {
      //m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima->Retain();
      return m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima;
    }


  if (!m_neighbour_LUT_calculated)//this boolean is used to mean that the ODF is resampled and resampled dirs are stored and everything (in this loop) is done
    {
      //cout << "size: " << ((int) ((ceil(2*M_PI/angle_step)+1)*ceil(2*M_PI/angle_step)*ceil(M_PI/angle_step))) << endl;
      m_ODF_sample_neighbours=(int *) malloc((theta_points+1)*theta_points*phi_points*sizeof(int) );//max size will be to hold all of a row of const. phi in the case of the poles.  the first entry is the number of neighbours, followed by their indices
      m_neighbour_LUT_calculated=true;



      //not sure whether this commits us to using all of these?  shouldn't? just need at least enough, I'm pretty sure
      //m_ODF_sample_dirs=new Directions(theta_points*phi_points);
      m_ODF_sample_dirs=new Directions(theta_points*((int) ceil(phi_points/2.0)+1));



      counter=0;
      //skipping second hemisphere except a bit to be safe: go just to ceil(phi_points/2.0)+1 points
      //this part below: && counter<theta_points*((int) ceil(phi_points/2.0)+1); isn't necessary
      for (phi=0;phi<angle_step*(ceil(phi_points/2.0)+1) && counter<theta_points*((int) ceil(phi_points/2.0)+1);phi+=angle_step)
	{
	  //cout << "phi " << phi << endl;


	  for (theta=0;theta<2*M_PI && counter<theta_points*((int) ceil(phi_points/2.0)+1);theta+=angle_step)
	    {

	      //assign neighbour indices:
	      if (phi<FLT_EPSILON)
		{
		  for (i=0;i<=theta_points;i++)
		    {
		      if (i==0)
			{

			  m_ODF_sample_neighbours[counter*(theta_points+1)+i]=theta_points;//number of neighbours
			}

		      //the whole row before the pole
		      else
			{
			  m_ODF_sample_neighbours[counter*(theta_points+1)+i]=(i-1)+theta_points;
			}

		    }
		}
	      /*
	      else if (phi<((float)(floor(M_PI/angle_step)))*angle_step+FLT_EPSILON && phi > ((float)(floor(M_PI/angle_step)))*angle_step-FLT_EPSILON)
		{
		  for (i=0;i<=ceil(2*M_PI/angle_step);i++)
		    {
		      if (i==0)
			{
			  m_ODF_sample_neighbours[(int) (counter*(ceil(2*M_PI/angle_step)+1)+i)]=ceil(2*M_PI/angle_step);//number of neighbours
			}
		      //the whole row before the pole
		      else
			{
			  //cout << (counter*(ceil(2*M_PI/angle_step)+1)+i) << endl;
			  m_ODF_sample_neighbours[(int) (counter*(ceil(2*M_PI/angle_step)+1)+i)]=(ceil(M_PI/angle_step)-2)*ceil(2*M_PI/angle_step)+(i-1);
			}
		    }
		}
	      */
	      else //using 8-neighbourhood.  could do 24 (need to implement, including all special cases near poles), but evn that won't work well for unfit ODFs - tiny maxima
		{
		  i=0;
		  m_ODF_sample_neighbours[(int) (counter*(theta_points+1)+i)]=8;

		  for (j=-1;j<=1;j++)
		    {
		      for (k=-1;k<=1;k++)
			{
			  if (!(j==0 && k==0))
			    {
			      i++;


			      if (counter%theta_points==0 && k==-1)
				{
				  m_ODF_sample_neighbours[(int) (counter*(theta_points+1)+i)]=counter+j*theta_points+theta_points-1;

				}
			      else if (counter%theta_points==theta_points-1 && k==1)
				{
				  m_ODF_sample_neighbours[(int) (counter*(theta_points+1)+i)]=counter+j*theta_points-theta_points+1;

				}
			      else
				{
				  m_ODF_sample_neighbours[(int) (counter*(theta_points+1)+i)]=counter+j*theta_points+k;

				}

			    }
			}
		    }

		}


	      m_ODF_sample_dirs->SetXComponent(counter, cos(theta)*sin(phi));
	      m_ODF_sample_dirs->SetYComponent(counter, sin(theta)*sin(phi));
	      m_ODF_sample_dirs->SetZComponent(counter, cos(phi));

	      //cout << counter << endl;

	      counter++;
	    }
	}


      //note this is where blurring occurs: play with amount: this is empirical for N~100, our acquisition protocol




      m_ODF_interpolator=new S2Interpolator(this, m_ODF_sample_dirs, 1.5*sqrt(2*M_PI/this->GetNSize()), 1.5*sqrt(2*M_PI/this->GetNSize()), Gaussian);
      //m_ODF_interpolator->CalculateLookupTable();


      //test:
      /*
      counter=0;
      i=0;
      for (phi=0;phi<M_PI;phi+=angle_step)
	{
	  //cout << "phi " << phi << endl;
	  for (theta=0;theta<2*M_PI;theta+=angle_step)
	    {

	      for (i=1;i<m_ODF_sample_neighbours[(int) (counter*(ceil(2*M_PI/angle_step)+1))];i++)
		{
		  cout << "dir: " << m_ODF_sample_dirs->GetXComponent(counter) << " " << m_ODF_sample_dirs->GetYComponent(counter) << " " << m_ODF_sample_dirs->GetZComponent(counter) << " neighbour: " << m_ODF_sample_dirs->GetXComponent(m_ODF_sample_neighbours[(int) (counter*(ceil(2*M_PI/angle_step)+1)+i)]) << " " << m_ODF_sample_dirs->GetYComponent(m_ODF_sample_neighbours[(int) (counter*(ceil(2*M_PI/angle_step)+1)+i)]) << " " << m_ODF_sample_dirs->GetZComponent(m_ODF_sample_neighbours[(int) (counter*(ceil(2*M_PI/angle_step)+1)+i)]) << endl;
		}
	      counter++;
	    }

	}
      */
    }



  //the rest of this is voxel specific:


  //resampled_ODF=gsl_vector_float_alloc((int)(ceil(2*M_PI/angle_step)*ceil(M_PI/angle_step)));

  //resample the ODF: currently done on the fly below
  /*
  counter=0;
  for (phi=0;phi<M_PI;phi+=angle_step)
    {
      for (theta=0;theta<2*M_PI;theta+=angle_step)
	{
	  gsl_vector_float_set(resampled_ODF,counter,m_ODF_interpolator->GetInterpolatedValue(counter,x,y,z));
	  counter++;
	}
    }
  */

  //now go back through all points and test for maxima: for each direction, check whether all neighbours are less.  if so, this is a maximum.

  //DO: to add, if they are close, take intermediate direction, and propagate this through many points so can identify a clump or line or similar values and take an intermediate direction.




  directions=new Directions(1);
  first_dir=true;
  maxmax=FLT_EPSILON;


  //use i for all dirs to check for maximum; j for neighbours of this dir:

  //for (i=0;i<m_ODF_sample_dirs->GetNumberOfDirections();i++)//for all resampled ODF directions
  //for (i=0;i<(ceil(2*M_PI/angle_step))*(ceil(M_PI/angle_step)/2+1);i++)//for all resampled ODF directions up to equator + a row; sample dirs different on other half??
  last_passed=false;
  for (i=0;i<theta_points*((int) ceil(phi_points/2.0));i++)
    //for (i=0;i<(ceil(2*M_PI/angle_step))*3;i++)
    {
      pass=true;
      if (m_ODF_interpolator->GetInterpolatedValue(i,x,y,z)>maxmax)
	{
	  maxmax=m_ODF_interpolator->GetInterpolatedValue(i,x,y,z);
	  maxmaxdir=m_ODF_sample_dirs->GetDirection(i);
	  maxmax_changed=true;
	  maxmax_passed=false;
	}
      else
	{
	  maxmax_changed=false;
	}

      j=1;
      while (pass && j<=m_ODF_sample_neighbours[i*(theta_points+1)])//while nothing bigger found and still looking at neighbours
	{
	  if (m_ODF_interpolator->GetInterpolatedValue(i,x,y,z)<m_ODF_interpolator->GetInterpolatedValue(m_ODF_sample_neighbours[i*(theta_points+1)+j],x,y,z))
	    {
	      //cout << "no pass\n";
	      pass=false;
	    }
	  j++;
	}

      //cout << "end pass: " << pass << endl;
      if (pass && last_passed)
	{
	  cout << "but last passed\n";

	}
      if (pass)
	{
	  if (maxmax_changed)
	    {
	      maxmax_passed=true;
	    }
	  if (!first_dir)
	    {
	      directions->IncreaseNumberOfDirections(1);
	    }
	  else//it was the first dir
	    {
	      first_dir=false;
	    }
	  last_passed=true;

	  directions->SetDirection(directions->GetNumberOfDirections()-1, m_ODF_sample_dirs->GetDirection(i));
	}
      else
	{
	  last_passed=false;
	}

    }

  //if we didn't get anything, add maxmax would, for instance, choose the max in an area of a broad maximum):  also add maxmax if it didn't get added, i.e., every time a maxmax gets changed, note whether it passed, and set a flag
  if (first_dir)
    {
      directions->SetDirection(0,maxmaxdir);
      //cout << "maxmax: " << maxmax << endl;
    }
  else if (!maxmax_passed)
    {
      directions->IncreaseNumberOfDirections(1);
      directions->SetDirection(directions->GetNumberOfDirections()-1,maxmaxdir);
    }

  m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima=directions;
  //m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima->Retain();
  m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated=true;



  return directions;

}






float ImageData::GetODFStdDev(int x,int y,int z)
{

  double *ODF_values;
  ODF_values=(double *) malloc(this->GetNSize()*sizeof(double));
  float std_dev;

  int n;
  for(n=0;n<this->GetNSize();n++)
    {
      ODF_values[n]=(double) this->GetValue(x,y,z,n);

    }

  std_dev = (float) pow(gsl_stats_variance(ODF_values,1,this->GetNSize()),0.5);
  free(ODF_values);
  return std_dev;
}

float ImageData::GetODFMean(int x,int y,int z)
{

  double *ODF_values;
  ODF_values=(double *) malloc(this->GetNSize()*sizeof(double));
  float mean;

  int n;
  for(n=0;n<this->GetNSize();n++)
    {
      ODF_values[n]=(double) this->GetValue(x,y,z,n);

    }

  mean = (float) gsl_stats_mean(ODF_values,1,this->GetNSize());
  free(ODF_values);
  return mean;
}



void ImageData::WriteODFMaxima(const char *filename, bool clobber, char *history_string)
{
  int x,y,z;
  int i;
  ofstream fd(filename);

  //error if ImageData isn't ODF...

  if (!fd)
    {
      Error("cannot create file");
    }



  //assign the maxima, else 0.  maxima written to file only if they are already calculated with GetODFMaxima()

  //write in zyx: (note z is not necessarily fastest varying in the ODF file

  for (x=0;x<this->GetXSize();x++)
    {
      for (y=0;y<this->GetYSize();y++)
	{

	  for (z=0;z<this->GetZSize();z++)

	    {
	      if (m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated)
		{
		  fd << this->GetODFMaxima(x,y,z)->GetNumberOfDirections() << " " ;

		  for (i=0;i<this->GetODFMaxima(x,y,z)->GetNumberOfDirections();i++)
		    {
		      fd << this->GetODFMaxima(x,y,z)->GetXComponent(i) << " " ;
		      fd << this->GetODFMaxima(x,y,z)->GetYComponent(i)<< " ";
		      fd << this->GetODFMaxima(x,y,z)->GetZComponent(i)<< " ";


		    }

		}
	      else //write all 0: confused because this code ssused to write 10 0s, should jsut be 1

		{

		  fd << 0 << " ";

		}
	    }
	}
    }



}





/*
  Parya has added this method to ImageData class to write the maxima and sigma_theta information obtained with probabilistic approach together with labeling information to some given files.
*/

void ImageData::WriteODFMaxima(ofstream& fd, bool clobber, char *history)
{
  int x,y,z;
  int i;

  for (x=0;x<this->GetXSize();x++)
    {
      for (y=0;y<this->GetYSize();y++)
	{
	  for (z=0;z<this->GetZSize();z++)
	    {
	      if (m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated)
		{
		  fd << "(x,y,z)=(" << x << "," << y << "," << z << ")\n";
		  fd << this->GetODFMaxima(x,y,z)->GetNumberOfDirections() << "\n" ;

		  for (i=0;i<this->GetODFMaxima(x,y,z)->GetNumberOfDirections();i++)
		    {
		      fd << this->GetODFMaxima(x,y,z)->GetXComponent(i) << " " ;
		      fd << this->GetODFMaxima(x,y,z)->GetYComponent(i)<< " ";
		      fd << this->GetODFMaxima(x,y,z)->GetZComponent(i)<< "\n";
		    }
		}
	      /*      else //write all 0: confused because used to write 10 0s, should jsut be 1
		{
		  fd << 0 << "\n";
		}
	      */
	    }
	}
    }
}

void ImageData::WriteSigmaTheta(const char *filename, bool clobber, char *history, bool dble)
{
	int x,y,z;
	int i;
	ofstream fd(filename);

	if (!fd)
    {
		Error("cannot create file");
    }

	for (x=0;x<this->GetXSize();x++)
	{
		for (y=0;y<this->GetYSize();y++)
		{
			for (z=0;z<this->GetZSize();z++)
			{
				if (m_sigma_theta_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated)
				{
					fd << this->GetSigmaTheta(x,y,z).number << " ";

					if (dble == false)
					{
						for (i=0;i<this->GetSigmaTheta(x,y,z).number;i++)
						{
							fd << this->GetSigmaTheta(x,y,z).sigma_theta[i] << " " ;
						}
					}
					else if (dble == true)
					{
						for (i=0;i<2*this->GetSigmaTheta(x,y,z).number;i++)
						{
							fd << this->GetSigmaTheta(x,y,z).sigma_theta[i] << " " ;
						}
					}
				}
				else //write 0:
				{
					fd << 0 << " ";
				}
			}
		}
    }
}



void ImageData::WriteSigmaTheta(ofstream& fd, bool clobber, char *history)
{

  int x,y,z;
  int i;

  for (x=0;x<this->GetXSize();x++)
    {
      for (y=0;y<this->GetYSize();y++)
	{
	  for (z=0;z<this->GetZSize();z++)
	    {
	      if (m_sigma_theta_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated)
		{
		  fd << this->GetSigmaTheta(x,y,z).number << " " ;
		  for (i=0;i<this->GetSigmaTheta(x,y,z).number;i++)
		    {
		      fd << this->GetSigmaTheta(x,y,z).sigma_theta[i] << " " ;
		    }
		  fd << "\n" ;
		}
	      /*	      else //write 0:
		{
		  fd << 0 << "\n";
		}
	      */
	    }
	}
    }
}

/*
  Here ends the method added by Parya.
*/


void ImageData::SetODFMaxima(const char *filename)
{
  int x,y,z;
  int i;

//error if ImageData isn't ODF...

  Directions *directions;
  ImageData *maxima_image_data;

  int num_max;
  float max_x,max_y, max_z;


  m_maxima_array=(ODF_MAXIMA *) malloc(this->GetXSize()*this->GetYSize()*this->GetZSize()*sizeof(ODF_MAXIMA));
  m_maxima_array_exists=true;


  std::string fname(filename);
  //check if a minc file: if so, assume it is a Maxima3 (should check) 


  if (fname.rfind(".mnc")==(fname.length()-4))//MINC file
  {

      maxima_image_data=new ImageData(filename);




      //assign value directly to m_maxima_array and set appropriate flags:


      for (x=0;x<this->GetXSize();x++)
	{
	  for (y=0;y<this->GetYSize();y++)
	    {

	      for (z=0;z<this->GetZSize();z++)

		{
		  //here, we consider 0 dirs to mean not calculated (yet)
		  if (maxima_image_data->GetValue(x,y,z,0)>FLT_EPSILON)
		    {
		      //because the file is limited, clamp it.  hardcoded at 3 dirs:
		      directions=new Directions((int) float_min(maxima_image_data->GetValue(x,y,z,0),3.0));

		      for (i=0;i<directions->GetNumberOfDirections()&&i<3;i++)
			{
			  directions->SetXComponent(i,maxima_image_data->GetValue(x,y,z,3*i+1));
			  directions->SetYComponent(i,maxima_image_data->GetValue(x,y,z,3*i+2));
			  directions->SetZComponent(i,maxima_image_data->GetValue(x,y,z,3*i+3));
			}

		      m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated=true;

		      m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima=directions;
		    }
		  //else as if never calculated.


		}
	    }
	}



      //delete the new ImageData (for now; may want to hang on to it?? doubtful will want to write it again, and if we do, maxima calculation won't have to be redone, but could save)

      maxima_image_data->Release();
  } else {//maxima format file
      ifstream fd(filename);
      if (!fd)
	{
	  Error("cannot open file");
	}
      for (x=0;x<this->GetXSize();x++)
	{
	  for (y=0;y<this->GetYSize();y++)
	    {
	      for (z=0;z<this->GetZSize();z++)
		{

		  //read number of maxima:

		  fd >> num_max;
		  if (num_max>FLT_EPSILON)
		    {
		      directions=new Directions(num_max);


		      for (i=0;i<num_max;i++)
			{
			  fd >> max_x >> max_y >> max_z;
			  directions->SetXComponent(i,max_x);
			  directions->SetYComponent(i,max_y);
			  directions->SetZComponent(i,max_z);
			}

		      m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated=true;

		      m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima=directions;
		    }


		}
	    }
	}
    }

}


void ImageData::SetODFMaxima(Directions *directions,int x,int y,int z)
{
  int i;
  int xx,yy,zz;

  if (!m_maxima_array_exists)
    {
      m_maxima_array=(ODF_MAXIMA *) malloc(this->GetXSize()*this->GetYSize()*this->GetZSize()*sizeof(ODF_MAXIMA));
      m_maxima_array_exists=true;


      //loop through and set .calculated to false:
      for (xx=0;xx<this->GetXSize();xx++)
	{
	  for (yy=0;yy<this->GetYSize();yy++)
	    {
	      for (zz=0;zz<this->GetZSize();zz++)
		{

		  m_maxima_array[zz*this->GetYSize()*this->GetXSize()+yy*this->GetXSize()+xx].calculated=false;
		}
	    }
	}
    }

  m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated=true;
  m_maxima_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].maxima=directions;

  directions->Retain();

}


void ImageData::SetSigmaTheta(SIGMA_THETA sigma_theta,int x,int y,int z)
{

  int xx,yy,zz;
  if (!m_sigma_theta_array_exists)
    {
      m_sigma_theta_array=(SIGMA_THETA *) malloc(this->GetXSize()*this->GetYSize()*this->GetZSize()*sizeof(SIGMA_THETA));
      //loop through and set .calculated to false:
      for (xx=0;xx<this->GetXSize();xx++)
	{
	  for (yy=0;yy<this->GetYSize();yy++)
	    {
	      for (zz=0;zz<this->GetZSize();zz++)
		{

		  m_sigma_theta_array[zz*this->GetYSize()*this->GetXSize()+yy*this->GetXSize()+xx].calculated=false;
		}
	    }
	}

      m_sigma_theta_array_exists=true;
    }

  m_sigma_theta_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x]=sigma_theta;


}

bool ImageData::SigmaThetaAvailable()
{
  return m_sigma_theta_array_exists;
}

SIGMA_THETA ImageData::GetSigmaTheta(int x,int y,int z)
{
  //if (m_sigma_theta_array_exists) && .calculated...
  return m_sigma_theta_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x];
}

void ImageData::WriteSigmaTheta(const char *filename, bool clobber, char *history)
{

  int x,y,z;
  int i;
  ofstream fd(filename);

  if (!fd)
    {
      Error("cannot create file");
    }

  for (x=0;x<this->GetXSize();x++)
    {
      for (y=0;y<this->GetYSize();y++)
	{

	  for (z=0;z<this->GetZSize();z++)

	    {
	      if (m_sigma_theta_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated)
		{
		  //test:
		  //cout << "write sigma theta: " << this->GetSigmaTheta(x,y,z).number << endl;
		  fd << this->GetSigmaTheta(x,y,z).number << " " ;
		  for (i=0;i<this->GetSigmaTheta(x,y,z).number;i++)
		    {
		      fd << this->GetSigmaTheta(x,y,z).sigma_theta[i] << " " ;
		    }
		}
	      else //write 0:
		{
		  fd << 0 << " ";

		}
	    }
	}
    }



}

void ImageData::SetSigmaTheta(const char *filename)
{
  int x,y,z;
  int i;

  ifstream fd(filename);
  int num_max;
  SIGMA_THETA sigma_theta;
  float st;

  m_sigma_theta_array=(SIGMA_THETA *) malloc(this->GetXSize()*this->GetYSize()*this->GetZSize()*sizeof(SIGMA_THETA));
  m_sigma_theta_array_exists=true;

  if (!fd)
    {
      Error("cannot open file");
    }
      for (x=0;x<this->GetXSize();x++)
	{
	  for (y=0;y<this->GetYSize();y++)
	    {
	      for (z=0;z<this->GetZSize();z++)
		{
		 fd >> num_max;
		 if (num_max>FLT_EPSILON)
		   {
		     sigma_theta.sigma_theta=(float *) malloc(num_max*sizeof(float));
		     for (i=0;i<num_max;i++)
		       {
			 fd >> st;
			 sigma_theta.sigma_theta[i]=st;
		       }
		     m_sigma_theta_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x]=sigma_theta;
		     m_sigma_theta_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated=true;
		   }
		 else
		   {
		     m_sigma_theta_array[z*this->GetYSize()*this->GetXSize()+y*this->GetXSize()+x].calculated=false;
		   }
		}
	    }
	}
}



IMAGE_GEOM ImageData::GetImageGeom()
{

  if (!m_image_geom_init)
    {
    
  //IMAGE_GEOM image_geom;

  

  m_image_geom.x_size=this->GetXSize();
  m_image_geom.y_size=this->GetYSize();
  m_image_geom.z_size=this->GetZSize();
  m_image_geom.x_step=this->GetXStep();
  m_image_geom.y_step=this->GetYStep();
  m_image_geom.z_step=this->GetZStep();
  m_image_geom.x_start=this->GetXStart();
  m_image_geom.y_start=this->GetYStart();
  m_image_geom.z_start=this->GetZStart();



  if (this->GetDimensions()>3)
    {
      m_image_geom.n_size=this->GetNSize();
      m_image_geom.n_start=this->GetNStart();
      m_image_geom.n_step=this->GetNStep();
    }
  m_image_geom.x_direction_cosine=this->GetXDirectionCosine();
  m_image_geom.y_direction_cosine=this->GetYDirectionCosine();
  m_image_geom.z_direction_cosine=this->GetZDirectionCosine();



  m_image_geom.point_inverse_xfm=NULL;
  m_image_geom.vector_inverse_xfm=NULL;
    
  m_image_geom_init=true;
    }
  return m_image_geom;
    

}


void ImageData::CreateNeighbourLUT()
{

  //test:
  int counter;

  int i,j;
  float radius;
  if (!(_image_type==DWIdata || _image_type==ODFdata))
    {
      Error("ImageData is not an ODF");
    }


  //just in case:
  this->GetDirections()->Normalize();



  radius=2.3*sqrt(2*M_PI/this->GetNSize()); //1.5-2.3 is visually nice. This has only been tested for N~100


  m_neighbour_LUT=(NeighbourLUTEntry **) malloc(this->GetNSize()*sizeof(NeighbourLUTEntry *));


  for (i=0;i<this->GetNSize();i++)
    {
      m_neighbour_LUT[i]=NULL;
    }

  for (i=0;i<this->GetNSize();i++)
    {
      counter=0;
      if ((pow(this->GetDirections()->GetXComponent(i),2)+pow(this->GetDirections()->GetYComponent(i),2)+pow(this->GetDirections()->GetZComponent(i),2))>FLT_EPSILON)
	{


	  for (j=0;j<this->GetNSize();j++)
	    {
	      //don't use "0" dirs:

	      if ((pow(this->GetDirections()->GetXComponent(j),2)+pow(this->GetDirections()->GetYComponent(j),2)+pow(this->GetDirections()->GetZComponent(j),2))>FLT_EPSILON && j!=i)
		{


		  if (float_abs((this->GetDirections()->GetXComponent(i)*this->GetDirections()->GetXComponent(j)+this->GetDirections()->GetYComponent(i)*this->GetDirections()->GetYComponent(j)+this->GetDirections()->GetZComponent(i)*this->GetDirections()->GetZComponent(j)))>float_abs(cos(radius)))
		  //test:
		  //if (float_abs((this->GetDirections()->GetXComponent(i)*this->GetDirections()->GetXComponent(j)+this->GetDirections()->GetYComponent(i)*this->GetDirections()->GetYComponent(j)+this->GetDirections()->GetZComponent(i)*this->GetDirections()->GetZComponent(j)))>0.85)
		    {
		      this->AddNeighbour(i,j);
		      counter++;
		      //cout << "dot prod " << this->GetDirections()->GetXComponent(i)*this->GetDirections()->GetXComponent(j)+this->GetDirections()->GetYComponent(i)*this->GetDirections()->GetYComponent(j)+this->GetDirections()->GetZComponent(i)*this->GetDirections()->GetZComponent(j) << endl;
		    }
		}

	    }
	  if (counter<1)
	    {
	      Error("radius too small");//DO: if ever get a 0, should rerun this with a larger radius.
	    }
	  //cout << counter << endl;
	}
    }
  m_neighbour_LUT_calculated=true;

}


void ImageData::AddNeighbour(int entry, int neighbour)
{



  //link to head, not end:
  NeighbourLUTEntry *new_entry=( NeighbourLUTEntry *) malloc(sizeof(NeighbourLUTEntry));
  new_entry->index=neighbour;
  new_entry->next=m_neighbour_LUT[entry];
  m_neighbour_LUT[entry]=new_entry;


  /*
    NeighbourLUTEntry *this_entry=&m_neighbour_LUT[entry];
    while ((*this_entry).next!=NULL)
    {
      this_entry=(*this_entry).next;
    }
  NeighbourLUTEntry *new_entry=( NeighbourLUTEntry *) malloc(sizeof(NeighbourLUTEntry));
  (*this_entry).next=new_entry;
  (*new_entry).index=neighbour;
  (*new_entry).next=NULL;
  */
}

void ImageData::NormalizeODF(int voxel_x, int voxel_y, int voxel_z)
{

  //DO: require this be ODF data!! Currently, only volume normalization is permitted

  float norm_factor;
  float scale;
  float value;


  int n;

  //normalize so that all ODFs from this acquisition have unit volume:

  //need to bring radius down to (eg, 1) first, or this blows up:

  norm_factor=0;
  scale =0;

  //get the value of the first direction and set a scale factor (only done if the first direction had a nonzero value - alt is to find the max value and use it...not implemented):
  scale=this->GetValue(voxel_x,voxel_y,voxel_z,0);
  if (scale>FLT_EPSILON)
    {
      for(n=0;n<this->GetDirections()->GetNumberOfDirections();n++)
	{
	  this->SetValue(voxel_x,voxel_y,voxel_z,n,this->GetValue(voxel_x,voxel_y,voxel_z,n)/scale);
	}
    }

  for(n=0;n<this->GetDirections()->GetNumberOfDirections();n++)
    {

      norm_factor+=(4*M_PI*pow(this->GetValue(voxel_x,voxel_y,voxel_z,n),3))/(3*this->GetDirections()->GetNumberOfDirections());
    }


      for(n=0;n<this->GetDirections()->GetNumberOfDirections();n++)
	{
	  this->SetValue(voxel_x,voxel_y,voxel_z,n,this->GetValue(voxel_x,voxel_y,voxel_z,n)/norm_factor);
	}

}

//make ODF out of nonzero directions/bvalues: currently has to do a deep copy anyway if no b=0; for speed, could return a flag if there are no b=0 and put the ouput as an argument.  need a Createb0Series() here too.
ImageData *ImageData::CreateDWIODF()
{
  int i;
  int counter=0;
  int x,y,z;

  //find out how many b=0 to chop off: look for direction 0 0 0 but could use b value instead
  for (i=0;i<this->GetNSize();i++)
    {
      if ((this->GetDirections()->GetXComponent(i)*this->GetDirections()->GetXComponent(i)+
           this->GetDirections()->GetYComponent(i)*this->GetDirections()->GetYComponent(i)+
           this->GetDirections()->GetZComponent(i)*this->GetDirections()->GetZComponent(i))>FLT_EPSILON)
	{
	  counter++;
	}
    }

  ImageData *image_data_no_zero=new ImageData(DWIdata,
                                              this->GetXSize(),this->GetYSize(),this->GetZSize(),counter,
                                              this->GetXStart(),this->GetYStart(),this->GetZStart(),0,
                                              this->GetXStep(),this->GetYStep(),this->GetZStep(),1,
                                              this->GetXDirectionCosine(),this->GetYDirectionCosine(),this->GetZDirectionCosine());

  Directions *directions=new Directions(counter);

  //go through and assign:
  counter=0;
  
  for (i=0;i<this->GetNSize();i++)
    {
      if ((this->GetDirections()->GetXComponent(i)*this->GetDirections()->GetXComponent(i)+
           this->GetDirections()->GetYComponent(i)*this->GetDirections()->GetYComponent(i)+
           this->GetDirections()->GetZComponent(i)*this->GetDirections()->GetZComponent(i))>FLT_EPSILON)
      
	{
 	  for (x=0;x<this->GetXSize();x++)
	    {
	      for (y=0;y<this->GetYSize();y++)
		{
		  for (z=0;z<this->GetZSize();z++)
		    {

		      image_data_no_zero->SetValue(x,y,z,counter,this->GetValue(x,y,z,i));


		    }
		}
	    }
	  directions->SetDirection(counter,this->GetDirections()->GetDirection(i));
	  image_data_no_zero->GetDiffusionParameters()->SetBValue(counter,this->GetDiffusionParameters()->GetBValue(i));
	  counter++;
	}
    }


  image_data_no_zero->SetDirections(directions);

  directions->Release();


  return image_data_no_zero;

}

int ImageData::FindClosestDirection(VECTOR3D TanVec)
{
	int i, ind;
	VECTOR3D dir;

	float max_dot = FLT_EPSILON;
	float dotprod;

	if (!_directions_are_set)
	{
		Error("No Directions set to compare to!!!");
	}

	for (i=0; i<_directions->GetNumberOfDirections(); i++)
	{
		dir = _directions->GetDirection(i);
		dotprod = float_abs(dot(TanVec, dir));
		dotprod = dotprod/(sqrt(dot(dir, dir))*sqrt(dot(TanVec, TanVec)));

		if (dotprod > max_dot)
		{
			max_dot = dotprod;
			ind = i;
		}
	}

	return ind;
}

