/************************************************************************

   File: misc.c++

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>

#include "misc.h"
#include "otherlibs_include.h"
#include "Directions.h"
#include "ImageData.h"
#include "S2Interpolator.h"

#include <boost/random.hpp>

//using namespace std;


void Error(const char *msg) 

{
  //  fprintf(stderr,"%s\n",msg);
  cerr << msg << endl;
  exit(1);
}

char *CreateHistoryString(int argc, char *argv[])
{
  int i;
  char *history;

  int l = 2;
  for(i = 0; i < argc; i++) 
    {
     l += strlen(argv[i]) + 1;
    }
  
  history = (char *) malloc(l*sizeof(char));

  strstream historystream(history,l*sizeof(char));
  
  historystream << "";
  
  for(i = 0; i < argc; i++) 
    {
      historystream << argv[i] << " ";
    }
  historystream << ends;
  
  
  return history;
  
}

void ResampleStandardDirCos(char *in_filename,char *out_filename,char
*tempdirname, bool nearest)
{
  char command[1000];
  char *tempfilename=(char *) malloc(200*sizeof(char));


  char *tempfilename2;
  char *tempxfmname=(char *) malloc(200*sizeof(char)); // IL added for xfm file


  int seed;

  float xstart, ystart, zstart, xstep, ystep, zstep;

  int xnelements, ynelements, znelements;

  VECTOR3D dircosx, dircosy, dircosz;


  seed = (int) time (NULL);
  srand48(seed);

  sprintf(tempfilename,"%s/tmp%f.mnc",tempdirname, (double) rand());
  sprintf(tempxfmname,"%s/tmpxfm%f.xfm",tempdirname, (double) rand()); // IL



  // check how the current direction cosines relate to std and resample if neg dir cosine:

  dircosx.x=1;
  dircosx.y=0;
  dircosx.z=0;

  dircosy.x=0;
  dircosy.y=1;
  dircosy.z=0;

  dircosz.x=0;
  dircosz.y=0;
  dircosz.z=1;

  ImageData *image= new ImageData(in_filename);




  tempfilename2=in_filename;
  // IL create indentity matrix with standard dircos directions
  sprintf(command,"param2xfm %s\n",tempxfmname);
  system(command);



  if (nearest)
    {
      // IL resample to new file using default direction cosine setup
      sprintf(command,"mincresample -nearest -tfm_input_sampling -xdircos 1 0 0 -ydircos 0 1 0 -zdircos 0 0 1 -transformation %s %s %s\n",tempxfmname,in_filename, tempfilename); //IL

    }
  else
    {
      // IL resample to new file using default direction cosine setup
      sprintf(command,"mincresample -tfm_input_sampling -xdircos 1 0 0 -ydircos 0 1 0 -zdircos 0 0 1 -transformation %s %s %s\n",tempxfmname,in_filename, tempfilename); //IL

    }

  system(command);

  strcpy(out_filename,tempfilename);

  image->Release();

  free(tempfilename);
  free(tempxfmname);


}



void ResampleNegStep(char *in_filename,char *out_filename,char *tempdirname)
{
  char command[1000];
  char *tempfilename=(char *) malloc(200*sizeof(char));
  int seed;
  float xstart,ystart,zstart,zstep,ystep,xstep;
  seed = (int) time (NULL);
  srand48(seed); 
  
  sprintf(tempfilename,"%s/tmp%f.mnc",tempdirname, (double) rand());
  //sprintf(tempfilename,"%s/tmp%f.mnc",tempdirname, (double) rand());
  
  //get info on infile:

  ImageData *image_data=new ImageData(in_filename);

  if (image_data->GetXStep()<0)
    {
      xstart=image_data->GetXStart();
      xstep=image_data->GetXStep();
      
    }
  else
    {
      xstart=image_data->GetXStart()+(image_data->GetXSize()-1)*image_data->GetXStep();
      xstep=-1.0*image_data->GetXStep();
    }

  if (image_data->GetYStep()>0)
    {
      ystart=image_data->GetYStart()+(image_data->GetYSize()-1)*image_data->GetYStep();
      ystep=-1.0*image_data->GetYStep();
    }
  else
    {
      ystart=image_data->GetYStart();
      ystep=image_data->GetYStep();

    }
  if (image_data->GetZStep()>0)
    {
      zstart=image_data->GetZStart()+(image_data->GetZSize()-1)*image_data->GetZStep();
      zstep=-1.0*image_data->GetZStep();
    }
  else
    {
      zstart=image_data->GetZStart();
      zstep=image_data->GetZStep();
    }


  sprintf(command,"mincresample -nearest -xstep %f -ystep %f -zstep %f -xstart %f -ystart %f -zstart %f %s %s\n",xstep,ystep,zstep,xstart,ystart,zstart,in_filename,tempfilename); 
  

  system(command);

  //sprintf(command,"mincreshape +direction %s %s\n",tempfilename,
  
  strcpy(out_filename,tempfilename);
  image_data->Release();

  free(tempfilename);

}

void ResamplePosStep(char *in_filename,char *out_filename,char *tempdirname)
{
  char command[1000];
  char *tempfilename=(char *) malloc(200*sizeof(char));
  int seed;
  float xstart,ystart,zstart,zstep,ystep,xstep;
  seed = (int) time (NULL);
  srand48(seed); 
  
  sprintf(tempfilename,"%s/tmp%f.mnc",tempdirname, (double) rand());
  //sprintf(tempfilename,"%s/tmp%f.mnc",tempdirname, (double) rand());
  
  //get info on infile:

  ImageData *image_data=new ImageData(in_filename);

  if (image_data->GetXStep()>0)
    {
      xstart=image_data->GetXStart();
      xstep=image_data->GetXStep();
      
    }
  else
    {
      xstart=image_data->GetXStart()+(image_data->GetXSize()-1)*image_data->GetXStep();
      xstep=-1.0*image_data->GetXStep();
    }
  if (image_data->GetYStep()<0)
    {
      ystart=image_data->GetYStart()+(image_data->GetYSize()-1)*image_data->GetYStep();
      ystep=-1.0*image_data->GetYStep();
    }
  else
    {
      ystart=image_data->GetYStart();
      ystep=image_data->GetYStep();

    }
  if (image_data->GetZStep()<0)
    {
      zstart=image_data->GetZStart()+(image_data->GetZSize()-1)*image_data->GetZStep();
      zstep=-1.0*image_data->GetZStep();
    }
  else
    {
      zstart=image_data->GetZStart();
      zstep=image_data->GetZStep();
    }


  sprintf(command,"mincresample -xstep %f -ystep %f -zstep %f -xstart %f -ystart %f -zstart %f %s %s\n",xstep,ystep,zstep,xstart,ystart,zstart,in_filename,tempfilename); 
  

  system(command);

  //sprintf(command,"mincreshape +direction %s %s\n",tempfilename,
  
  strcpy(out_filename,tempfilename);
  image_data->Release();


}


void Resample(char* in_filename,char* out_filename,char* tempdirname,float step,bool use_input_sampling,char*templatefilename, bool& templatefilename_set,bool nearest)
{
  std::string _in_filename=in_filename;
  std::string _out_filename=out_filename;
  std::string _tempdirname=tempdirname;
  std::string _templatefilename=templatefilename;
  Resample(_in_filename,_out_filename,_tempdirname,step,use_input_sampling,_templatefilename,templatefilename_set,nearest);
  
  strcpy(in_filename,_in_filename.c_str());
  strcpy(out_filename,_out_filename.c_str());
  strcpy(tempdirname,_tempdirname.c_str());
  strcpy(templatefilename,_templatefilename.c_str());
  
}

void Resample(std::string& in_filename,std::string& out_filename,std::string& tempdirname,float step,bool use_input_sampling,std::string& templatefilename, bool& templatefilename_set, bool nearest)
{
  //nomatter what you do here, the file should be resampled to be zyx,because it appears that volume_io isn't handling xyz?  no, it should be handling it: but -transverse is still here: can remove...


  char command[1000];
  char *tempfilename=(char *) malloc(200*sizeof(char));
  ImageData *image_data=NULL;
  int seed;
  int xnelements,ynelements,znelements;
  float xstart, ystart, zstart;

  //if (templatefilename==NULL)
  if (templatefilename_set==false)
    {
      seed = (int) time (NULL);
      srand48(seed); 
    }

  //else 
  //{
  // cout << "Resample: templatefilename " << templatefilename << endl; 
  //}
  
  //this appears to be unique, but check it: (could just use a counter..) this gives the same series of numbers every time program is run, regardless of seed(=time)
  sprintf(tempfilename,"%s/tmp%f.mnc",tempdirname.c_str(), (double) rand());

  
  
  //cout << "tempfilename: " <<  tempfilename << endl;
 
      //if (templatefilename==NULL)
  if (templatefilename_set==false)
    {
      //Q: do I now need to allocate enough space for templatefilename?
      //templatefilename=(char *) malloc(200*sizeof(char));
      templatefilename=tempfilename;


      templatefilename_set=true;
      //cout << "now set: templatefilename: " << templatefilename << endl;
      
      if (use_input_sampling)
	{

    sprintf(command,"cp %s %s\n",in_filename.c_str(),tempfilename);
	  //sprintf(command,"mincresample -transverse %s %s\n",in_filename,tempfilename);
	  //cout << command;
	  system(command);
	}
      else 
	{
	  //mincresample it if abs value of step differs by >0.1; get step from mincinfo: have to write to file and read back in, or create an ImageData and delete it.  

    image_data=new ImageData(in_filename.c_str());

	  if (float_abs(float_abs(step)-float_abs(image_data->GetXStep()))>0.1 || float_abs(float_abs(step)-float_abs(image_data->GetYStep()))>0.1 || float_abs(float_abs(step)-float_abs(image_data->GetZStep()))>0.1)
	    {
	      xnelements=(int) float_abs((image_data->GetXSize()*(image_data->GetXStep()/step)));
	      ynelements=(int) float_abs((image_data->GetYSize()*(image_data->GetYStep()/step)));
	      znelements=(int) float_abs((image_data->GetZSize()*(image_data->GetZStep()/step)));
	      if (step<0) 
		{
		  Error("step must be positive");//IRL Error for gcc4
		}
	      if (image_data->GetXStep()<0)
		{
		  xstart=image_data->GetXStart()+(image_data->GetXSize()-1)*image_data->GetXStep();
		}
	      else
		{
		  xstart=image_data->GetXStart();
		}
	      if (image_data->GetYStep()<0)
		{
		  ystart=image_data->GetYStart()+(image_data->GetYSize()-1)*image_data->GetYStep();
		  }
	      else
		{
		  ystart=image_data->GetYStart();
		}
	      if (image_data->GetZStep()<0)
		{
		  zstart=image_data->GetZStart()+(image_data->GetZSize()-1)*image_data->GetZStep();
		}
	      else
		{
		  zstart=image_data->GetZStart();
		}

	      //sprintf(command,"mincresample -transverse  -xstep %f -ystep %f -zstep %f -xstart %f -ystart %f -zstart %f -xnelements %i -ynelements %i -znelements %i %s %s\n",step,step,step,xstart, ystart, zstart, xnelements,ynelements,znelements,in_filename,tempfilename);
	      if (nearest)
		{
      sprintf(command,"mincresample -nearest -xstep %f -ystep %f -zstep %f -xstart %f -ystart %f -zstart %f -xnelements %i -ynelements %i -znelements %i %s %s\n",step,step,step,xstart, ystart, zstart, xnelements,ynelements,znelements,in_filename.c_str(),tempfilename);
		}
	      else
		{
      sprintf(command,"mincresample -xstep %f -ystep %f -zstep %f -xstart %f -ystart %f -zstart %f -xnelements %i -ynelements %i -znelements %i %s %s\n",step,step,step,xstart, ystart, zstart, xnelements,ynelements,znelements,in_filename.c_str(),tempfilename);
		}
	      //cout << command;
	      system(command);
	    }
	  else
	    {
        sprintf(command,"cp %s %s\n",in_filename.c_str(),tempfilename);
	      //sprintf(command,"mincresample -transverse  %s %s\n",in_filename,tempfilename);

	      system(command);

	    }
	}
    }
  else 
    {
      //resample -like templatefilename IF the sampling is different (at a certain precision or not transverse, here I do it in all cases):
      //DO if checking: transverse and then create a new ImageData...
      if (true)
	{
	  if (nearest)
	    {
        sprintf(command,"mincresample -nearest -like %s %s %s\n", 
                templatefilename.c_str(),in_filename.c_str(), tempfilename);
	    }
	  else
	    {
        sprintf(command,"mincresample -like %s %s %s\n", 
                templatefilename.c_str(),in_filename.c_str(), tempfilename);
	    }
	  //cout << command;
	  system(command);
	}
      else
	{
	  //just copy it:
    sprintf(command,"cp %s %s\n",
            in_filename.c_str(),tempfilename);
	  //cout << command; 
	  system(command);
	}
    }

  if (image_data!=NULL)
    {
      image_data->Release();
    }

  out_filename=tempfilename;

  free(tempfilename);

}


float float_min3(float a,float b,float c) 
{

  if (a<b && a<c) return a;
  else if (b<a && b<c) return b;
  else return c;

}


float float_max3(float a,float b,float c) 
{

  if (a>b && a>c) return a;
  else if (b>a && b>c) return b;
  else return c;

}

int int_min3(int a,int b,int c) 
{

  if (a<b && a<c) return a;
  else if (b<a && b<c) return b;
  else return c;

}


int int_max3(int a,int b,int c) 
{

  if (a>b && a>c) return a;
  else if (b>a && b>c) return b;
  else return c;

}

float float_max(float a, float b)
{
  if (a>b)
    {
      return a;
      
    }
  else 
    {
      return b;
      
    }
  
 
}

float float_min(float a, float b)
{
  if (a<b)
    {
      return a;
      
    }
  else 
    {
      return b;
      
    }
  
 
}


float float_abs(float a)
{
  if (a<0)
    {
      return -a;
      
    }
 else
   {
     return a;
     
   }
  
}

int int_abs(int a)
{
  if (a<0)
    {
      return -a;
      
    }
 else
   {
     return a;
     
   }
}

float flt_norm(VECTOR3D vector)
{
  return sqrt(vector.x*vector.x+vector.y*vector.y+vector.z*vector.z);
}

VECTOR3D Normalize(VECTOR3D vector)
{
  VECTOR3D out_vect;
  float len;

  len=sqrt(vector.x*vector.x+vector.y*vector.y+vector.z*vector.z);

  out_vect.x=vector.x/len;
  out_vect.y=vector.y/len;
  out_vect.z=vector.z/len;

  return out_vect;
}

void MakeTwoOrthogonalVectors(float vector[3], float *basis1, float *basis2)
{
  float len;
  

  //all planes will contain an orth. vector.  we will choose one:

 

  if (!(vector[1]==0 && vector[2]==0)) { //equivalently, vector[0]!=1
    basis1[0]=1;
    if (!(vector[0]==0 && vector[2]==0)) { //equivalently, vector[1]!=1
      if (vector[2]!=0)
	{
	  
	  basis1[1]=1;
	  basis1[2]=-1*(vector[0]+vector[1])/vector[2];
	}
      else
	{
	  basis1[2]=1;
	  basis1[1]=-1*vector[1];
	}
    }
    else { //vector[2]!=1 because (0,1,0)
      basis1[2]=1;
      basis1[1]=0;      
    }
  }
  
  else { 
    basis1[0]=0;
    basis1[1]=1;
    basis1[2]=1;
  }


  //normalize it:

  len=sqrt(basis1[0]*basis1[0]+basis1[1]*basis1[1]+basis1[2]*basis1[2]);
  basis1[0]=basis1[0]/len;
  basis1[1]=basis1[1]/len;
  basis1[2]=basis1[2]/len;
  
  //take the cross product:

  basis2[0]=vector[1]*basis1[2]-vector[2]*basis1[1];
  basis2[1]=-1*(vector[0]*basis1[2]-vector[2]*basis1[0]);
  basis2[2]=vector[0]*basis1[1]-vector[1]*basis1[0];

  //normalize it:

  len=sqrt(basis2[0]*basis2[0]+basis2[1]*basis2[1]+basis2[2]*basis2[2]);
  basis2[0]=basis2[0]/len;
  basis2[1]=basis2[1]/len;
  basis2[2]=basis2[2]/len;
  
	  

  return;
  


}


void MakeTwoOrthogonalVectors(VECTOR3D vector, VECTOR3D &basis1, VECTOR3D &basis2)
{
  float len;
  

  //all planes will contain an orth. vector.  we will choose one:

 

  if (!(vector.y==0 && vector.z==0)) { //equivalently, vector.x!=1
    basis1.x=1;
    if (!(vector.x==0 && vector.z==0)) { //equivalently, vector.y!=1
      if (vector.z!=0)
	{
	  basis1.y=1;
	  basis1.z=-1*(vector.x+vector.y)/vector.z;
	}
      else
	{
	  basis1.z=1;
	  basis1.y=-1*vector.y;
	}
    }
    else { //vector.z!=1 because (0,1,0)
      basis1.z=1;
      basis1.y=0;      
    }
  }
  
  else { 
    basis1.x=0;
    basis1.y=1;
    basis1.z=1;
  }


  //normalize it:

  len=sqrt(basis1.x*(basis1.x)+basis1.y*(basis1.y)+basis1.z*(basis1.z));
  basis1.x=basis1.x/len;
  basis1.y=basis1.y/len;
  basis1.z=basis1.z/len;
  
  //take the cross product:

  basis2.x=vector.y*(basis1.z)-vector.z*(basis1.y);
  basis2.y=-1*(vector.x*(basis1.z)-vector.z*(basis1.x));
  basis2.z=vector.x*(basis1.y)-vector.y*(basis1.x);

  //normalize it:

  len=sqrt(basis2.x*(basis2.x)+basis2.y*(basis2.y)+basis2.z*(basis2.z));
  basis2.x=basis2.x/len;
  basis2.y=basis2.y/len;
  basis2.z=basis2.z/len;
  
	  

  return;
  


}



int int_min(int a, int b)
{
  if (a<b)
    {
      return a;
      
    }
  else 
    {
      return b;
      
    }
  
}

int int_max(int a, int b)
{
  if (a>b)
    {
      return a;
      
    }
  else 
    {
      return b;
      
    }
  
}

VECTOR3D cross(VECTOR3D v1, VECTOR3D v2)
{
  VECTOR3D result;
  result.x=v1.y*v2.z-v1.z*v2.y;
  result.y=-v1.x*v2.z+v1.z*v2.x;
  result.z=v1.x*v2.y-v1.y*v2.x;

  return result;

}

float dot(VECTOR3D v1, VECTOR3D v2)
{
  return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;


}

VECTOR3D neg(VECTOR3D v)
{
  VECTOR3D v2;
  v2.x=-1.0*v.x;
  v2.y=-1.0*v.y;
  v2.z=-1.0*v.z;
}


POINT3D VoxelToWorld(float x,float y,float z, IMAGE_GEOM image_geom)
{


  POINT3D point;

  
  point.x=(image_geom.x_start+image_geom.x_step*x)*image_geom.x_direction_cosine.x+(image_geom.y_start+image_geom.y_step*y)*image_geom.y_direction_cosine.x+(image_geom.z_start+image_geom.z_step*z)*image_geom.z_direction_cosine.x;

  point.y=(image_geom.x_start+image_geom.x_step*x)*image_geom.x_direction_cosine.y+(image_geom.y_start+image_geom.y_step*y)*image_geom.y_direction_cosine.y+(image_geom.z_start+image_geom.z_step*z)*image_geom.z_direction_cosine.y;

  point.z=(image_geom.x_start+image_geom.x_step*x)*image_geom.x_direction_cosine.z+(image_geom.y_start+image_geom.y_step*y)*image_geom.y_direction_cosine.z+(image_geom.z_start+image_geom.z_step*z)*image_geom.z_direction_cosine.z;

  


  return point;

}


POINT3D VoxelToWorld(POINT3D inpoint, IMAGE_GEOM image_geom)
{
  return VoxelToWorld(inpoint.x,inpoint.y,inpoint.z,image_geom);

}


VECTOR3D VoxelToWorld(VECTOR3D input_vector, IMAGE_GEOM image_geom)
{
  /*IRL just rotate, do not scale*/
  //POINT3D r=VoxelToWorld(input_vector.x,input_vector.y,input_vector.z,image_geom);
  //return VECTOR3D(r.x,r.y,r.z)
  //rotation takes into account the voxel steps, possibly negative or anisotropic


  VECTOR3D v,vt;

  //account for steps:

  vt.x=image_geom.x_step*input_vector.x;
  vt.y=image_geom.y_step*input_vector.y;
  vt.z=image_geom.z_step*input_vector.z;

  vt=Normalize(vt);

  //JC changed this to transpose of what it was before

  v.x=vt.x*image_geom.x_direction_cosine.x+
           vt.y*image_geom.y_direction_cosine.x+
           vt.z*image_geom.z_direction_cosine.x;
           
  v.y=vt.x*image_geom.x_direction_cosine.y+
           vt.y*image_geom.y_direction_cosine.y+
           vt.z*image_geom.z_direction_cosine.y;
           

  v.z=vt.x*image_geom.x_direction_cosine.z+
           vt.y*image_geom.y_direction_cosine.z+
           vt.z*image_geom.z_direction_cosine.z;





  return v;
}



POINT3D WorldToVoxel(float x,float y,float z, IMAGE_GEOM *image_geom)
{

  POINT3D point;

  INVERSE *inv;
 
  gsl_vector *inter;
  //inter=gsl_vector_alloc(3);
  //gsl_vector *gslpoint;
  //gslpoint=gsl_vector_alloc(3);
  gsl_matrix *xform;
  gsl_permutation *permutation;
  gsl_matrix *inverse;

  if (image_geom->point_inverse_xfm==NULL)
    {
      image_geom->point_inverse_xfm=new INVERSE();
     
    }

 
  inv=image_geom->point_inverse_xfm;
  
  if (inv->calculated==true)
    {
      inter=inv->inter;
    }
  else
    {
    
      inter=gsl_vector_alloc(3);
     
    }

  //inter was saved just so that we can skip alloc/free each call
  /*
  gsl_vector_set(new_inv->inter,0,(double) x-(image_geom->x_start*image_geom->x_direction_cosine.x+image_geom->y_start*image_geom->y_direction_cosine.x+image_geom->z_start*image_geom->z_direction_cosine.x));
  gsl_vector_set(new_inv->inter,1,(double) y-(image_geom->x_start*image_geom->x_direction_cosine.y+image_geom->y_start*image_geom->y_direction_cosine.y+image_geom->z_start*image_geom->z_direction_cosine.y));
  gsl_vector_set(new_inv->inter,2,(double) z-(image_geom->x_start*image_geom->x_direction_cosine.z+image_geom->y_start*image_geom->y_direction_cosine.z+image_geom->z_start*image_geom->z_direction_cosine.z));
  */

    gsl_vector_set(inter,0,(double) x-(image_geom->x_start*image_geom->x_direction_cosine.x+image_geom->y_start*image_geom->y_direction_cosine.x+image_geom->z_start*image_geom->z_direction_cosine.x));
  gsl_vector_set(inter,1,(double) y-(image_geom->x_start*image_geom->x_direction_cosine.y+image_geom->y_start*image_geom->y_direction_cosine.y+image_geom->z_start*image_geom->z_direction_cosine.y));
  gsl_vector_set(inter,2,(double) z-(image_geom->x_start*image_geom->x_direction_cosine.z+image_geom->y_start*image_geom->y_direction_cosine.z+image_geom->z_start*image_geom->z_direction_cosine.z));
  


  if (inv->calculated==false)
    {

      permutation=gsl_permutation_alloc(3);
      int p;

      //create xform 
  
      xform=gsl_matrix_alloc(3,3);

      gsl_matrix_set(xform,0,0,(double) image_geom->x_step*image_geom->x_direction_cosine.x);
      gsl_matrix_set(xform,0,1,(double) image_geom->y_step*image_geom->y_direction_cosine.x);
      gsl_matrix_set(xform,0,2,(double) image_geom->z_step*image_geom->z_direction_cosine.x);

      gsl_matrix_set(xform,1,0,(double) image_geom->x_step*image_geom->x_direction_cosine.y);
      gsl_matrix_set(xform,1,1,(double) image_geom->y_step*image_geom->y_direction_cosine.y);
      gsl_matrix_set(xform,1,2,(double) image_geom->z_step*image_geom->z_direction_cosine.y);

      gsl_matrix_set(xform,2,0,(double) image_geom->x_step*image_geom->x_direction_cosine.z);
      gsl_matrix_set(xform,2,1,(double) image_geom->y_step*image_geom->y_direction_cosine.z);
      gsl_matrix_set(xform,2,2,(double) image_geom->z_step*image_geom->z_direction_cosine.z);
  
      inverse=gsl_matrix_alloc(3,3);
  

      gsl_linalg_LU_decomp(xform,permutation,&p);
      gsl_linalg_LU_invert (xform,permutation,inverse);
      //gsl_matrix_free(xform);
      //gsl_permutation_free(permutation);
    
    }
  else
    {
      inverse=inv->xform;
      //permutation=image_geom->point_inverse_xfm.permutation;
    }

  //gsl_linalg_LU_solve(xform,permutation,inter,gslpoint);



  /*
  gsl_blas_dgemv(CblasNoTrans,1.0,inverse,inter,0.0,gslpoint);
  
  //change format of voxel coordinates:

  point.x=(float) gsl_vector_get(gslpoint,0);
  point.y=(float) gsl_vector_get(gslpoint,1);
  point.z=(float) gsl_vector_get(gslpoint,2);
  */
 
  point.x=(float) (gsl_matrix_get(inverse,0,0)*gsl_vector_get(inter,0)+gsl_matrix_get(inverse,0,1)*gsl_vector_get(inter,1)+gsl_matrix_get(inverse,0,2)*gsl_vector_get(inter,2));
  point.y=(float) (gsl_matrix_get(inverse,1,0)*gsl_vector_get(inter,0)+gsl_matrix_get(inverse,1,1)*gsl_vector_get(inter,1)+gsl_matrix_get(inverse,1,2)*gsl_vector_get(inter,2));
  point.z=(float) (gsl_matrix_get(inverse,2,0)*gsl_vector_get(inter,0)+gsl_matrix_get(inverse,2,1)*gsl_vector_get(inter,1)+gsl_matrix_get(inverse,2,2)*gsl_vector_get(inter,2));
  
  //gsl_vector_free(gslpoint);
  //gsl_vector_free(inter);

  if (inv->calculated==false)
    {
      //new_inv->calculated=true;
      inv->calculated=true;
      // new_inv->xform=inverse;
      inv->xform=inverse;
      //new_inv->inter=inter;
      inv->inter=inter;
      //inv=new_inv;
     
    }


  return point;

}


POINT3D WorldToVoxel(POINT3D inpoint, IMAGE_GEOM *image_geom)
{

  return WorldToVoxel(inpoint.x,inpoint.y,inpoint.z,image_geom);




 

}


VECTOR3D WorldToVoxel(VECTOR3D input_vector, IMAGE_GEOM *image_geom)
{
  /*IRL just rotate, do not scale.*/
  //rotation takes into account the steps, possibly negative or anisotropic

  VECTOR3D v;


  INVERSE *inv;
  gsl_vector *inter;
  //inter=gsl_vector_alloc(3);
  //gsl_vector *gslvector;
  //gslvector=gsl_vector_alloc(3);
  gsl_permutation *permutation;
  gsl_matrix *xform;
  gsl_matrix *inverse;

  if (image_geom->vector_inverse_xfm==NULL)
    {
      image_geom->vector_inverse_xfm=new INVERSE();
      
    }
  inv=image_geom->vector_inverse_xfm;

  if (inv->calculated==true)
    {
      inter=inv->inter;
     
    }
  else
    {
      inter=gsl_vector_alloc(3);
      
    }

  
  gsl_vector_set(inter,0,(double) input_vector.x);
  gsl_vector_set(inter,1,(double) input_vector.y);
  gsl_vector_set(inter,2,(double) input_vector.z);
  
  if (inv->calculated==false)
    {

  

      permutation=gsl_permutation_alloc(3);
      
      int p;

      //create xform 
      
      xform=gsl_matrix_alloc(3,3);

      gsl_matrix_set(xform,0,0,(double) image_geom->x_direction_cosine.x);
      gsl_matrix_set(xform,0,1,(double) image_geom->y_direction_cosine.x);
      gsl_matrix_set(xform,0,2,(double) image_geom->z_direction_cosine.x);

      gsl_matrix_set(xform,1,0,(double) image_geom->x_direction_cosine.y);
      gsl_matrix_set(xform,1,1,(double) image_geom->y_direction_cosine.y);
      gsl_matrix_set(xform,1,2,(double) image_geom->z_direction_cosine.y);

      gsl_matrix_set(xform,2,0,(double) image_geom->x_direction_cosine.z);
      gsl_matrix_set(xform,2,1,(double) image_geom->y_direction_cosine.z);
      gsl_matrix_set(xform,2,2,(double) image_geom->z_direction_cosine.z);



      
      inverse=gsl_matrix_alloc(3,3);
	
     

      gsl_linalg_LU_decomp(xform,permutation,&p);
      gsl_linalg_LU_invert (xform,permutation,inverse);
      gsl_matrix_free(xform);
      gsl_permutation_free(permutation);
     
    }
  else
    {
      //inverse=image_geom->vector_inverse_xfm.xform;
      inverse=inv->xform;
    }

  //gsl_linalg_LU_solve(xform,permutation,inter,gslvector);
  /*
  gsl_blas_dgemv(CblasNoTrans,1.0,inverse,inter,0.0,gslvector);


  
  //change format of voxel coordinates:


  vector.x=(float) gsl_vector_get(gslvector,0);
  vector.y=(float) gsl_vector_get(gslvector,1);
  vector.z=(float) gsl_vector_get(gslvector,2);
  */

  //just mult here:

  v.x=(float) (gsl_matrix_get(inverse,0,0)*gsl_vector_get(inter,0)+gsl_matrix_get(inverse,0,1)*gsl_vector_get(inter,1)+gsl_matrix_get(inverse,0,2)*gsl_vector_get(inter,2));
  v.y=(float) (gsl_matrix_get(inverse,1,0)*gsl_vector_get(inter,0)+gsl_matrix_get(inverse,1,1)*gsl_vector_get(inter,1)+gsl_matrix_get(inverse,1,2)*gsl_vector_get(inter,2));
  v.z=(float) (gsl_matrix_get(inverse,2,0)*gsl_vector_get(inter,0)+gsl_matrix_get(inverse,2,1)*gsl_vector_get(inter,1)+gsl_matrix_get(inverse,2,2)*gsl_vector_get(inter,2));
 
  
  //gsl_vector_free(gslvector);
  //gsl_vector_free(inter);


  if (inv->calculated==false)
    {
      inv->calculated=true;
      inv->xform=inverse;
      inv->inter=inter;

      //image_geom->vector_inverse_xfm.gslpoint=gslvector;
      
    }

  //account for steps:


  v.x=v.x/image_geom->x_step;
  v.y=v.y/image_geom->y_step;
  v.z=v.z/image_geom->z_step;

  v=Normalize(v);
 

  return v;
}



//this is normal curve centered at zero (could overload with 3-input nonzero centre)
float normal(float x, float sigma)
{
  return (1.0/(sigma*sqrt(2*M_PI)))*exp(-x*x/(2*sigma*sigma));
}

//this means max is one always
float normal_notnorm(float x, float sigma)
{
  return exp(-x*x/(2*sigma*sigma));
}

//this gives normal distribution value, normalized to unit volume assuming *truncated* at limit
float normal(float x, float sigma, float limit)
{
  //will get "nan" if sigma_theta==0 - then check for that later

  float index_inf=FLT_MAX;

  if (sigma<FLT_EPSILON)
    {
      return index_inf;
    }

  return (1.0/(sigma*sqrt(M_PI)*gsl_sf_erf(M_PI/2)))*exp(-x*x/(2*sigma*sigma));



}

float Factorial(int n) {
  float res = 1;
  for(int i=2; i<=n; i++)
    res *= i;

  return res;
}


PolarVector3D Cartesian2Spherical(Vector3D v) {
  float r, phi, theta;
  PolarVector3D s;

  r = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
  /*
    theta computation.
  */
  if(r == 0) {
    cout << "Cannot have a 0 radius in spherical coordinates!\n";
    exit(1);
  }

  theta = acos(v.z/r);  /*this value is between 0 and PI*/

  /*
    phi computation.

    atan2 is the same as atan(y/x) but looks at the sign of x and y
    to determine the quadrant of the answer. It returns a value between
    -PI and PI.
  */
  phi = atan2(v.y, v.x);  /*this value is between 0 and PI*/

  s.r = r;
  s.phi = phi;
  s.theta = theta;

  return s;
}

std::complex<float> GetSH(int _l, int _m, float theta, float phi) {
  int absm = std::abs(_m);
  float sign;

  if (absm%2==1)
    sign=-1.0;
  else
    sign=1.0;

  std::complex<float> retval(0.0,(float)(absm*phi));
  retval = std::exp(retval);
  
  float factor = sqrt(((float)(2*_l+1) / (4.0*M_PI))*(Factorial(_l - absm) / Factorial(_l + absm)))*gsl_sf_legendre_Plm(_l,absm,cos(theta));

  retval = factor*retval;

  if (_m<0) {
    retval = std::conj(retval);
    retval = sign*retval;
  }

  return retval;
}



gsl_matrix *ComputeSHMatrix(Directions *directions, int order) {
  //cout << "Computing SH basis of ";    
	

  int rank=(order+1)*(order+2)/2;
  int n_s=directions->GetNumberOfDirections();

  /*
    We declare the Bmatrix of size rank x n_s.  (????)
  */
  //gsl_matrix *B = gsl_matrix_alloc(rank,n_s);
  gsl_matrix *B = gsl_matrix_alloc(n_s,rank);
  
  float phi, theta;
  std::complex<float> cplx, cplx_1, cplx_2;

  float p[3];

  

  Vector3D d;
  PolarVector3D d_spherical;

  for(int i = 0; i < n_s; i++) {
    int j = 0;   //counter for the j dimension of B
        
    //get directions in Cartesian coordinates

 
    d.x = directions->GetXComponent(i);
    d.y = directions->GetYComponent(i);
    d.z = directions->GetZComponent(i);

    //get the direction in spherical coordinates
    d_spherical = Cartesian2Spherical(d);

    //DO:check that theta and phi correct, not swapped:    
    phi   = d_spherical.phi;
    theta = d_spherical.theta;
        
    /*
      Generate the modified Spherical Harmonic basis
      It is designed to be Real with standard SH properties such
      as orthonormality
      See Research Report Descoteaux RR 5681      
    */
    int sign;
    for(int l = 0; l <= order; l+=2)
      for(int m = 0; m <= l; m++) {
                
	/* positive "m" SH */
	if(m == 0) {
	  cplx_1 = GetSH(l,m,theta, phi);

	  if(fabs(imag(cplx_1)) > 0.0001) {
	    cout << "Modified spherical harmonic basis must be REAL!!!\n";
	    exit(0);
	  }
					
	  //B(j,i) = real(cplx_1);
	  //B[j][i] = real(cplx_1);
	  gsl_matrix_set(B,i,j,real(cplx_1));
	}
	else{
	  /* negative "m" SH  */
	  //get the corresponding spherical harmonic
	  cplx_1 = GetSH(l, m, theta, phi);
	  cplx_2 = GetSH(l,-m, theta, phi);

	  /* (-1)^m */
	  if(m % 2 == 1)
	    sign = -1;
	  else
	    sign = 1;
	  {
	    std::complex<float> s(sign, 0.0);
	    std::complex<float> n(sqrt((float)2)/2, 0.0);
	    cplx = n*(s*cplx_1 + cplx_2);
	  }
	  /* (-1)^(m+1) */
	  if(m % 2 == 1)
	    sign = 1;
	  else
	    sign = -1;
	  {
	    std::complex<float> n(0.0, (float)sqrt((float)2)/2);
	    std::complex<float> s(sign, 0.0);
	    cplx_2 = n*(s*cplx_1 + cplx_2);
	  }
	  if(fabs(imag(cplx)) > 0.0001 ||
	     fabs(imag(cplx_2)) > 0.0001) {
	    cout << "Modified spherical harmonic basis must be REAL!!!\n";
	    exit(0);
	  }
	  //B(j,i) = real(cplx);
	  //B[j][i] = real(cplx);
	  gsl_matrix_set(B,i,j,real(cplx_1));//JC changed i,j

	  j+=1;
	  //B(j,i) = real(cplx_2);
	  //B[j][i] = real(cplx_2);
	  gsl_matrix_set(B,i,j,real(cplx_2));//JC changed i,j
	}
	j++;
      }
  }

  
  //cout << "Done.\n";	
  /* pointer to first element of the matrix */
  //return (float *) gsl_matrix_ptr(B,0,0);	
  return B;
}


float Legendre0(const int & order) {
  if(order == 0)
    return 1.0;
  if(order % 2 != 0)
    return 0.0;
  else {
    float odd = 1;
    float even = 2;

    for(int i = 3; i <= order - 1; i+=2) {
      odd = odd * i;
    }
    for(int i = 4; i <= order; i+=2) {
      even = even * i;
    }
		
    if((order / 2) % 2 == 0)
      return odd / even;
    else
      return -1 * odd / even;
  }
}


double DeltaAngle(double sigma)
{
  // using namespace boost::random;

  boost::mt19937 rngn(42u);
  boost::normal_distribution<double> gaussian(0.0, sigma);

  int seed;
  seed = (int) time (NULL);
  srand48(seed);
  rngn.seed(static_cast<unsigned int>((double) rand()));

  boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> >  gauss_num_gen(rngn, gaussian);

  double rangle;
  rangle=gauss_num_gen();

  return rangle;
}

double double_abs(double a)
{
  if (a<0)
    {
      return -a;
      
    }
 else
   {
     return a;
     
   }
  
}

float ScalarForRGB(float* RGB,vtkLookupTable *lut)
{
  double rgb[3];
  rgb[0]=(double) RGB[0];
  rgb[1]=(double) RGB[1];
  rgb[2]=(double) RGB[2];


  return ScalarForRGB(rgb,lut);

  
}

float ScalarForRGB(VECTOR3D RGB,vtkLookupTable *lut)
{
  double rgb[3];
  rgb[0]=(double) RGB.x;
  rgb[1]=(double) RGB.y;
  rgb[2]=(double) RGB.z;


  return ScalarForRGB(rgb,lut);

  
}


float ScalarForRGB(double* RGB,vtkLookupTable *lut)
{
  //go through and find what scalar maps to this RGB value
 
  double scalarstep;
  double s;
  double rgb[3];
  //get the range of the lut and determine scalarstep from it:
  double *range;
  double norm;
  float max_dot;
  double max_s;
  double RGBnew[3];

  //consider only abs val of components:
	  if (RGB[0]<0)
	    {
	      RGBnew[0]=-RGB[0];
	    }
	  else
	    {
	      RGBnew[0]=RGB[0];
	    }
	  if (RGB[2]<0)
	   
	    {
	      RGBnew[2]=-RGB[2];
	    }
	  else
	    {
	      RGBnew[2]=RGB[2];
	    }
	  if (RGB[1]<0)
	    
	    {
	      RGBnew[1]=-RGB[1];
	    }
	  else
	    {
	      RGBnew[1]=RGB[1];
	    }
	  


  range=lut->GetRange();
  scalarstep=(range[1]-range[0])/lut->GetNumberOfTableValues();

  
  max_dot=0;
  for (s=range[0];s<=range[1];s+=scalarstep)
    {
      lut->GetColor(s,rgb);
      norm=sqrt(rgb[0]*rgb[0]+rgb[1]*rgb[1]+rgb[2]*rgb[2]);
      rgb[0]=rgb[0]/norm;
      rgb[1]=rgb[1]/norm;
      rgb[2]=rgb[2]/norm;
      
      //look for closest:
      if (float_abs((float) (rgb[0]*RGBnew[0]+rgb[1]*RGBnew[1]+rgb[2]*RGBnew[2]))>max_dot)
	{
	  max_dot=float_abs((float) (rgb[0]*RGBnew[0]+rgb[1]*RGBnew[1]+rgb[2]*RGBnew[2]));
	  max_s=s;
	}
      

    }


  return max_s;


}

ImageData *Create100DirInput(ImageData *lowdirinput)
{
  //interpolate the old input data.  we will output the higher res ODF.
  Directions *directions=new Directions("100dirs.txt");
  ImageData *newdata=new ImageData(DWIdata,lowdirinput->GetXSize(),lowdirinput->GetYSize(),lowdirinput->GetZSize(),100,lowdirinput->GetXStart(),lowdirinput->GetYStart(),lowdirinput->GetZStart(),0,lowdirinput->GetXStep(),lowdirinput->GetYStep(),lowdirinput->GetZStep(),1,lowdirinput->GetXDirectionCosine(),lowdirinput->GetYDirectionCosine(),lowdirinput->GetZDirectionCosine());
  newdata->SetDirections(directions);

  S2Interpolator *interpolator=new S2Interpolator(lowdirinput, directions, 1.5*sqrt(2*M_PI/lowdirinput->GetNSize()), 1.5*sqrt(2*M_PI/lowdirinput->GetNSize()), Gaussian);

  int dir,x,y,z;

  for (dir=0;dir<newdata->GetNSize();dir++)
    {
      for (x=0;x<newdata->GetXSize();x++)
	{
	  for (y=0;y<newdata->GetYSize();y++)
	    {
	      for (z=0;z<newdata->GetZSize();z++)
		{
		  newdata->SetValue(x,y,z,dir,interpolator->GetInterpolatedValue(dir,x,y,z));
		}
	    }
	}
    }

 

  DiffusionParameters *diffusion_parameters=new DiffusionParameters(newdata->GetNSize());
  newdata->SetDiffusionParameters(diffusion_parameters);
  newdata->GetDiffusionParameters()->SetConstb(lowdirinput->GetDiffusionParameters()->GetConstb());
  newdata->GetDiffusionParameters()->SetDirections(directions);
  return newdata;

}
