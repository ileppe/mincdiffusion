/************************************************************************

   File: misc.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __misc_h
#define __misc_h

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <complex>
#include <float.h>
#include <math.h>
#include <string>

#include "otherlibs_include.h"

/*
note: the WorldToVoxel / VoxelToWorld tranforms for vectors take into account the possibility of negative and/or anisotropic steps

note that the pointers in INVERSE and IMAGE_GEOM never get freed automatically: user has to once done  In these programs there is generally only one place where this needs to be done.

*/

class ImageData;
class Directions;

struct Point3D 
{
  float x;
  float y;
  float z;
  
  Point3D(float _x,float _y,float _z):
   x(_x),y(_y),z(_z)
  {}
  
  Point3D()
  {}
};


inline double euclidian_distance(const Point3D& a,const Point3D& b)
{
  return sqrt((a.x-b.x)*(a.x-b.x) + 
              (a.y-b.y)*(a.y-b.y) + 
              (a.z-b.z)*(a.z-b.z));
}

struct Vector3D 
{
  float x;
  float y;
  float z;
  
  Vector3D(float _x,float _y,float _z):
   x(_x),y(_y),z(_z)
  {}
  
  Vector3D()
  {}
  
} ;

struct PolarVector3D 
{
  float r;
  float phi;
  float theta;
} ;

struct Index3D 
{
  int x;
  int y;
  int z;
  
  Index3D() {}
  
  Index3D(int _x,int _y,int _z):x(_x),y(_y),z(_z)
  {}
} ;

inline bool operator==(const Index3D& a,const Index3D& b)
{
  return a.x==b.x &&
         a.y==b.y &&
         a.z==b.z;
}

struct Inverse
{
  
  gsl_matrix *xform;
  bool calculated;

  //these are for reuse, so we don't alloc/free every time
  gsl_vector *inter;
  //gsl_vector *gslpoint;
  Inverse()
  {
    calculated=false;
  }
  
};

typedef Inverse INVERSE;

struct image_geom 
{
  int x_size;
  int y_size; 
  int z_size; 
  int n_size; 
  
  float x_start; 
  float y_start; 
  float z_start; 
  float n_start; 
  
  float x_step; 
  float y_step; 
  float z_step; 
  float n_step; 
  
  Vector3D x_direction_cosine; 
  Vector3D y_direction_cosine; 
  Vector3D z_direction_cosine;

  INVERSE *vector_inverse_xfm;
  INVERSE *point_inverse_xfm;




};




typedef image_geom IMAGE_GEOM;
typedef Point3D POINT3D;
typedef Index3D INDEX3D;
typedef PolarVector3D POLARVECTOR3D;
typedef Vector3D VECTOR3D;


void Error(const char *msg);


char *CreateHistoryString(int argc, char *argv[]); //ilana modify this to use sstream

void ResampleStandardDirCos(char *in_filename,char *out_filename,char *tempdirname,bool nearest);

void Resample(std::string & in_filename,std::string & out_filename,std::string & tempdirname,float step,bool use_input_sampling,std::string &templatefilename, bool& templatefilename_set,bool nearest);

void Resample(char* in_filename,char* out_filename,char* tempdirname,float step,bool use_input_sampling,char*templatefilename, bool& templatefilename_set,bool nearest);


void ResamplePosStep(char *in_filename,char *out_filename,char *tempdirname);

void ResampleNegStep(char *in_filename,char *out_filename,char *tempdirname);

float float_min3(float a,float b,float c);

float float_max3(float a,float b,float c); 

float float_max(float a, float b);

float float_min(float a, float b);

float float_abs(float a);

void MakeTwoOrthogonalVectors(float vector[3], float *basis1, float *basis2);

void MakeTwoOrthogonalVectors(VECTOR3D vector, VECTOR3D &basis1, VECTOR3D &basis2);

int int_max3(int a, int b, int c);

int int_min3(int a,int b,int c);

int int_max(int a, int b);

int int_min(int a, int b);

int int_abs(int a);

VECTOR3D cross(VECTOR3D v1, VECTOR3D v2);

float dot(VECTOR3D v1, VECTOR3D v2);

POINT3D VoxelToWorld(float x,float y,float z, IMAGE_GEOM image_geom);

POINT3D VoxelToWorld(POINT3D inpoint, IMAGE_GEOM image_geom);

VECTOR3D VoxelToWorld(VECTOR3D input_vector, IMAGE_GEOM image_geom);

//POINT3D WorldToVoxel(float x,float y,float z, IMAGE_GEOM &image_geom);

/*
POINT3D WorldToVoxel(float x,float y,float z, IMAGE_GEOM image_geom);
POINT3D WorldToVoxel(float x,float y,float z, IMAGE_GEOM image_geom, INVERSE *inv);


POINT3D WorldToVoxel(POINT3D inpoint, IMAGE_GEOM &image_geom);
POINT3D WorldToVoxel(POINT3D inpoint, IMAGE_GEOM &image_geom, INVERSE *inv);


VECTOR3D WorldToVoxel(VECTOR3D input_vector, IMAGE_GEOM &image_geom);
VECTOR3D WorldToVoxel(VECTOR3D input_vector, IMAGE_GEOM &image_geom, INVERSE *inv);
*/
POINT3D WorldToVoxel(float x,float y,float z, IMAGE_GEOM *image_geom);
POINT3D WorldToVoxel(POINT3D inpoint, IMAGE_GEOM *image_geom);
VECTOR3D WorldToVoxel(VECTOR3D input_vector, IMAGE_GEOM *image_geom);

float normal(float theta, float sigma_theta);

float normal(float x, float sigma, float limit);

float normal_notnorm(float x, float sigma);

float Factorial(int n);

PolarVector3D Cartesian2Spherical(Vector3D v);

gsl_matrix *ComputeSHMatrix(Directions *directions, int order);

std::complex<float> GetSH(int _l, int _m, float theta, float phi);

float Legendre0(const int & order);

float flt_norm(VECTOR3D vector);

VECTOR3D Normalize(VECTOR3D vector);


double DeltaAngle(double sigma);

float ScalarForRGB(float* RGB,vtkLookupTable *lut);
float ScalarForRGB(double* RGB,vtkLookupTable *lut);
float ScalarForRGB(VECTOR3D RGB,vtkLookupTable *lut);
 
double double_abs(double a);

VECTOR3D neg(VECTOR3D v);

ImageData *Create100DirInput(ImageData *lowdirinput);

#endif
