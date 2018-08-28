/************************************************************************

   File: ImageData.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   DescriptI think sigma_theta.number is never set and that is whyion: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __ImageData_h
#define __ImageData_h 

#include "otherlibs_include.h"
#include "DiffusionParameters.h"
#include "Directions.h"
#include "SphericalHarmonics.h"
#include "misc.h"

#include <minc_1_simple_rw.h>

/*

the ProbLabel file gets read in as a Maxima3 because we haven't written an attribute to distinguish it

the constructors need to be cleaned up so that all common stuff appears once in code; began with Construct(): more to be done

(t),z,y,x so that minc utilities (vis) such as xdisp work with the outputted volume

currently float only

as is, creating from volume means a volume copy must be performed (I think we want it this way...)

array saved as a gsl vector (handy for file IO and finding min/max, error handling, copying, math ops (***), initializing, etc

a 3rd constructor that has a parent ImageData as a template (for certain things but not values in the array) could be called for; also need a copy constuctor

DO: copy constructor!!!!!!

DO: volume_exists and all other public variables should be private and their values should be acccessible by functions!!!!  as is, can really mess things up by changing them
when the data structures don''t change

should check for diffusion_direction_x(y,z,) attributes, but for backwards compatibility should allow for an extra file to be specified and for no diffusion information to be present in the minc file 

check bayes.c for how to have full word options in command line (use ParseArgv).  then specify all input files for a program with, e.g., -e1x -e1y -e1z, etc.

for eigenvectors:should really have a vector ImageData with a 4D array.  but then FileIO gets tricky: require multiple files or vector files for a single ImageData.  As is, keeping eigenvector components separate, as separate ImageDatas, and ProcessedDiffusion.... has them all.


need to think about how masking, etc, will work.  checking for whether the minc coords are the same?  Or would user sometimes not want that?  Safe to assume the files and therefore the ImageDatas match?  the calling program can take care of that with mincresample.  

should have backdoor access to the array (and any pertinent information) prototype there but not implemented

should have a voxel-loop type method for doing something to all voxels

DWIdata image type could be 3D or 4D.  requires directions and bvalues
information. ODFdata requires directions (only). spherical harmonics
may require other things.  should perhaps inherit wherever these objects start to differ in their attributes.  not sure.  set all attributes to something like NULL before
(if inherit - make a subclass DWIImageData that reads diffusion parameters.  curently not doing this)

could be smart to have a lock on certain members (that are pointers) to say "someone else is using this" e.g. ....

can do reference counting

add: if a ODFData, can have just Directions member (instead of DiffusionParameters)

could have GetVolume() if want to use other volume_io features and explicitly have a volume.  currently the volume is hidden (could be stored in any way)

constructors for creating from scratch will have to handle telling it to create a Directions or a DiffusionParameters, etc.  image type passed as arg.  Means you can''t change the image_type after.  must have this info at time of creation.

currently will require steps and starts (and ...(?)) to be given at time of creation from scratch.  could have some defaults and add Set methods for these.  
prototype for default constructor is here...

could Set starts and steps but must change the volume itself if change sizes (prob can''t easily)

need a method for outputing a modified volume with output_modified_volume: so have the header information.  currently lose that because that info is not written out.  would have to OutputMincFileWithOrigHeaderInfo (if read in OR copied from another ImageData).  as is, an ImageData does have a header so expect nothing.  maybe a minc utility could add the header info afterward?  nah.  so need to have a header associated with a copied or read-in ImageData.

mean want other stats such as mean, etc.  maybe these things should start being in another object that calculates stats, in an roi if desired.

   not necessary if using volume structure:  might want anyway?
  int _xsize;
  int _ysize;
  int _zsize;
  int _nsize;
  float _xstep;
  float _ystep;
  float _zstep;
  float _nstep;//irrelevant for ODF data but relevant for time series
  float _xstart;
  float _ystart;
  float _zstart;
  float _nstart;//irrelevant for ODF data but relevant for time series
  

DO: there are way too many flags here.  get rid of them and assume things went right.  or actually use them. i.e. don''t GetXComponent if not allocated...
get rid of all except need array_exists : others?  to inconsistant: does exist mean exist or initialized? etc

DO: need to figure out how to get specific header info, but not all, from a 4D file (eg DWIdata or ODFdata) to a 3D file (eg connectivity values)

currently no set method for image type: do we want it?  prob not.

could force positive step signs by mincreshape +direction -> temp dir beforehand

Image3D_with_diffusion_parameters is not implemented

DO: add a method to output geom parameters and can use this struct when initializing a new ImageData, or in FibreTractSet, etc.

INitialize to 0 should be here!

set all voxels is only done for 3D

maybe adds methods to init with IMAGE_GEOM and return IMAGE_GEOM, but issues with how many dims, etc: think about it.

may need to write a wrapper around output min cfile so that it outputs slice by slice and concats for 4D files: getting strange behaviour

should ImageData be more arbitrary? hold any sort of struct?

GetMinimumNonZeroValue(int n); is specifically for 4D.. and it assumes positive values!!!!!!

and getmax is too.

-------

WriteODFMaxima: file format:  

(binary): in order zyx always!!!  sampling is from the ODF file sampling (the maxima file accompanies it; should be attached somehow!!

for each voxel:
#maxima max_x max_y max_z. 

The maxima are in voxel space (use direction cosines of ODF file to get world space vectors)

-------
    for set/get ODFMaxima, must have calculated already to write: else write 0 dirs... 

an ODFMaxima file is linked specifically to the ODF ImageData''s file, but out in the file system, it is up to the user to keep this straight.  There is no way to resample the maxima to go along with a resampling of an ODF, so 

ODF maxima are in *voxel space*

Maxima3 is here for backward compatibility only; phase out

 */


class S2Interpolator;
class SphericalHarmonics;

typedef enum image_type {DWIdata, ODFdata, Image3D_vanilla, Image3D_with_diffusion_parameters, Maxima3, ProbLabel} IMAGE_TYPE;


struct eNeighbourLUTEntry
{
  int index;
  struct eNeighbourLUTEntry *next;

} ;

typedef eNeighbourLUTEntry NeighbourLUTEntry;


struct ODF_Maxima
{
  bool calculated;
  Directions *maxima;
} ;

typedef ODF_Maxima ODF_MAXIMA;

//typically in degrees
struct Sigma_theta
{
  float *sigma_theta;
  int number;
  bool calculated;//this probably isn't necessary: just have number be 0
} ;

typedef Sigma_theta SIGMA_THETA;

struct Confid_Value
{
  float *confid;
  int number;
  bool calculated;
};

typedef Confid_Value CONFID_VALUE;


class ImageData
{
  typedef minc::fixed_vec<3,double> _vector;

 public:
  //ImageData::ImageData(); 
  #if 0
  ImageData(IMAGE_TYPE image_type, int x_size, int y_size, int z_size);
  
  ImageData(IMAGE_TYPE image_type, int x_size, int y_size, int z_size, 
            int n_size);
            
            
            
  #endif           
  //ImageData::ImageData(ImageData* template_image_data);
  //ImageData::ImageData(const ImageData &image_data);
  //ImageData::ImageData(volume_io::Volume image_minc_volume);
  //ImageData::ImageData(volume_io::Volume *image_minc_volume, DiffusionParameters *diffusion_parameters);
  ImageData();
  
  ImageData(const char* minc_filename);
  
  ImageData(IMAGE_TYPE image_type, 
            int x_size,    int y_size,    int z_size,    int n_size, 
            float x_start, float y_start, float z_start, float n_start, 
            float x_step,  float y_step,  float z_step,  float n_step, 
            const VECTOR3D& x_direction_cosine,  const VECTOR3D& y_direction_cosine, const VECTOR3D& z_direction_cosine);
            
  ImageData(IMAGE_TYPE image_type, int x_size, int y_size, int z_size, 
            float x_start, float y_start, float z_start, 
            float x_step, float y_step, float z_step, 
            const VECTOR3D& x_direction_cosine,  const VECTOR3D& y_direction_cosine, const VECTOR3D& z_direction_cosine);
            
  ImageData(IMAGE_TYPE image_type, const IMAGE_GEOM& image_geom);
  //ImageData::ImageData(char* minc_filename,char* diffusion_parameters_filename);

  float GetValue(int x,int y,int z);
  float GetValue(int x,int y,int z,int n);
  void  SetValue(int x,int y,int z,float value);
  void  SetValue(int x,int y,int z,int n,float value);
  
  void LoadMincFile(const char *minc_filename);
  void WriteMincFile(const char* minc_filename, bool clobber);
  
  //void ImageData::WriteMincFile(char* minc_filename, bool clobber, char *history_string);
  void WriteMincFile(const char* minc_filename, bool clobber, const char *history_string, bool bbox); //IL added option for bounding box
  void WriteMincFileNoRange(const char* minc_filename, bool clobber, const char *history_string);
  
  int GetDimensions();
  int GetXSize();
  int GetYSize();
  int GetZSize();
  int GetNSize();
  double GetXStep();
  double GetYStep();
  double GetZStep();
  double GetNStep();//irrelevant for ODF data but relevant for time series
  int GetXStepSign();
  int GetYStepSign();
  int GetZStepSign();
  int GetNStepSign();//irrelevant period
  double GetXStart();
  double GetYStart();
  double GetZStart();
  double GetNStart();//irrelevant for ODF data but relevant for time series
  VECTOR3D GetXDirectionCosine();
  VECTOR3D GetYDirectionCosine();
  VECTOR3D GetZDirectionCosine();
  VECTOR3D GetVoxelFromWorldCoordinates(double x_world, double y_world, double z_world);
  DiffusionParameters *GetDiffusionParameters(); 
  void SetDiffusionParameters(DiffusionParameters *diffusion_parameters);//use with caution
  Directions *GetDirections();
  void SetDirections(Directions *directions);
  IMAGE_TYPE GetImageType();  
  float GetMaximumValue();
  float GetMaximumValue(float exclude);
  float GetMaximumValue(int x,int y, int z);
  float GetMaximumValue(int x, int y, int z, int &index);
  float GetMaximumValue(int x, int y, int z, int &index, int start, int stop);
  float GetMinimumValue();
  float GetMinimumNonZeroValue(int x, int y, int z);
  float GetMinimumNonZeroValue(int x, int y, int z, int &index);
  float GetNumberOfVoxelsAbove(float threshold=0.0);
  void SetAllVoxels(float value);
  void Retain();
  void Release();
  Directions *GetODFMaxima(int x,int y,int z);
  float GetODFStdDev(int x,int y,int z);
  float GetODFMean(int x,int y,int z);
  void WriteODFMaxima(const char *filename, bool clobber, char *history);
  void SetODFMaxima(const char *filename);
  void SetODFMaxima(Directions *directions,int x,int y,int z);
  void SetSigmaTheta(SIGMA_THETA sigma_theta,int x,int y,int z);
  SIGMA_THETA GetSigmaTheta(int x,int y,int z);
  void SetSigmaTheta(const char *filename);
  void WriteSigmaTheta(const char *filename, bool clobber, char *history);
  void NormalizeODF(int voxel_x, int voxel_y, int voxel_z);
  IMAGE_GEOM GetImageGeom();

  //DO: Set method where appropriate, prob to set data after get it (?)

  //void ImageData::InitializeToZero(); 
  
  ImageData *CreateDWIODF();
  Directions *GetODFMaximaFast(int x,int y,int z);
  Directions *GetODFMaximaHR(int x,int y,int z);


  bool SigmaThetaAvailable();


  void WriteSigmaTheta(const char *filename, bool clobber, char *history, bool dble);
  void WriteSigmaTheta(ofstream& fd, bool clobber, char *history);
  void WriteODFMaxima(ofstream& fd, bool clobber, char *history);

  int FindClosestDirection(VECTOR3D TanVec);

 protected:
  ~ImageData();

 private:    

  bool _image_array_is_current; 
  bool _image_array_exists;
  gsl_vector_float *_image_array;  
  minc::simple_4d_volume<float> _image_minc_volume;
  IMAGE_TYPE _image_type;
  DiffusionParameters *_diffusion_parameters; //for use if _image_type==DWIdata
  Directions *_directions; //for use if _image_type==ODFdata
  bool     _direction_cosines_exist;
  VECTOR3D _x_direction_cosine;
  VECTOR3D _y_direction_cosine;
  VECTOR3D _z_direction_cosine;
  int _reference_count;
  S2Interpolator *_ODF_max_interpolator;
  bool _ODF_max_interpolator_set;
  ODF_MAXIMA *m_maxima_array;
  bool m_maxima_array_exists;
  SIGMA_THETA *m_sigma_theta_array;
  bool m_sigma_theta_array_exists;
  NeighbourLUTEntry **m_neighbour_LUT;
  bool m_neighbour_LUT_calculated;
  bool _directions_are_set;
  Directions *m_ODF_sample_dirs;
  int *m_ODF_sample_neighbours;
  S2Interpolator *m_ODF_interpolator;
  SphericalHarmonics *m_sphere_harm;
  bool m_sphere_harm_constructed;
  ImageData *m_SH_fit;

  void _PutVolumeInArray();
  bool _bvaluesAttributeExists(minc::minc_1_reader& rdr);
  bool _directionsAttributeExists(minc::minc_1_reader& rdr);
  void CreateNeighbourLUT();
  void AddNeighbour(int entry, int neighbour);
  
  void Construct(IMAGE_TYPE image_type, 
                 int x_size,    int y_size,    int z_size,    int n_size, 
                 float x_start, float y_start, float z_start, float n_start, 
                 float x_step,  float y_step,  float z_step,  float n_step, 
                 const VECTOR3D& x_direction_cosine,  const VECTOR3D& y_direction_cosine, const VECTOR3D& z_direction_cosine);
                 
  float m_maximum_cutoff_fraction;
  
// method added by Parya
  bool _LabelProbAttributeExists(char *minc_filename);

  IMAGE_GEOM m_image_geom;
  bool m_image_geom_init;

};

#endif //__ImageData_h
