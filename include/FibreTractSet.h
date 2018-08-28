/************************************************************************

   File: FibreTractSet.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __FibreTractSet_h
#define __FibreTractSet_h 

#include <vector>

#include "otherlibs_include.h"
#include "FibreTract.h" 
#include "ImageData.h"
#include "misc.h"
#include <stdio.h>

/*

currently attempting JUST CumMinOTT, (w/ tracts debatable)

DOTT needs an overhaul: DOTT itself could be useful (not implemented); minDOTT, as minOTT, probably not, and would need CumMinDOTT...

for template file: need to implement copy constructor of ImageData.  making it without all the header info for now: just geometrical info

FibreTractSet copy constructor would be used by ...

SetScalarOutputMode assumes that scalar exists for all tracts in set

DO: instead of writing "index" generically to the file, write exactly what the index is: connectivity index, minDOTT, meanDOTT


DO: change so can write as many scalars as you choose. ie. scalar_ouput_mode is not exclusive

NOTE: the new FibreTractSet created with Prune() has pointers to the original tracts: they are not copied, and hence if changed, are changed in all Sets that point to them..

add method to set connectivity index type should this become necessary

in the Get..Map methods, could compute if not done already (have a flag)

change so scalar output mode can be a list and write to both the tracts file and minc files for each: for now, just tracts file

need to finish the misc scalars option by adding it to FibreTract and to outputvtk..

Q allow multiple scalars to be written at once?

only cum..minOTT is implemented for Prune..

this should use a different sort of data structure that allows removal of any tract easily (instead of creating a new set when prune)

NOTE!! destructor deletes the imagedatas: shouldn't assign other pointers to point to them: if you do: they're gone when this is: use Retain() where necessary (probably need to change this)

connectivity index isn't done in most cases

  //do this only if it makes sense (e.g., have all subsets):
  //fibre_tract_set->SetScalarOutputMode(CI [or min OTT,DOTT..] only makes sense for an entire tract, and as part of a set, only makes sense if have the subsets.  so not really using it much);

NOTE: can't prune, and hence can't do multi-roi, if have tracts off.  makes sense though: can we really do multi-roi for probabilistic?  thinking about what we really want here



*/

enum tract_scalar_output_mode {no_tract_scalar, CI, minDOTT, meanDOTT, minAI, meanAI, minOTT, misc_tract_scalar} ;
enum point_scalar_output_mode {no_point_scalar, cumulative_CI, cumulative_min_DOTT, DOTT, cumulative_min_OTT, OTT, misc_point_scalar} ;

typedef tract_scalar_output_mode TRACT_SCALAR_OUTPUT_MODE;
typedef point_scalar_output_mode POINT_SCALAR_OUTPUT_MODE;

class FibreTractSet
{
 public:
  FibreTractSet();
  FibreTractSet(char *filename);//doesn't seem to work at this point
  //copy constructor 
  //FibreTractSet(const FibreTractSet &fibre_tract_set);
  void ReInitialize();
  ~FibreTractSet();
  //void AddTract(FibreTract *fibre_tract);
  void AddTract(FibreTract *fibre_tract);
  
  int GetNumberOfTracts()
  {
    return _fibre_tract_array.size();
  }

  
  void OutputVTKAscii(); // IL added
  void OutputVTKFile(const char *filename);
  void OutputVTKFile(const char *filename, bool clobber);
  void OutputOBJFile(const char *filename, bool clobber);
  void SetImageDataTemplate(ImageData *image_data);
  
  void SetImageDataTemplate(const char *image_filename);
  
  void SetImageDataParameters(int x_size,    int y_size,    int z_size, 
                              float x_start, float y_start, float z_start,
                              float x_step,  float y_step,  float z_step, 
                              const VECTOR3D& x_direction_cosine, 
                              const VECTOR3D& y_direction_cosine, 
                              const VECTOR3D& z_direction_cosine);  
                              
  void SetImageDataParameters(IMAGE_GEOM image_geom);



  //these all need an overhaul: -------
  //need to amalgamate and treat as same here.  could set type in FibreTracking instead.  except in cases where might want more than one of them: could allow for that in amalgamation

  //this is a currently undefined Connectivity Index; using CumMinOTT from edge of ROI1 right now.  
  //void InitConnectivityIndexMap();
  //void ComputeConnectivityIndexMap();
  //ImageData *GetConnectivityIndexMap();
  //void SetConnectivityIndexMap();

  //in addition to these, could do MinAI, MeanAI
  //void InitMinDOTTMap();
  //void ComputeMinDOTTMap();
  //ImageData *GetMinDOTTMap();
  //void SetMinDOTTMap(ImageData *min_DOTT_map);

  //void InitMeanDOTTMap();
  //void ComputeMeanDOTTMap();
  //ImageData *GetMeanDOTTMap();
  //void SetMeanDOTTMap();

  //these are the most likely to be useful, although MinOTT is doubtful without predinfed edges, in which case, CumMinOTT has all the info
  //void InitMinOTTMap();
  //void ComputeMinOTTMap();
  //ImageData *GetMinOTTMap();
  //void SetMinOTTMap();

  void InitCumMinOTTMap();
  void ComputeCumMinOTTMap();
  ImageData *GetCumMinOTTMap();
  void SetCumMinOTTMap(ImageData *cum_min_OTT_map);

  //void InitOTTMap();
  //void ComputeOTTMap();
  //ImageData *GetOTTMap();
  //void SetOTTMap();



  //--------

  //these are "done":
  void InitBinaryConnectivityMap();
  void ComputeBinaryConnectivityMap(); //slow: usually done elsewhere
  ImageData *GetBinaryConnectivityMap();
  void SetBinaryConnectivityMap(ImageData *binary_connectivity_map);


  void SetTractScalarOutputMode(TRACT_SCALAR_OUTPUT_MODE tract_scalar_output_mode);
  void SetPointScalarOutputMode(POINT_SCALAR_OUTPUT_MODE point_scalar_output_mode);


  FibreTractSet *Prune(float threshold);//not implemented yet; for display purposes
  FibreTractSet *PruneSameIndices(bool compare_last_only);//works (IL)
  FibreTractSet *Prune(ImageData *roi, bool exclude, bool compare_last_only);


  //void PruneThisSet(ImageData *roi, bool exclude, bool compare_last_only);
  //FibreTractSet *Prune(bool compare_last_only, TRACT_SCALAR_OUTPUT_MODE tract_scalar_output_mode);

  FibreTract *GetTract(int n)
  {
    return _fibre_tract_array[n];

  }
  //SetTract(int n, FibreTract * fibre_tract); //for changing not at end

  void TractsOff();
  bool SaveTractsOn();
  void RemoveTract(int n);
  void SetReferenceROI(int reference_roi);


 private:
  std::vector<FibreTract *> _fibre_tract_array; //VF: it is 21st century
  //int _number_of_tracts;
  IMAGE_GEOM _output_geometry_parameters;
  ImageData *_binary_connectivity_map;
  //ImageData *_connectivity_index_map;
  //ImageData *_min_DOTT_map;
  //ImageData *_mean_DOTT_map;
  ImageData *_cum_min_OTT_map;
  //int _tracts_allocated;
  //int _tract_allocation_step;
  TRACT_SCALAR_OUTPUT_MODE _tract_scalar_output_mode;
  POINT_SCALAR_OUTPUT_MODE _point_scalar_output_mode;
  bool _binary_connectivity_map_init;
  //bool _min_DOTT_map_init;
  bool _cum_min_OTT_map_init;
  //void _ReallocateTractArray();
  bool _output_geometry_parameters_set;
  bool _save_tracts;
  bool _ASCII; // IL added for ASCII option
  FILE *_prune_dump;
};

#endif
