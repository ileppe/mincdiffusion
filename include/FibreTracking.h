/************************************************************************

   File: FibreTracking.h

   Author: Jennifer Campbell

   Created: 

   Revisions: 

   Description: 

   Copyright (c) 2004 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.  

***********************************************************************/
#ifndef __FibreTracking_h
#define __FibreTracking_h 

#include "otherlibs_include.h"
#include "FibreTractSet.h" 
#include "ImageData.h"
#include "misc.h"
#include "ProcessedDTIData.h"
#include "Directions.h"
#include "DiffusionODFCalculator.h"
#include "DisplayWindow.h"
#include "FibreTract.h" 


/*
bug alert: we are getting nan here  in AssignCumulativeMinOTTs sometimes ???? voxel lineup problem??
		       (float_min(CMOTT,tract->GetOTT(k))


note that all label tracking is prob. now, and that checking for m_ProbLabels!=NULL is the same as checking the algorithm: could clean up

if the requested curvature constraint is greater than 90, m_big_voxels gets set to true and the multi-step check is not performed.  this is a nonintuitive hack for a test of mine an isn't intended for anyone else's use

there is always a FibreTractSet regardless of whether we will save the tracts, ever.  it is used for the binary maps in all cases, and other maps in certain cases

the ability to do cumulative (rel to first ROI) or total min diffusion ODF tangent to tract is here for PDD/FACT, but not in the options of FibreTrack at this point: can add later and implement it given a *diffusion* ODF input for Maxima algorithm.


as is, prop. everywhere in cone with equal certainty; can implement e.g. Gaussian and have a non-binary connectivity profile (probably should)

(maybe different algorithms will end up subclassed...)


.Track: run the tracking, output the tracts and Lc (anything else?  speed function values...should be as minc file...)

eventually: ask Peter for estimates of how many curves/what curvature is in each voxel

->SetInterpolationType

->SetInputProcessedDTIData


I don't know how to determine curvature constraint as a function of the segments, which can have varying length, in FACT.  This is a problem.

should prob (yes) try Mori's implementation



mask can be anything: give a threshold and stop tracking if above/below the threshold.  so brain mask, AI mask, trace (e.g., <1, to exclude CSF)?.  this takes care of as many as you want, so you don't have to create a composite mask beforehand.  

allocating the roi_array one by one: can be changed

the curvature constraint is given by the radius of curvature, in mm.

as far as curvature constraint goes, voxels are assumed to be isotropic.  if that isn't the case, x is used.  the default is radius of curvature equal to the x step size

don't have to calculate max_angle every time. prob faster to do on the fly for variable, but could do ahead for whole volume.  do only once if not variabel

could have private pointer _this_fibre_tract instead of passing it around

DO: check StopTracking() again: when can check the angles and etc.  also check Prune() in FibreTrack: make sure that don't assign too many/few voxels given pointson edge

currently, _seeds_per_voxel is float so that divide ops return float: check for other such errors


minDOTTmapon requires all the elements of the tensor.  DOTT per point written to the FibreTractSet automatically in this case

ODF gives, by default, the min ODF value tangent to tract as a minc file and as tract point scalars Distinguished from diffusion ODF case by OTT, not DOTT.  somewhat (!) confusing

need checking to see whether ODF data is an ODF!

ODF_MC option: currently, ODF must already have been calculated everywhere, although could calculate it on the fly as track instead

in FibreTractSet, FibreTract: should change the sort of scalar to be gerneal and have a pair of type-string that can be calculated and string gets written to the vtk file.

in this, completey different code, so should specify.  should allow writing of multiple scalars if desired

currently have MC, but could also have X splits per step (more like Tournier...).  might make more sense: why overdo it on the proximal tracts?  need a fibretract deep copy to restart at each point...

DO: add algorithm ODF_levelset: try it with interplation into precalculated ODF

Q do I free up the ODF_sampler array at the end??

max curvature is set to 90: could in theory have greater but hardcoded as 90 now
OTT would in every case imaginable be fibre ODF

the current minOTT assignments for the maxima are useful just in figuring out whether you were along the maximum or not. However, maxima are guaranteed to lie along sampling directions, so there is always a bias.  Is this a problem??



*/


enum eLabel {unknown, single_curve, fan, crossing, branch, bottleneck, doublecross,triplecross,none} ;

enum eTrackingAlgorithm {PDD, MultiMax, ODF_MC, Cone_MC, Labels} ;

enum eComparison {ge,gt,lt,le} ;

enum eIntegrationType {FACT} ;

struct sPropagationInfo
{
  INDEX3D this_index;
  INDEX3D last_index;
  POINT3D this_point;
  POINT3D last_point;
  VECTOR3D this_vector;
  VECTOR3D last_vector;
  float this_angle;
  float last_angle;
  float this_step;
  float last_step;
  int this_plane; 
  bool start;
  //POINT3D this_voxel_point;
} ;


typedef eLabel Label;
typedef eTrackingAlgorithm TrackingAlgorithm;
typedef eComparison Comparison;
typedef eIntegrationType IntegrationType;
typedef sPropagationInfo PropagationInfo;

class FibreTracking
{
 public:
  FibreTracking();
  ~FibreTracking();
  void SetTrackingAlgorithm(TrackingAlgorithm tracking_algorithm);
  void Track();
  void SetROI(const char *roi_filename);
  void SetROI(ImageData *roi);
  void SetExclusionMask(const char *exclusion_mask_filename);
  void SetMask(ImageData *mask_input, Comparison comparison, float threshold);
  void SetMask(const char *mask_filename, Comparison comparison, float threshold);
  void SetBruteForceOn();
  void SetProcessedDTIData(ProcessedDTIData *processed_DTI_data);
  void SetODFData(ImageData *ODF_data);
  void SetSeedFrequency(int seed_frequency);
  void SetCurvatureConstraint(const char *curvature_constraint_filename);
  void SetCurvatureConstraint(ImageData *curvature_constraint);
  void SetCurvatureConstraint(float curvature_constraint);
  void SetCurvatureConstraintAngle(float curvature_constraint_angle);
  void SetMaxTractLength(float max_tract_length);
  float GetMaxTractLength(); /*IL*/
  FibreTractSet *GetFibreTractSet();
  void SetIntegrationType(IntegrationType integration_type);
  void BinaryConnectivityMapOn();
  //void MinDOTTMapOn();
  void CumMinOTTMapOn();
  void SetPruneOn();
  void SetMCIterations(int iterations);
  void SetLabelsOn(char *labels_filename);
  //void SetLabelsOn();
  ImageData *GetBinaryTrackingMask();
  float GetMaxAngle(PropagationInfo &propagation_info);
  void SetSigmaTheta();
  //void SetSigmaTheta(ImageData *sigma_theta_image);replacing this with an ODFdata that has sigma theta as a member
  void SetSigmaTheta(ImageData *ODFimage_with_sigma_theta);
  void SetSigmaTheta(float sigma_theta_const);
  void SetReferenceROI(int reference_roi);
  //void ReadROI(const char *roi_txt_filename, char *out_filename, const char *tempdirname, const char *template_filename, bool vox);/*use std::string*/
  void ReadROI(std::string & roi_txt_filename,std::string & out_filename,std::string & tempdirname,std::string & template_filename,bool vox);
  void SetProbLabels(ImageData *prob_labels0, ImageData *prob_labels1, ImageData *prob_labels2, ImageData *prob_labels3, ImageData *prob_labels4, ImageData *prob_labels5);
  VECTOR3D Jitter(VECTOR3D v);

 private:
  TrackingAlgorithm m_tracking_algorithm;
  IntegrationType m_integration_type;
  int m_number_of_rois;
  int m_number_of_exclusion_masks;
  FibreTractSet *m_fibre_tract_set;
  ImageData **m_roi_array;
  ImageData **m_exclusion_mask_array;
  int m_mask_exists;
  ImageData *m_mask;
  ImageData *m_curvature_constraint_image;
  float m_curvature_constraint_constant;
  float m_curvature_constraint_angle;
  bool m_variable_curvature_constraint;
  bool m_curvature_constraint_exists; 
  bool m_curvature_constraint_angle_exists; 
  float m_max_angle;
  bool m_max_tract_length_exists; //IL
  float m_max_tract_length; //IL 
  bool m_max_angle_calculated;
  ProcessedDTIData *m_processed_DTI_data;
  ImageData *m_ODF_data;
  bool m_brute_force;
  float m_seeds_per_voxel; //this means per dim per voxel

  bool m_binary_connectivity_map_on;
  //bool m_min_DOTT_map_on;
  //bool m_cum_min_DOTT_map_on;
  //bool m_min_OTT_map_on;
  bool m_cum_min_OTT_map_on;
  IMAGE_GEOM m_image_geom;
  Directions **m_ODF_sampler;
  bool m_prune;
  int m_MC_iterations;
  bool m_ODF_sampler_allocated;
  bool m_labels_on;
  ImageData *m_labels;
  bool m_labels_exist;
  int m_it;  
  float m_sigma_theta_const;
  ImageData *m_sigma_theta_image;  
  float m_OTT;
  int m_reference_roi;
  bool m_cone;
  bool StopTracking(FibreTract *fibre_tract, PropagationInfo &propagation_info);/*IL change to fibre_tract input*/
 // bool StopTracking(int putative_point_number, PropagationInfo &propagation_info);
  bool m_big_voxels;
  ImageData *m_ProbLabels0, *m_ProbLabels1, *m_ProbLabels2, *m_ProbLabels3, *m_ProbLabels4, *m_ProbLabels5;
  float m_scale;
  Label m_label;
  float m_index_max;
  ImageData *m_fullpdf;
  ImageData *m_polarity_vector;

  bool Propagate(FibreTractSet *fibre_tract_set, PropagationInfo &propagation_info);
  void RunStreamline();
  VECTOR3D GetNextVectorFromODF(int x, int y, int z);
  VECTOR3D GetNextVector(int x, int y, int z, bool swap, PropagationInfo propagation_info);
  VECTOR3D GetNextVectorFromMultiMax(int x, int y, int z, PropagationInfo propagation_info);
  VECTOR3D GetNextVectorForLabelTracking(int x, int y, int z, PropagationInfo propagation_info);
  VECTOR3D GetNextVectorFromMultiMaxCone(int x, int y, int z, PropagationInfo propagation_info);
  VECTOR3D GetNextVectorFromPDDCone(int x, int y, int z);
  Label GetLabel(int x, int y, int z, PropagationInfo propagation_info);  
  void PruneBranchedSets(FibreTractSet *tract_set_array[2], ImageData *roi, bool exclude);
  //void CopySetsToFinalSet(FibreTractSet *tract_set_array[2], FibreTractSet *(*final_set), bool prune_same_indices);
  void CopySetsToFinalSet(FibreTractSet *tract_set_array[2], bool prune_same_indices);
  SIGMA_THETA GetSigmaTheta(int x, int y, int z);
  void AssignCumulativeMinOTTs(FibreTractSet *tract_set_array[2], ImageData *roi_ref);
  bool CheckPolarity(VECTOR3D vector, PropagationInfo propagation_info);
  INVERSE *m_point_inverse;
  INVERSE *m_vector_inverse;

  bool m_fanning;
};
#endif
