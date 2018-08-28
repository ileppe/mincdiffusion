/************************************************************************

   File: FibreTracking.c++

   Author: Jennifer Campbell
   Created:
   Revisions:

   Description:

   Copyright (c) 2009 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.

***********************************************************************/
#include <float.h>
#include "FibreTracking.h"
#include "otherlibs_include.h"
#include "FibreTractSet.h"
#include "ImageData.h"
#include "misc.h"
#include "ProcessedDTIData.h"
#include "Directions.h"
#include "DiffusionODFCalculator.h"
#include "DisplayWindow.h"
#include "FibreTract.h"
#include <assert.h>
#include <fstream>

FibreTracking::FibreTracking()
{

  m_seeds_per_voxel=1;
  m_brute_force=false;
  //default is 90 degree which also is what we will get if radius of curvature is set to 0.
  m_curvature_constraint_angle=90;
  m_variable_curvature_constraint=false;
  m_curvature_constraint_angle_exists=true;

  m_curvature_constraint_exists=false;
  m_variable_curvature_constraint=false;
  m_mask_exists=false;
  m_tracking_algorithm=PDD;
  m_number_of_rois=0;
  m_roi_array=(ImageData **) malloc(1*sizeof(ImageData *));
  m_exclusion_mask_array=(ImageData **) malloc(1*sizeof(ImageData *));
  m_integration_type=FACT;
  m_binary_connectivity_map_on=false;
  //m_min_DOTT_map_on=false; //all of the DOTT stuff from here on should be removed
  //m_cum_min_DOTT_map_on=false;
  //m_min_OTT_map_on=false;
  m_cum_min_OTT_map_on=false;
  m_prune=false;//this isn't default because we don't want to prune if calculating cumulative minOTT.  (we also don't want to save tracts!!  and shouldn't...this is temporary
  m_labels_exist=false;
  m_number_of_exclusion_masks=0;

  m_fibre_tract_set=new FibreTractSet();

  m_max_angle_calculated=false;
  m_MC_iterations=1; //this gets upped if algorithm requires more or user asks for more

  //initialize random number generator (for MC algorithms only: could put this elsewhere):

  int seed;
  seed = (int) time(NULL);
  srand48(seed);

  m_sigma_theta_const=0.0;
  m_sigma_theta_image=NULL;
  m_reference_roi=0;
  m_big_voxels=false;
  //stays so unless in ProbLabels, where it should always be updated:
  m_label=single_curve;
  m_index_max=FLT_MAX;
  m_ProbLabels0=NULL;
  m_max_tract_length_exists=false;
  //m_image_geom.vector_inverse_xfm.calculated=false;
  //m_image_geom.point_inverse_xfm.calculated=false;

  //not nec. - done in GetImageGeom() the first time: 
  m_image_geom.vector_inverse_xfm=NULL;
  m_image_geom.point_inverse_xfm=NULL;


  //m_point_inverse=new INVERSE();
  //m_vector_inverse=new INVERSE();

  m_fanning=false;

}


FibreTracking::~FibreTracking()
{
  int i;

  if (m_mask_exists) {
    m_mask->Release();
  }
  //free roi_array: not the roi ImageData pointers themselves.  although some of them should be freed if those files were not assigned to ImageDatas elsewhere!! (should have taken care of that at another point?)
  for (i=0;i<m_number_of_rois;i++) {
    m_roi_array[i]->Release();

  }

  free(m_roi_array);

  for (i=0;i<m_number_of_exclusion_masks;i++) {
    m_exclusion_mask_array[i]->Release();

  }

  free(m_exclusion_mask_array);

  if (m_ODF_sampler_allocated) {
    for (i=0;i<m_image_geom.x_size*m_image_geom.y_size*m_image_geom.z_size;i++) {
      m_ODF_sampler[i]->Release();
    }

    free(m_ODF_sampler);

  }

  if (m_labels_exist) {
    m_labels->Release();
  }

  //the image option is for PDD/eigenvector
  if (m_sigma_theta_image!=NULL) {
    m_sigma_theta_image->Release();
  }

  //DO: for mostlikely, need to Release m_ProbLabels (6 of them)
  //DO: clean up inverses in m_image_geom structure
}

void FibreTracking::SetTrackingAlgorithm(TrackingAlgorithm tracking_algorithm)
{
  m_tracking_algorithm=tracking_algorithm;

  if (tracking_algorithm==ODF_MC || tracking_algorithm==Labels) {

    m_binary_connectivity_map_on=true;
    m_MC_iterations=10000;

  }

  if (tracking_algorithm==Labels && m_fanning)
    {
      m_fullpdf=new ImageData("full_pdf.mnc");
      m_polarity_vector=new ImageData("polarity.mnc");
    }
}




void FibreTracking::Track()
{
  int i,j;
  FibreTractSet *fibre_tract_set_tmp;
  int x,y,z;
  float roi_value;

  if (m_number_of_rois<0 && !m_brute_force) {
    Error("must give at least one ROI unless call SetBruteForceOn()");
  }


  //VF:debug
  /*
  char _tmp[2000];
  for(int i=0;i<m_number_of_rois;i++)
  {
    sprintf(_tmp,"Roi_%d.mnc",i);
    m_roi_array[i]->WriteMincFile(_tmp,true);
  }*/


  switch (m_tracking_algorithm) {
  case PDD:

    m_fibre_tract_set->SetImageDataTemplate(m_processed_DTI_data->Gete1x());


    //if (m_min_DOTT_map_on)
    //{
    //  m_fibre_tract_set->InitMinDOTTMap();
    //}
    //if (m_cum_min_DOTT_map_on)
    //{
    //  m_fibre_tract_set->InitCumMinDOTTMap();
    //}




    break;
  case ODF_MC:

    //DO:can do this with m_image_geom
    m_fibre_tract_set->SetImageDataParameters(m_ODF_data->GetXSize(),m_ODF_data->GetYSize(),m_ODF_data->GetZSize(),m_ODF_data->GetXStart(),m_ODF_data->GetYStart(),m_ODF_data->GetZStart(),m_ODF_data->GetXStep(),m_ODF_data->GetYStep(),m_ODF_data->GetZStep(),m_ODF_data->GetXDirectionCosine(),m_ODF_data->GetYDirectionCosine(),m_ODF_data->GetZDirectionCosine());



    break;
  case Labels:


   
    m_fibre_tract_set->SetImageDataParameters(m_image_geom);



    break;
  case MultiMax:

    m_fibre_tract_set->SetImageDataParameters(m_ODF_data->GetXSize(),m_ODF_data->GetYSize(),m_ODF_data->GetZSize(),m_ODF_data->GetXStart(),m_ODF_data->GetYStart(),m_ODF_data->GetZStart(),m_ODF_data->GetXStep(),m_ODF_data->GetYStep(),m_ODF_data->GetZStep(),m_ODF_data->GetXDirectionCosine(),m_ODF_data->GetYDirectionCosine(),m_ODF_data->GetZDirectionCosine());



    //for this option, we calculate the maxima if they were not input: loop through all voxels and ImageData will calc if not there already (else GetODFMaxima exits)

    //why do we do this?  could just run when needed

    /*
       for (x=0;x<m_ODF_data->GetXSize();x++)
    {
    for (y=0;y<m_ODF_data->GetYSize();y++)
      {
        for (z=0;z<m_ODF_data->GetZSize();z++)
    {
     if (m_mask->GetValue(x,y,z)>FLT_EPSILON)
       {
         m_ODF_data->GetODFMaxima(x,y,z)->Release();
       }

    }
      }
    }
    */



    break;

  }

  if (m_binary_connectivity_map_on) {
    m_fibre_tract_set->InitBinaryConnectivityMap();
  }

  if (m_cum_min_OTT_map_on) {
    m_fibre_tract_set->InitCumMinOTTMap();
  }



  RunStreamline();

  //setting value in ROI to 2.0*the max in the index volume.  this is arbitrary, but will be the highest and will be able to view with reasonable windowing
  //also settting voxels with inf to this value; that happened if sigma_theta was 0 with higher confidence than 0.2 (all hardcoded in various places)

  if (m_cum_min_OTT_map_on) {

    roi_value=2.0*m_fibre_tract_set->GetCumMinOTTMap()->GetMaximumValue(FLT_MAX);

    for (x=0;x<m_fibre_tract_set->GetCumMinOTTMap()->GetXSize();x++) {
      for (y=0;y<m_fibre_tract_set->GetCumMinOTTMap()->GetYSize();y++) {
        for (z=0;z<m_fibre_tract_set->GetCumMinOTTMap()->GetZSize();z++) {
          if (m_fibre_tract_set->GetCumMinOTTMap()->GetValue(x,y,z)!=0) {
            //cout << "index was " << m_fibre_tract_set->GetCumMinOTTMap()->GetValue(x,y,z) << endl;
          }

          if (m_fibre_tract_set->GetCumMinOTTMap()->GetValue(x,y,z)==-10.0 || m_fibre_tract_set->GetCumMinOTTMap()->GetValue(x,y,z)==FLT_MAX) {

            m_fibre_tract_set->GetCumMinOTTMap()->SetValue(x,y,z,roi_value);
          }
        }
      }
    }
  }



}


void FibreTracking::SetROI(const char *roi_filename)
{
  //open and assign to roi_array

  m_roi_array[m_number_of_rois]=new ImageData(roi_filename);
  m_number_of_rois+=1;

  //realloc the roi_array so it can hold the next roi
  m_roi_array=(ImageData **) realloc(m_roi_array,(m_number_of_rois+1)*sizeof(ImageData *));


}

void FibreTracking::SetROI(ImageData *roi)
{
  m_roi_array[m_number_of_rois]=roi;
  m_number_of_rois+=1;

  //realloc the roi_array so it can hold the next roi
  m_roi_array=(ImageData **) realloc(m_roi_array,(m_number_of_rois+1)*sizeof(ImageData *));

}

/*Add this function to be able to read points from input text file
Use a template to create the ROI   ilana*/
void FibreTracking::ReadROI(std::string & roi_txt_filename,std::string & out_filename,std::string & tempdirname,std::string & template_filename,bool vox)
//ReadROI(const char *roi_txt_filename, char *out_filename, const char *tempdirname,const char *template_filename, bool vox)
{
  ImageData *temp;
  ImageData *roi;
  char *tempfilename=(char *) malloc(200*sizeof(char));
  char in[30];
  double x,y,z;

  //sprintf(tempfilename,"%s/tmp%f.mnc",tempdirname, (double) rand());
  sprintf(tempfilename,"%s/tmp%f.mnc",tempdirname.c_str(), (double) rand());

  /*create roi using template*/
  //temp = new ImageData(template_filename);
  temp = new ImageData(template_filename.c_str());

  roi = new ImageData(Image3D_vanilla,temp->GetXSize(),temp->GetYSize(),temp->GetZSize(),temp->GetXStart(),temp->GetYStart(),temp->GetZStart(),temp->GetXStep(),temp->GetYStep(),temp->GetZStep(),temp->GetXDirectionCosine(),temp->GetYDirectionCosine(),temp->GetZDirectionCosine());

  roi -> SetAllVoxels(0);

  /*read input file string by string*/
  //ifstream pFile(roi_txt_filename);
  ifstream pFile(roi_txt_filename.c_str());
  pFile >> in;

  while (pFile.good()) {

    x=atof(in);
    pFile >> in;
    y=atof(in);
    pFile >> in;
    z=atof(in);
    pFile >> in;


    if (!vox) {/*have to make world to voxel conversion if the roi was inputted that way*/
      VECTOR3D voxel;
      voxel = temp -> GetVoxelFromWorldCoordinates(x, y, z);
      /*have to round to integer*/
      roi -> SetValue((int)floor(voxel.x+0.5),(int)floor(voxel.y+0.5),(int)floor(voxel.z+0.5),1); 
      //cout <<"World values: "<<x<<" "<<y<<" "<<z<<endl;
      //cout <<"Voxel values: "<<voxel.z<<" "<<voxel.y<<" "<<voxel.x<<endl;
    } else {
      roi -> SetValue((int)floor(x+0.5),(int)floor(y+0.5),(int)floor(z+0.5),1);
      //cout <<"World values: "<<x<<" "<<y<<" "<<z<<endl;
    }
  }

  /*create temp minc roi file*/
  roi ->WriteMincFile(tempfilename,true); //test to see roi

  //strcpy(out_filename,tempfilename);
  out_filename=tempfilename;
  free(tempfilename);
}


void FibreTracking::SetExclusionMask(const char *exclusion_mask_filename)
{
  //open and assign to exclusion_mask_array

  m_exclusion_mask_array[m_number_of_exclusion_masks]=new ImageData(exclusion_mask_filename);
  m_number_of_exclusion_masks+=1;

  //realloc the exclusion_mask_array so it can hold the next exclusion_mask
  m_exclusion_mask_array=(ImageData **) realloc(m_exclusion_mask_array,(m_number_of_exclusion_masks+1)*sizeof(ImageData *));
}


void FibreTracking::SetMask(ImageData *mask_input, Comparison comparison, float threshold)
{
  //this assumes mask has same start,step,size as all other data files to be used in this FibreTracking.
  int x,y,z;

  if (!m_mask_exists) {

    //create mask ImageData
    m_mask=new ImageData(Image3D_vanilla,
                         mask_input->GetXSize(), mask_input->GetYSize(), mask_input->GetZSize(),
                         mask_input->GetXStart(),mask_input->GetYStart(),mask_input->GetZStart(),
                         mask_input->GetXStep(), mask_input->GetYStep(), mask_input->GetZStep(),
                         mask_input->GetXDirectionCosine(),mask_input->GetYDirectionCosine(),mask_input->GetZDirectionCosine());
  }

  //loop through each voxel and do appropriate comparison
  for (x=0;x<m_mask->GetXSize();x++) {
    for (y=0;y<m_mask->GetYSize();y++) {
      for (z=0;z<m_mask->GetZSize();z++) {
        switch (comparison) {
        case ge:

          if (!m_mask_exists) {
            if (mask_input->GetValue(x,y,z)>=threshold) {
              m_mask->SetValue(x,y,z,1);
            } else {
              m_mask->SetValue(x,y,z,0);
            }
          } else { //m_mask_exists
            if (mask_input->GetValue(x,y,z)<threshold && m_mask->GetValue(x,y,z)>FLT_EPSILON) {
              m_mask->SetValue(x,y,z,0);
            }

            //alt to above is to set all to zero without checking
          }

          break;

        case gt:

          if (!m_mask_exists) {
            if (mask_input->GetValue(x,y,z)>threshold) {
              m_mask->SetValue(x,y,z,1);
            } else {
              m_mask->SetValue(x,y,z,0);
            }
          } else { //m_mask_exists
            if (mask_input->GetValue(x,y,z)<=threshold && m_mask->GetValue(x,y,z)>FLT_EPSILON) {
              m_mask->SetValue(x,y,z,0);
            }

            //alt to above is to set all to zero without checking
          }

          break;

        case lt:

          if (!m_mask_exists) {
            if (mask_input->GetValue(x,y,z)<threshold) {
              m_mask->SetValue(x,y,z,1);
            } else {
              m_mask->SetValue(x,y,z,0);
            }
          } else { //m_mask_exists
            if (mask_input->GetValue(x,y,z)>=threshold && m_mask->GetValue(x,y,z)>FLT_EPSILON) {
              m_mask->SetValue(x,y,z,0);
            }

            //alt to above is to set all to zero without checking
          }

          break;

        case le:

          if (!m_mask_exists) {
            if (mask_input->GetValue(x,y,z)<=threshold) {
              m_mask->SetValue(x,y,z,1);
            } else {
              m_mask->SetValue(x,y,z,0);
            }
          } else { //m_mask_exists
            if (mask_input->GetValue(x,y,z)>threshold && m_mask->GetValue(x,y,z)>FLT_EPSILON) {
              m_mask->SetValue(x,y,z,0);
            }

            //alt to above is to set all to zero without checking
          }

          break;

        }
      }
    }
  }


  m_mask_exists=true;

}

void FibreTracking::SetMask(const char* mask_filename, Comparison comparison, float threshold)
{
  ImageData *mask_input;

  mask_input=new ImageData(mask_filename);

  this->SetMask(mask_input,comparison, threshold);


}

void FibreTracking::SetBruteForceOn()
{
  m_brute_force=true;

}

void FibreTracking::SetProcessedDTIData(ProcessedDTIData *processed_DTI_data)
{
  m_processed_DTI_data=processed_DTI_data;

  m_image_geom.x_size=m_processed_DTI_data->GetXSize();
  m_image_geom.y_size=m_processed_DTI_data->GetYSize();
  m_image_geom.z_size=m_processed_DTI_data->GetZSize();
  m_image_geom.x_start=m_processed_DTI_data->GetXStart();
  m_image_geom.y_start=m_processed_DTI_data->GetYStart();
  m_image_geom.z_start=m_processed_DTI_data->GetZStart();
  m_image_geom.x_step=m_processed_DTI_data->GetXStep();
  m_image_geom.y_step=m_processed_DTI_data->GetYStep();
  m_image_geom.z_step=m_processed_DTI_data->GetZStep();
  m_image_geom.x_direction_cosine=m_processed_DTI_data->GetXDirectionCosine();
  m_image_geom.y_direction_cosine=m_processed_DTI_data->GetYDirectionCosine();
  m_image_geom.z_direction_cosine=m_processed_DTI_data->GetZDirectionCosine();

}


void FibreTracking::SetODFData(ImageData *ODF_data)
{
  m_ODF_data=ODF_data;

  m_image_geom.x_size=m_ODF_data->GetXSize();
  m_image_geom.y_size=m_ODF_data->GetYSize();
  m_image_geom.z_size=m_ODF_data->GetZSize();
  m_image_geom.x_start=m_ODF_data->GetXStart();
  m_image_geom.y_start=m_ODF_data->GetYStart();
  m_image_geom.z_start=m_ODF_data->GetZStart();
  m_image_geom.x_step=m_ODF_data->GetXStep();
  m_image_geom.y_step=m_ODF_data->GetYStep();
  m_image_geom.z_step=m_ODF_data->GetZStep();
  m_image_geom.x_direction_cosine=m_ODF_data->GetXDirectionCosine();
  m_image_geom.y_direction_cosine=m_ODF_data->GetYDirectionCosine();
  m_image_geom.z_direction_cosine=m_ODF_data->GetZDirectionCosine();

  m_ODF_sampler=(Directions **) malloc(m_image_geom.x_size*m_image_geom.y_size*m_image_geom.z_size*sizeof(Directions *));
  m_ODF_sampler_allocated=true;



}


void FibreTracking::SetSeedFrequency(int seeds_per_voxel)
{
  m_seeds_per_voxel=(float) seeds_per_voxel;

}

void FibreTracking::SetCurvatureConstraint(const char *curvature_constraint_filename)
{
  m_curvature_constraint_image = new ImageData(curvature_constraint_filename);
  m_variable_curvature_constraint=true;
  m_curvature_constraint_exists=true;

}



void FibreTracking::SetCurvatureConstraint(ImageData *curvature_constraint)
{
  m_curvature_constraint_image=curvature_constraint;
  m_variable_curvature_constraint=true;
  m_curvature_constraint_exists=true;


}

void FibreTracking::SetCurvatureConstraint(float curvature_constraint)
{


  m_curvature_constraint_constant=curvature_constraint;
  m_variable_curvature_constraint=false;
  m_curvature_constraint_exists=true;

}

void FibreTracking::SetCurvatureConstraintAngle(float curvature_constraint_angle)
{

  if (curvature_constraint_angle>90) {
    m_big_voxels=true;
  }

  m_curvature_constraint_angle=curvature_constraint_angle;

  m_variable_curvature_constraint=false;
  m_curvature_constraint_angle_exists=true;
}

void FibreTracking::SetMaxTractLength(float max_tract_length) /*IL*/
{
  m_max_tract_length=max_tract_length;
  m_max_tract_length_exists=true;
}

float FibreTracking::GetMaxTractLength()  /*IL*/
{
  return m_max_tract_length;
}

FibreTractSet *FibreTracking::GetFibreTractSet()
{
  return m_fibre_tract_set;

}

void FibreTracking::SetIntegrationType(IntegrationType integration_type)
{
  m_integration_type=integration_type;

}

void FibreTracking::BinaryConnectivityMapOn()
{
  m_binary_connectivity_map_on=true;

}

//void FibreTracking::MinDOTTMapOn()
//{
// m_min_DOTT_map_on=true;

//}

void FibreTracking::CumMinOTTMapOn()
{
  m_cum_min_OTT_map_on=true;
}


float FibreTracking::GetMaxAngle(PropagationInfo &propagation_info)
{
  if (m_max_angle_calculated && !m_variable_curvature_constraint) {
    return m_max_angle;
  }

  if (m_curvature_constraint_angle_exists) {
    m_max_angle=(M_PI/180.0)*m_curvature_constraint_angle;
  } else if (m_variable_curvature_constraint) {
    if (m_image_geom.x_step/2.0>m_curvature_constraint_image->GetValue(propagation_info.this_index.x,propagation_info.this_index.y,propagation_info.this_index.z)) {
      m_max_angle=M_PI/2.0;

    } else {

      m_max_angle=2*asin(float_abs(m_image_geom.x_step)/(2.0*m_curvature_constraint_image->GetValue(propagation_info.this_index.x,propagation_info.this_index.y,propagation_info.this_index.z)));


    }
    //note this assumes isotropic voxel size
  } else {
    if (float_abs(m_image_geom.x_step)/2.0>m_curvature_constraint_constant) {
      m_max_angle=M_PI/2.0;

    } else {

      m_max_angle=2*asin(float_abs(m_image_geom.x_step)/(2.0*m_curvature_constraint_constant));
    }


  }

  if (m_max_angle>M_PI/2.0) {
    m_max_angle=M_PI/2.0;

  }

  m_max_angle_calculated=true;

  return m_max_angle;
}

bool FibreTracking::StopTracking(FibreTract *fibre_tract, PropagationInfo &propagation_info)
{
  //this function returns true if the voxel we're entering next is out of bounds or the last angle was not allowed (both mean the point just calculated should not be assigned).
  // this function will also return true if the tract length exceeds a facultative user defined maximum length in mm  IL
  bool result=true;
  float max_angle;
  int putative_point_number = fibre_tract->GetNumberOfPoints()+1;

  //cout << "StopTracking\n";

  //check if in volume:

  if (propagation_info.this_index.x>=0 && propagation_info.this_index.x<m_image_geom.x_size && propagation_info.this_index.y>=0 && propagation_info.this_index.y<m_image_geom.y_size && propagation_info.this_index.z>=0 && propagation_info.this_index.z<m_image_geom.z_size) {
    result=false;


  }
  //debug:
  else
  {
    //cout << "outside image\n";
  }
  //end debug


  //check mask:

  if (m_mask_exists && result==false) {

    if (m_mask->GetValue(propagation_info.this_index.x,propagation_info.this_index.y,propagation_info.this_index.z)>FLT_EPSILON) {
      result=false;
    }

    else {
      //debug
      //cout << "outside mask\n";
      //end debug
      result=true;
    }


  }


  //check curvature constraint

  //not doing radius of curvature and consideration of last angle here for now:
  if (0)
    {
  //check the last angle:
  if (putative_point_number>2 && result==false) {

    propagation_info.last_angle=propagation_info.this_angle;

    //propagation_info.this_angle=acos(float_abs(propagation_info.this_vector.x*propagation_info.last_vector.x+propagation_info.this_vector.y*propagation_info.last_vector.y+propagation_info.this_vector.z*propagation_info.last_vector.z));

    //cout << "Stop: this vector " <<  propagation_info.this_vector.x << " " << propagation_info.this_vector.y << " " << propagation_info.this_vector.z << endl;

    if (propagation_info.this_vector.x*propagation_info.last_vector.x+propagation_info.this_vector.y*propagation_info.last_vector.y+propagation_info.this_vector.z*propagation_info.last_vector.z>1.0+FLT_EPSILON) { 
      propagation_info.this_angle=0;
      //debug
      //cout << "length above 1 - " << propagation_info.this_vector.x*propagation_info.last_vector.x+propagation_info.this_vector.y*propagation_info.last_vector.y+propagation_info.this_vector.z*propagation_info.last_vector.z<< " this " << propagation_info.this_vector.x*propagation_info.this_vector.x+propagation_info.this_vector.y*propagation_info.this_vector.y+propagation_info.this_vector.z*propagation_info.this_vector.z << " last " << propagation_info.last_vector.x*propagation_info.last_vector.x+propagation_info.last_vector.y*propagation_info.last_vector.y+propagation_info.last_vector.z*propagation_info.last_vector.z << endl;

      //end debug
    } 
  }
    } //end of radius of curvature, etc.
  //else 
  if (putative_point_number>2 && result==false) 
    {
      propagation_info.this_angle=acos(propagation_info.this_vector.x*propagation_info.last_vector.x+propagation_info.this_vector.y*propagation_info.last_vector.y+propagation_info.this_vector.z*propagation_info.last_vector.z); 
      //returns nan if same vector!
      if (isnan(propagation_info.this_angle))
	{ 
	    propagation_info.this_angle=0;
	}
    max_angle=this->GetMaxAngle(propagation_info);

    if (max_angle>propagation_info.this_angle) {
      result=false;
    } else {
      //debug
      //cout << "this angle: " << (180.0/M_PI)*propagation_info.this_angle << endl;
      //cout << "max angle " << (180.0/M_PI)*max_angle << " exceeded\n";
      //cout << "this vector: " <<  propagation_info.this_vector.x << " " << propagation_info.this_vector.y << " " << propagation_info.this_vector.z << endl;
      //cout << "last vector: " <<  propagation_info.last_vector.x << " " << propagation_info.last_vector.y << " " << propagation_info.last_vector.z << endl;
      //end debug
      result=true;
    }
    }

    //}

  //check the total angle per arc length in the last two angles:
    if (0)
      {
  //if (putative_point_number>3 && result==false && !m_curvature_constraint_angle_exists)
  if (putative_point_number>3 && result==false) {
    //cout << "this angle " <<  propagation_info.this_angle << " last angle " <<  propagation_info.last_angle << " this step " << propagation_info.this_step << " last step " << propagation_info.last_step << endl;

    //I'm afraid this won't work if the voxels are "big" relative to the anatomy.  hence, allowing user to shut it off

    if (!m_big_voxels) {
      if (((propagation_info.this_angle+propagation_info.last_angle)/(propagation_info.this_step+propagation_info.last_step))<max_angle/float_abs(m_image_geom.x_step)) {
        result=false;


      }

      else {
        //debug:
        //cout << "max angle exceeded 2\n";
        ////result=false; //test
        //end debug

        result=true;
      }
    }

  
  }
      }

  /*Check if length exceeds user defined maximum*/ /*IL*/
  if (m_max_tract_length_exists && result==false) {
    if (fibre_tract->GetTractLength() < this->GetMaxTractLength()) {
      result = false; /*it's shorter than the max*/
    } else {
      //debug:
      //cout << "max length "<<  this->GetMaxTractLength() <<" exceeded: "<< fibre_tract->GetTractLength()<<"\n"<<endl;
      //end debug
      result = true;
    }
  }


  return result;


}


bool FibreTracking::Propagate(FibreTractSet *fibre_tract_set, PropagationInfo &propagation_info)

{
  bool result=true;
  float a;
  float x,y,z;
  bool intercept_found=false;
  Directions *directions;
  DiffusionODFCalculator *DODFCalculator;

  int number_of_directions;



  PropagationInfo *propagation_info_copy;
  bool propagation_info_copy_allocated=false;

  FibreTract *fibre_tract=fibre_tract_set->GetTract(fibre_tract_set->GetNumberOfTracts()-1);

  FibreTract *new_tract;
  FibreTract *prev_tract;


  VECTOR3D voxel_vector;
  POINT3D voxelstar_point;
  POINT3D voxel_point;

  int i;

  

  //note that need to check vector in m_branch case, else GetNextVector does it
 



  if (fibre_tract->GetNumberOfPoints()>1) //this was commented but I don't know why...it is needed..

  {

 
    propagation_info.last_step=propagation_info.this_step;


    

    if (propagation_info.this_vector.x*propagation_info.last_vector.x+propagation_info.this_vector.y*propagation_info.last_vector.y+propagation_info.this_vector.z*propagation_info.last_vector.z<0)
    {

      propagation_info.this_vector.x=-propagation_info.this_vector.x;
      propagation_info.this_vector.y=-propagation_info.this_vector.y;
      propagation_info.this_vector.z=-propagation_info.this_vector.z;

    }
    
  }


  //cout << "vector " << propagation_info.this_vector.x << " " << propagation_info.this_vector.y << " " << propagation_info.this_vector.z << endl;


  propagation_info.last_point=propagation_info.this_point;

  propagation_info.last_index=propagation_info.this_index;




  //set the OTTs here: 
  //JC removed doing this if following the maxima only.  could, but not.  if following the maxima and the cone around them, OTT will be assigned using gaussian with sigma_theta.

  if (m_cum_min_OTT_map_on && m_tracking_algorithm==ODF_MC) { 
    directions=new Directions(1);
    directions->SetXComponent(0,propagation_info.this_vector.x);
    directions->SetYComponent(0,propagation_info.this_vector.y);
    directions->SetZComponent(0,propagation_info.this_vector.z);

    //make an interpolator and find ODF value in the dir of this vector:


    S2Interpolator *interpolator=new S2Interpolator(m_ODF_data,directions,1.5*sqrt(2*M_PI/m_ODF_data->GetNSize()), 1.5*sqrt(2*M_PI/m_ODF_data->GetNSize()), Gaussian);

    fibre_tract->SetOTT(interpolator->GetInterpolatedValue(0, propagation_info.last_index.x,propagation_info.last_index.y,propagation_info.last_index.z),fibre_tract->GetNumberOfPoints()-1);


    directions->Release();
    delete(interpolator);


  

  }

  else if (m_cum_min_OTT_map_on && m_cone) {
   
    if (m_OTT>=FLT_MAX) { //nan - DO: figure out where that is coming from.  one case is the sigma_theta==0 case: OK.
      m_OTT=m_index_max;//which is FLT_MAX right now so this isn't necessary
    }

    fibre_tract->SetOTT(m_OTT,fibre_tract->GetNumberOfPoints()-1);
  }




  switch (m_integration_type) {
  case FACT:


    voxel_vector=WorldToVoxel(propagation_info.this_vector,&m_image_geom);

    
    //voxelstar coords; we use an anisotropic voxel, but sides are voxel x,y,z for the FACT integration

    
    //voxel_point=WorldToVoxel(propagation_info.this_point,*m_image_geom);
    
    voxel_point=WorldToVoxel(propagation_info.this_point,&m_image_geom);
    


    voxelstar_point.x=voxel_point.x-(float)propagation_info.this_index.x+0.5;
    voxelstar_point.y=voxel_point.y-(float)propagation_info.this_index.y+0.5;
    voxelstar_point.z=voxel_point.z-(float)propagation_info.this_index.z+0.5;
    
    if (voxelstar_point.x<0) // || voxelstar_point.x>float_abs(m_image_geom.x_step))
      {
	voxelstar_point.x=0;
      }

    
    if (voxelstar_point.y<0) // || voxelstar_point.y>float_abs(m_image_geom.y_step))
      {
	voxelstar_point.y=0;
      }

    
    if (voxelstar_point.z<0) // || voxelstar_point.z>float_abs(m_image_geom.z_step))
      {
	voxelstar_point.z=0;
      }
    
    
    //POINT3D test=VoxelToWorld(voxel_point,m_image_geom);
    //VECTOR3D testv=VoxelToWorld(voxel_vector,m_image_geom);
    //cout << "world point " << propagation_info.this_point.x << " " << propagation_info.this_point.y << " " << propagation_info.this_point.z << "\nvoxel_point " << voxel_point.x << " " << voxel_point.y << " " << voxel_point.z << "\nindex " << propagation_info.this_index.x << " " << propagation_info.this_index.y << " " << propagation_info.this_index.z << "\nvoxelstar_point " << voxelstar_point.x << " " << voxelstar_point.y << " " << voxelstar_point.z <<  endl;
    //cout << "world vector " << propagation_info.this_vector.x << " " << propagation_info.this_vector.y << " " << propagation_info.this_vector.z << "\nvoxel vector " << voxel_vector.x << " " << voxel_vector.y << " " << voxel_vector.z <<  endl;

    //"\ntest back to world " << test.x << " " << test.y << " " << test.z << 
    //"\nback to world " << testv.x << " " << testv.y << " " << testv.z <<
    //cout << "plane " << propagation_info.this_plane << endl;
    

    //note: I don't check for exactly zero vector components (divide by zero)

    // z -



    if (propagation_info.this_plane!=1 && voxel_vector.z<0)
      {
	//cout << "z -" << endl;




      a=-voxelstar_point.z/voxel_vector.z;
      
      
      x=voxelstar_point.x+a*voxel_vector.x;
      y=voxelstar_point.y+a*voxel_vector.y;
  

      if (x<=1.0 && x>=0.0 && y<=1.0 && y>=0.0 ) {

	propagation_info.this_point=VoxelToWorld(x+(float)propagation_info.this_index.x-0.5,y+(float)propagation_info.this_index.y-0.5,(float)propagation_info.this_index.z-0.5,m_image_geom);
      
	
        propagation_info.this_plane=6;

        propagation_info.this_index.z-=1;

        //cout << "yes: z-" << endl;
        intercept_found=true;


      }


      }

    // y -

    
    if (intercept_found==false && propagation_info.this_plane!=2 && voxel_vector.y<0) {
      //cout << "y -" << endl;




      
      a=-voxelstar_point.y/voxel_vector.y;

      x=voxelstar_point.x+a*voxel_vector.x;
      z=voxelstar_point.z+a*voxel_vector.z;
 
      if (x<=1.0 && x>=0.0 && z<=1.0 && z>=0.0 )  {
 



	propagation_info.this_point=VoxelToWorld(x+(float)propagation_info.this_index.x-0.5,(float)propagation_info.this_index.y-0.5,z+(float)propagation_info.this_index.z-0.5,m_image_geom);
      
	//cout << "voxel " << x+(float)propagation_info.this_index.x-0.5 << " " << (float)propagation_info.this_index.y-0.5 << " " << z+(float)propagation_info.this_index.z-0.5 << " world " << propagation_info.this_point.x << " " << propagation_info.this_point.y << " " << propagation_info.this_point.z << endl;



      propagation_info.this_plane=4;

 
      propagation_info.this_index.y-=1;

      //cout << "yes: y-" << endl;

      intercept_found=true;



      }


    }    

    // x -



    if (intercept_found==false && propagation_info.this_plane!=3 && voxel_vector.x<0) {
      //cout << "x -" << endl;



      
      a=-voxelstar_point.x/voxel_vector.x;

      //cout << "a " << a << endl;
      
      z=voxelstar_point.z+a*voxel_vector.z;
      y=voxelstar_point.y+a*voxel_vector.y;



      //cout << "z " << z << " y " << y << endl;

      if (z<=1.0 && z>=0.0 && y<=1.0 && y>=0.0 )  {

	propagation_info.this_point=VoxelToWorld((float)propagation_info.this_index.x-0.5,y+(float)propagation_info.this_index.y-0.5,z+(float)propagation_info.this_index.z-0.5,m_image_geom);
      
  

        propagation_info.this_plane=5;


        propagation_info.this_index.x-=1;

        //cout << "yes: x-" << endl;
        intercept_found=true;

      }


    }



    // z +



    if (intercept_found==false && propagation_info.this_plane!=6 && voxel_vector.z>=0) {


      //cout << "z +" << endl;

      a=(1.0-voxelstar_point.z)/voxel_vector.z;



      x=voxelstar_point.x+a*voxel_vector.x;
      y=voxelstar_point.y+a*voxel_vector.y;
  

      if (x<=1.0 && x>=0.0 && y<=1.0 && y>=0.0 ) {

	propagation_info.this_point=VoxelToWorld(x+(float)propagation_info.this_index.x-0.5,y+(float)propagation_info.this_index.y-0.5,(float)propagation_info.this_index.z+0.5,m_image_geom);
        

        propagation_info.this_plane=1;



        propagation_info.this_index.z+=1;

        //cout << "yes: z+" << endl;
        intercept_found=true;

      }


    }

    // y +



    if (intercept_found==false && propagation_info.this_plane!=4 && voxel_vector.y>=0) {
      //cout << "y +" << endl;


     
    
      a=(1.0-voxelstar_point.y)/voxel_vector.y;

      x=voxelstar_point.x+a*voxel_vector.x;
      z=voxelstar_point.z+a*voxel_vector.z;
 
      if (x<=1.0 && x>=0.0 && z<=1.0 && z>=0.0 )  {
 
	propagation_info.this_point=VoxelToWorld(x+(float)propagation_info.this_index.x-0.5,(float)propagation_info.this_index.y+0.5,z+(float)propagation_info.this_index.z-0.5,m_image_geom);
      

        propagation_info.this_plane=2;


        propagation_info.this_index.y+=1;

        //cout << "yes: y+" << endl;
        intercept_found=true;

      }


    }

    // x +



    if (intercept_found==false && propagation_info.this_plane!=5 && voxel_vector.x>=0) {
      //cout << "x +" << endl;


     
      a=(1.0-voxelstar_point.x)/voxel_vector.x;

      //cout << "a " << a << endl;
      
      z=voxelstar_point.z+a*voxel_vector.z;
      y=voxelstar_point.y+a*voxel_vector.y;

      //cout << "z " << z << " y " << y << endl;

      if (z<=1.0 && z>=0.0 && y<=1.0 && y>=0.0 )  {
	
	propagation_info.this_point=VoxelToWorld((float)propagation_info.this_index.x+0.5,y+(float) propagation_info.this_index.y-0.5,z+(float)propagation_info.this_index.z-0.5,m_image_geom);
      


        propagation_info.this_plane=3;


        propagation_info.this_index.x+=1;

        //cout << "yes: x+" << endl;
        intercept_found=true;


      }


    }




    if (intercept_found==false) { //can happen if 90 degrees exceeded, which should never happen...
      //cout << "intercept not found" << endl;

      //Error("intercept not found");
      result=false;

    }

    break;

  }



 

  if (result==true) {

    propagation_info.this_step=sqrt(pow(propagation_info.this_point.x-propagation_info.last_point.x,2)+pow(propagation_info.this_point.y-propagation_info.last_point.y,2)+pow(propagation_info.this_point.z-propagation_info.last_point.z,2));

    //cout << "Prop: this step: " << propagation_info.this_step << endl;
    //cout << "found intercept\n";

    if (StopTracking(fibre_tract,propagation_info)) {
      //if (fibre_tract->GetNumberOfUniquePoints()<1)
      //{
      //  fibre_tract_set->RemoveTract(fibre_tract_set->GetNumberOfTracts()-1);
      //}
      //cout << "StopTracking\n";
      result=false;
    }



    else {

      
      //cout << propagation_info.this_point.x << " " << propagation_info.this_point.y << " " << propagation_info.this_point.z << endl;


      fibre_tract->AddPoint(propagation_info.this_point.x,propagation_info.this_point.y,propagation_info.this_point.z);
      fibre_tract->SetIndex(propagation_info.this_index.x,propagation_info.this_index.y,propagation_info.this_index.z,fibre_tract->GetNumberOfPoints()-1);



      if (m_binary_connectivity_map_on) {
	//not doing this any more:
        //m_fibre_tract_set->GetBinaryConnectivityMap()->SetValue(propagation_info.last_index.x,propagation_info.last_index.y,propagation_info.last_index.z,1);
        //test: not nec.
        //fibre_tract_set->GetBinaryConnectivityMap()->SetValue(propagation_info.last_index.x,propagation_info.last_index.y,propagation_info.last_index.z,1);
      }



      //one last check: did the tract make a loop?

      //test:
      //if (fibre_tract->GetNumberOfPoints()>1000)
      //{
      //  result=false;
      //}
      if (fibre_tract->CheckForDuplicateIndex(true)) {
        result=false;

      } else {
        result=true;
      }



      //now get vector heading for next point:

     

      propagation_info.last_vector.x=propagation_info.this_vector.x;

      propagation_info.last_vector.y=propagation_info.this_vector.y;

      propagation_info.last_vector.z=propagation_info.this_vector.z;

      //let everybody know we're not on the first anymore:
      propagation_info.start=false;


      //figure out whether to branch or not: the labels algorithm is necessarily iterative to handle the fannings, so will branch that way instead

      if (m_tracking_algorithm==MultiMax && m_labels_exist) {

        if (GetLabel(propagation_info.this_index.x,propagation_info.this_index.y,propagation_info.this_index.z, propagation_info)!=branch && result==true) { //don't branch

          propagation_info.this_vector=GetNextVector(propagation_info.this_index.x,propagation_info.this_index.y,propagation_info.this_index.z,true,propagation_info);




          Propagate(fibre_tract_set,propagation_info);

        } else if (result==true) { //branch


          directions=m_ODF_data->GetODFMaxima(propagation_info.this_index.x,propagation_info.this_index.y,propagation_info.this_index.z);
          directions->Retain();
          number_of_directions=directions->GetNumberOfDirections();
          //cout << number_of_directions << endl;
          //DO: check the origin of this error:

          if (number_of_directions<0) {

            result=false;
          }

          if (number_of_directions==1) {
            //cout << "don't branch\n";
            propagation_info.this_vector=directions->GetDirection(0);
            Propagate(fibre_tract_set,propagation_info);
          } else { //number_of_directions>1
            propagation_info_copy_allocated=true;
            //cout << "branch, " << number_of_directions  << endl;
            //make a copy of the prop info for each direction, and copy the tracts:
            prev_tract=fibre_tract_set->GetTract(fibre_tract_set->GetNumberOfTracts()-1);
            prev_tract->Retain();
            fibre_tract_set->RemoveTract(fibre_tract_set->GetNumberOfTracts()-1);


            //make copy of prop. info here:
            propagation_info_copy=(PropagationInfo *) malloc(number_of_directions*sizeof(PropagationInfo));

            for (i=0;i<number_of_directions && i<number_of_directions;i++) {
              propagation_info_copy[i]=propagation_info;


              //this will create a new tract that points to the previous segment, which does not change:
              new_tract=new FibreTract(prev_tract);



              propagation_info_copy[i].this_vector=directions->GetDirection(i);




              fibre_tract_set->AddTract(new_tract);

              new_tract->Release();

              Propagate(fibre_tract_set,propagation_info_copy[i]);
            }

            prev_tract->Release();

          }

          directions->Release();
        }

      }


      else { //labels aren't on so no chance of branch:
        if (result==true) {
	  
          propagation_info.this_vector=GetNextVector(propagation_info.this_index.x,propagation_info.this_index.y,propagation_info.this_index.z,true,propagation_info);

          if (flt_norm(propagation_info.this_vector)>FLT_EPSILON) {
	    
            Propagate(fibre_tract_set,propagation_info);
          }
	  
        }
      }



    }

  }



  if (propagation_info_copy_allocated) {
    delete propagation_info_copy;
  }

  return result;


}


void FibreTracking::RunStreamline()
{


  int x,y,z;
  FibreTract *fibre_tract;
  int i,j,k;
  int forward;
  int iterations;
  FibreTractSet *fibre_tract_set_tmp;
  bool start_at_this_voxel;
  Directions *directions;
  PropagationInfo propagation_info;
  POINT3D voxel_point;
  POINT3D world_point;
  int start_roi;
  int dir_counter;
  int dir_counter_max=1;
  FibreTractSet *tract_set_array[2];

  bool include;
  int t;

  int tract_counter=0;

  int debug_counter=0;

  int r;



  iterations=m_MC_iterations;

  if (!m_brute_force) {
    start_roi=1;
  } else {
    start_roi=0;
  }





  //initialize the sets.   tracts are saved by default.  scalar maps depend on final set:

  for (forward=0;forward<2;forward++) {
    tract_set_array[forward]=new FibreTractSet();
    tract_set_array[forward]->SetImageDataParameters(m_image_geom);

    //if (m_min_DOTT_map_on)
    //{
    //  tract_set_array[forward]->InitMinDOTTMap();
    //}

    if ((m_cum_min_OTT_map_on && m_brute_force) || (m_cum_min_OTT_map_on && m_reference_roi>0))//isn't used for non brute force case: just use master. but do use for case reference_roi>0
      tract_set_array[forward]->InitCumMinOTTMap();
  }

  double total_number_of_points=m_image_geom.x_size*m_image_geom.y_size*m_image_geom.z_size;

  double last_progress=0.0;
  double progress_step=5.0; //report every 5%

  //std::cout.precision(1);//for progress reporting

  for (x=0;x<m_image_geom.x_size;x++) {
    for (y=0;y<m_image_geom.y_size;y++) {
      for (z=0;z<m_image_geom.z_size;z++) {

        //cout << m_image_geom.x_size << " " << m_image_geom.y_size << " " << m_image_geom.z_size << endl;
        if (m_mask_exists) {
          if (m_mask->GetValue(x,y,z)>FLT_EPSILON) {
            start_at_this_voxel=true;

          } else {
            start_at_this_voxel=false;

          }
        } else {
          start_at_this_voxel=true;

        }


        if (start_at_this_voxel) {

          if (m_brute_force) {
            start_at_this_voxel=true;
          } else if (m_roi_array[0]->GetValue(x,y,z)>FLT_EPSILON) {
            start_at_this_voxel=true;
          } else {
            start_at_this_voxel=false;

          }

        }

        

        if (start_at_this_voxel) {

        debug_counter++;

	//progress has to have the total number of start voxels in start ROI.
        //double progress=100.0*debug_counter/total_number_of_points;

        //if ((progress-last_progress)>=progress_step || last_progress==0.0) {
          //std::cout<<floor(progress)<<"%\t"<<std::flush;
          //last_progress=progress;
        //}

          cout << "\r" << debug_counter;
          //calculate approximate %% of completeon

          //cout << x<<" " << y<< " " << z<< endl;
          //determine number of iterations for MultiMax: number of maxima
          //also if Labels and not ProbLabels, but that is not being maintained
          if (m_tracking_algorithm==MultiMax) {
            directions=m_ODF_data->GetODFMaxima(x,y,z);
            directions->Retain();
            dir_counter_max=directions->GetNumberOfDirections();
          }

          //get the number of starts for Labels tracking: will sample all possible geometries: set directions here
          if (m_ProbLabels0!=NULL) {
            dir_counter_max=1;  //through iteration, will sample everything
          }


          for (i=0;i<m_seeds_per_voxel;i++) {
            for (j=0;j<m_seeds_per_voxel;j++) {
              for (k=0;k<m_seeds_per_voxel;k++) {

                for (m_it=0;m_it<iterations;m_it++) {
                  //cout << "it " << m_it << endl;
                  for (dir_counter=0;dir_counter<dir_counter_max;dir_counter++) {



                    //cout << "it " << it << " dir " <<  dir_counter
                    //make forward and backward tract:
                    for (forward=0;forward<2;forward++) {

                      propagation_info.this_index.x=x;
                      propagation_info.this_index.y=y;
                      propagation_info.this_index.z=z;

                      propagation_info.start=true;

		      /*delete
                      propagation_info.this_point.x=m_image_geom.x_start+m_image_geom.x_step*(-0.5+1.0/(2.0*m_seeds_per_voxel)+i/m_seeds_per_voxel+x);
                      propagation_info.this_point.y=m_image_geom.y_start+m_image_geom.y_step*(-0.5+1.0/(2.0*m_seeds_per_voxel)+j/m_seeds_per_voxel+y);
                      propagation_info.this_point.z=m_image_geom.z_start+m_image_geom.z_step*(-0.5+1.0/(2.0*m_seeds_per_voxel)+k/m_seeds_per_voxel+z);
		      */

		      
		      //0 index is centre of voxel. assigning points from voxel -0.5 to voxel 0.5.
		      propagation_info.this_point=VoxelToWorld(x-0.5+1.0/(2.0*m_seeds_per_voxel)+i/m_seeds_per_voxel, y-0.5+1.0/(2.0*m_seeds_per_voxel)+j/m_seeds_per_voxel,z-0.5+1.0/(2.0*m_seeds_per_voxel)+k/m_seeds_per_voxel,m_image_geom);



                      if (m_tracking_algorithm==MultiMax) {
                        propagation_info.this_vector=directions->GetDirection(dir_counter);
                        propagation_info.start=false;
                      }



                      else { //one vector or probabilistic
                        propagation_info.this_vector=GetNextVector(propagation_info.this_index.x,propagation_info.this_index.y,propagation_info.this_index.z,true,propagation_info);
                      }

                      if (forward==1) {
                        propagation_info.this_vector.x=-propagation_info.this_vector.x;
                        propagation_info.this_vector.y=-propagation_info.this_vector.y;
                        propagation_info.this_vector.z=-propagation_info.this_vector.z;

                      }



                      fibre_tract=new FibreTract();

                      tract_set_array[forward]->AddTract(fibre_tract);
                      fibre_tract->Release();

		      /*delete
                      voxel_point.x=(propagation_info.this_point.x-m_image_geom.x_start)/m_image_geom.x_step;
                      voxel_point.y=(propagation_info.this_point.y-m_image_geom.y_start)/m_image_geom.y_step;
                      voxel_point.z=(propagation_info.this_point.z-m_image_geom.z_start)/m_image_geom.z_step;
		      

                      world_point=VoxelToWorld(voxel_point.x,voxel_point.y,voxel_point.z,m_image_geom);
		      */
		      world_point=propagation_info.this_point;

                      //for FACT integration (currently FACT is the only int. type so always done):
                      propagation_info.this_plane=0;


                      tract_set_array[forward]->GetTract(0)->AddPoint(world_point);

                      tract_set_array[forward]->GetTract(0)->SetIndex(propagation_info.this_index,
                          tract_set_array[forward]->GetTract(0)->GetNumberOfPoints()-1);





                      if (flt_norm(propagation_info.this_vector)>FLT_EPSILON) {
                        Propagate(tract_set_array[forward],propagation_info);
                      }

                    } //forward/backward loop


                    //cout << "forward tracts: " << tract_set_array[0]->GetNumberOfTracts() << " forward points: " << tract_set_array[0]->GetTract(0)->GetNumberOfPoints() << endl;
                    //cout << "backward tracts: " << tract_set_array[1]->GetNumberOfTracts() << " backward points: " << tract_set_array[1]->GetTract(0)->GetNumberOfPoints() << endl;

                    //PruneBranchedSets just figures out whether to include/exclude tract based on include/exclude ROIs.
                    for (r=start_roi;r<m_number_of_rois;r++) {
                      PruneBranchedSets(tract_set_array,m_roi_array[r],false);
                    }



                    for (r=0;r<m_number_of_exclusion_masks;r++) {

                      PruneBranchedSets(tract_set_array,m_exclusion_mask_array[r],true);
                    }


                    //assign the cumulative min OTTs to the *tracts* by looking at their OTTs:
                    if (m_cum_min_OTT_map_on) {
                      AssignCumulativeMinOTTs(tract_set_array, m_roi_array[m_reference_roi]);
                    }

                    //CopySetsToFinalSet(tract_set_array,&m_fibre_tract_set,m_prune);
                    CopySetsToFinalSet(tract_set_array,m_prune);




                    //DO: note that reinit is not nec. if it already happened because the tracts were rejected in PruneBranchedSets...
                    //cout << "Final set tracts " << m_fibre_tract_set->GetNumberOfTracts() << endl;
                    //cout << "\nfinal reinit:\n";
                    //cout << tract_set_array[0]->GetNumberOfTracts() << " tracts\n";
                    tract_set_array[0]->ReInitialize();

                    //cout << tract_set_array[1]->GetNumberOfTracts() << " tracts\n\n";
                    tract_set_array[1]->ReInitialize();


                  }//dir counter loop
                } //iterations loop

              }//seeds per voxel
            }//seeds per voxel


          }//seeds per voxel

          if (m_tracking_algorithm==MultiMax) {
            directions->Release();
          }


        }//voxel
      }//voxel
    }//voxel
  }//voxel

  std::cout<<std::endl;

  delete(tract_set_array[0]);

  delete(tract_set_array[1]);
}


/*
check:  also need to implement the m_OTT for sigma_theta.

also can streamline the code so not freeing the basis vectors at every return
*/

VECTOR3D FibreTracking::GetNextVectorForLabelTracking(int x, int y, int z, PropagationInfo propagation_info)
{

 
//previous: NOTE: to track through the fan, can choose either (1) directions within a fan defined by two vectors, or (2) propagate everywhere in the input fibre ODF.  Pending is a third approach that will replace approach (1) by a fan of ellipsoidal cross section, the parameters input (not sure exactly yet).  Peter Savadjiev will provide the input for this. It can be looked at as a regularized (2).  doing (2) is somewhat problematic, as the fODFs are rarely 0 and we don't have a sense of what the small nonzero values mean.  Just threshold at 50%.  Peter's threshold takes care of these.  The rest of the values are considered legit.  Note that the planar segment could be considered the skeleton/medial surface of the fan.  
//Now that if we are using texture flow based fitting to a helicoid, we will switch to (2).

//the current implementation is (2):

  //cout << "GetNextVectorForLabelTracking\n";

  VECTOR3D vector;
  Directions *directions;
  float dotprod;
  float max_dot=0;
  int i;


  float cone_theta0;
  float cone_theta2;
  float cone_r;
  float r,theta;
  float xx,yy;
  bool mergeonly=false;



  if (m_ProbLabels0==NULL) { //does this ever happen?  if labels, it now is always prob.
    directions=m_ODF_data->GetODFMaxima(x,y,z);
  }


  //basis orthogonal to the direction:

  //float *basis1=(float *) malloc(3*sizeof(float));

  //float *basis2=(float *) malloc(3*sizeof(float));

  float orig_vector[3];

  float len;



  VECTOR3D crossprod;

  VECTOR3D projection;

  float anglewidth0,anglewidth2,angle;

  VECTOR3D vector0;

  VECTOR3D vector2;

  VECTOR3D polarity;



  m_label=GetLabel(x, y, z, propagation_info);
  //cout << "label " << m_label << endl;
 
  if (m_label==none) //need to stop tracking.  means no data exists.
    {
    vector0.x=0;
    vector0.y=0;
    vector0.z=0;
    return vector0;
    

  }

  //assert m_label= single, cross, or fan
  
  return GetNextVectorFromMultiMaxCone(x, y, z,propagation_info);


}



VECTOR3D FibreTracking::GetNextVectorFromMultiMaxCone(int x, int y, int z, PropagationInfo propagation_info)
{

  //cout << "GetNextVectorFromMultiMaxCone\n";

  VECTOR3D vector1,vector2,vector3;
  Directions *directions;
  float dotprod;
  float max_dot=0;
  int i;
  float cone_theta;

  int closest_index;


  float random;
  float ODF_mag=0;



    if (m_ProbLabels0!=NULL) {

    //test with fullpdf
      if (0)
	//if (m_it<4)
      {
	
	directions=m_fullpdf->GetODFMaxima(x,y,z);
	
	if (directions->GetNumberOfDirections()>4 || directions->GetNumberOfDirections()<1)
	  {
	    vector2.x=0;
	    vector2.y=0;
	    vector2.z=0;
	    return vector2;
	  }

      }
    else //assign the directions depending on the label
      
      {
      
	

      switch (m_label) {
      case single_curve:
        directions=new Directions(1);
	vector1.x=m_ProbLabels0->GetValue(x,y,z,3);
        vector1.y=m_ProbLabels0->GetValue(x,y,z,4);
        vector1.z=m_ProbLabels0->GetValue(x,y,z,5);
        directions->SetDirection(0,vector1);
        break;
      case doublecross:
        directions=new Directions(2);
        vector1.x=m_ProbLabels1->GetValue(x,y,z,3);
        vector1.y=m_ProbLabels1->GetValue(x,y,z,4);
        vector1.z=m_ProbLabels1->GetValue(x,y,z,5);
        directions->SetDirection(0,vector1);
        vector1.x=m_ProbLabels2->GetValue(x,y,z,3);
        vector1.y=m_ProbLabels2->GetValue(x,y,z,4);
        vector1.z=m_ProbLabels2->GetValue(x,y,z,5);
        directions->SetDirection(1,vector1);
        break;

      case triplecross:
        directions=new Directions(3);
        vector1.x=m_ProbLabels3->GetValue(x,y,z,3);
        vector1.y=m_ProbLabels3->GetValue(x,y,z,4);
        vector1.z=m_ProbLabels3->GetValue(x,y,z,5);
        directions->SetDirection(0,vector1);
        vector1.x=m_ProbLabels4->GetValue(x,y,z,3);
        vector1.y=m_ProbLabels4->GetValue(x,y,z,4);
        vector1.z=m_ProbLabels4->GetValue(x,y,z,5);
        directions->SetDirection(1,vector1);
        vector1.x=m_ProbLabels5->GetValue(x,y,z,3);
        vector1.y=m_ProbLabels5->GetValue(x,y,z,4);
        vector1.z=m_ProbLabels5->GetValue(x,y,z,5);
        directions->SetDirection(2,vector1);
        break;
      case fan:
        directions=new Directions(1);
	vector1.x=m_ProbLabels0->GetValue(x,y,z,3);
        vector1.y=m_ProbLabels0->GetValue(x,y,z,4);
        vector1.z=m_ProbLabels0->GetValue(x,y,z,5);
        directions->SetDirection(0,vector1);
        break;

      }


      }
    }

    ////////////////////////////////////////////////////
    else { //not prob - just Maxima.  (with cone - not being updated...)


      directions=m_ODF_data->GetODFMaxima(x,y,z);
      //cout << directions->GetDirection(0).x << " " << directions->GetDirection(0).y << " " << directions->GetDirection(0).z << endl;
    }
    ////////////////////////////////////////////////////

    //find closest direction to the ingoing one:
    //DO: can modify to branch
    //if the number of directions is 1, we skip this:
    
    if (directions->GetNumberOfDirections()==1) 
      {

      closest_index=0;


      } 
    else if (propagation_info.start==true) //this means there is no last_vector to consult, so randomize
      { 
	//for iterations 0-3(4), we hardcode to choose different maxima
	if (m_it==0)
	  {
	    closest_index=0;
	  }
	else if (m_it==1 && directions->GetNumberOfDirections()>1)
	  {
	    closest_index=1;
	  }
	else if (m_it==2 && directions->GetNumberOfDirections()>2)
	  {
	    closest_index=2;
	  }
	else if (m_it==3 && directions->GetNumberOfDirections()>3)//we never have 4 for prob
	  {
	    closest_index=3;
	  }
	else //random start in one of the orientations
	  {
	    random=(float)(rand() / ((double)(RAND_MAX)+(double)(1.0))); 


	    closest_index=0;

	    for (i=0;i<directions->GetNumberOfDirections();i++) {
	      if ((i/directions->GetNumberOfDirections())<random && random<=((i+1)/directions->GetNumberOfDirections())) {
		closest_index=i;
	      }
	    
	    }
	  }
      }
    else { //find closest:
      //cout << "options for closest:\n";
      for (i=0;i<directions->GetNumberOfDirections();i++) {
        dotprod=float_abs(propagation_info.last_vector.x*directions->GetXComponent(i)+propagation_info.last_vector.y*directions->GetYComponent(i)+propagation_info.last_vector.z*directions->GetZComponent(i));
	//cout << directions->GetXComponent(i) << " " << directions->GetYComponent(i) << " " << directions->GetZComponent(i) << endl;

        if (dotprod>max_dot) {
          max_dot=dotprod;
          closest_index=i;

        }

      }



    }

    if (m_ProbLabels0==NULL) //this is just multimax case; shouldn't end up here anymore (?) 
      {
      vector2=directions->GetDirection(closest_index); 
      m_OTT=FLT_MAX;
      }
    
    //take a random vector from the pdf for the fibre direction:


    if (m_ProbLabels0!=NULL) {
  
      //cout << m_label << endl;
      float cutoff=FLT_EPSILON;//this may need to change (?)
      int rand_ODF_dir;
      int max_index;
      //currently chooses a direction at random and then checks whether nonzero prob and keeps choosing until it finds one.
      //DO? should be random but based on the ODF? choose higher directions more often?? could create a histogram with dirs resampled proportional to confidence in them.  I think better if just random
      //DO? extract the nonzero dirs and then choose randomly from them: saves time ...
    

      while (ODF_mag<cutoff)
	{
	 

	  rand_ODF_dir=(rand()%(m_ProbLabels0->GetNSize()-6))+6;
	  
	  //test:
	  if (0)
	    //if (m_it<4)
	    {
	      vector2=directions->GetDirection(closest_index);
	      ODF_mag=1;//test: need to get value of max of pdf...
	    }
	  else 
	    {
	  switch (m_label) {
	  case single_curve:
	    //cone_theta=m_ProbLabels->GetValue(x,y,z,1);//won't get used
	    m_scale=m_ProbLabels0->GetValue(x,y,z,1);

	    if (m_it<4)//we take the maximum the first three its
	      // if (0)
	      {
		
	

		vector2.x=m_ProbLabels0->GetValue(x,y,z,3);
		vector2.y=m_ProbLabels0->GetValue(x,y,z,4);
		vector2.z=m_ProbLabels0->GetValue(x,y,z,5);
		//ODF_mag=1;

		
		ODF_mag=m_ProbLabels0->GetMaximumValue(x,y,z,max_index,6,m_ProbLabels0->GetNSize());
		//vector2=m_ProbLabels0->GetDirections()->GetDirection(max_index);
	      }
	    else 
	      {
		ODF_mag=m_ProbLabels0->GetValue(x,y,z,rand_ODF_dir);
		vector2=m_ProbLabels0->GetDirections()->GetDirection(rand_ODF_dir);
	      }
	    break;
	  case doublecross:

	    switch (closest_index) {
	    case 0:
	      //cone_theta=m_ProbLabels->GetValue(x,y,z,6);

	      m_scale=m_ProbLabels1->GetValue(x,y,z,1);
	
	      if (m_it<4)
		//if (0)
		{
		  
		  
		  vector2.x=m_ProbLabels1->GetValue(x,y,z,3);
		  vector2.y=m_ProbLabels1->GetValue(x,y,z,4);
		  vector2.z=m_ProbLabels1->GetValue(x,y,z,5);
		  //ODF_mag=1;

		 
		  ODF_mag=m_ProbLabels1->GetMaximumValue(x,y,z,max_index,6,m_ProbLabels0->GetNSize());
		  //vector2=m_ProbLabels1->GetDirections()->GetDirection(max_index);
		}
	      else 
		{
		  ODF_mag=m_ProbLabels1->GetValue(x,y,z,rand_ODF_dir);
		  vector2=m_ProbLabels1->GetDirections()->GetDirection(rand_ODF_dir);
		}
	      break;
	    case 1:
	      //cone_theta=m_ProbLabels->GetValue(x,y,z,11);
	      m_scale=m_ProbLabels2->GetValue(x,y,z,1);
	
	      if (m_it<4)
		//if (0)
		{
		  
		  
		  vector2.x=m_ProbLabels2->GetValue(x,y,z,3);
		  vector2.y=m_ProbLabels2->GetValue(x,y,z,4);
		  vector2.z=m_ProbLabels2->GetValue(x,y,z,5);
		  //ODF_mag=1;

		  
		  ODF_mag=m_ProbLabels2->GetMaximumValue(x,y,z,max_index,6,m_ProbLabels0->GetNSize());
		  //vector2=m_ProbLabels2->GetDirections()->GetDirection(max_index);
		  
		}
	      else 
		{
		  ODF_mag=m_ProbLabels2->GetValue(x,y,z,rand_ODF_dir);
		  vector2=m_ProbLabels2->GetDirections()->GetDirection(rand_ODF_dir);
		}
	      
	      break;
	    }

	    break;

	  case triplecross:


	    switch (closest_index) {
	    case 0:
	      //cone_theta=m_ProbLabels->GetValue(x,y,z,16);
	      m_scale=m_ProbLabels3->GetValue(x,y,z,1);
	   
	      if (m_it<4)
		//if (0)
		{
		  
		  
		  vector2.x=m_ProbLabels3->GetValue(x,y,z,3);
		  vector2.y=m_ProbLabels3->GetValue(x,y,z,4);
		  vector2.z=m_ProbLabels3->GetValue(x,y,z,5);
		  //ODF_mag=1;

		 
		  ODF_mag=m_ProbLabels3->GetMaximumValue(x,y,z,max_index,6,m_ProbLabels0->GetNSize());
		  //vector2=m_ProbLabels3->GetDirections()->GetDirection(max_index);
		}
	      else 
		{
		  ODF_mag=m_ProbLabels3->GetValue(x,y,z,rand_ODF_dir);
		  vector2=m_ProbLabels3->GetDirections()->GetDirection(rand_ODF_dir);
		}
	      
		//if (directions->GetNumberOfDirections...only if track w/ full pdf: means can't find closest so well
	      break;
	    case 1:
	      //cone_theta=m_ProbLabels->GetValue(x,y,z,21);
	      m_scale=m_ProbLabels4->GetValue(x,y,z,1);
	  
	      if (m_it<4)
		//if (0)
		{
		  
		  
		  vector2.x=m_ProbLabels4->GetValue(x,y,z,3);
		  vector2.y=m_ProbLabels4->GetValue(x,y,z,4);
		  vector2.z=m_ProbLabels4->GetValue(x,y,z,5);
		  //ODF_mag=1;

		 
		  ODF_mag=m_ProbLabels4->GetMaximumValue(x,y,z,max_index,6,m_ProbLabels0->GetNSize());
		  //vector2=m_ProbLabels4->GetDirections()->GetDirection(max_index);
		}
	      else 
		{
		  ODF_mag=m_ProbLabels4->GetValue(x,y,z,rand_ODF_dir);
		  vector2=m_ProbLabels4->GetDirections()->GetDirection(rand_ODF_dir);
		}

	      break;
	    case 2:
	      //cone_theta=m_ProbLabels->GetValue(x,y,z,26);
	      m_scale=m_ProbLabels5->GetValue(x,y,z,1);
	  
	      if (m_it<4)
		//if (0)
		{
		  
		  
		  vector2.x=m_ProbLabels5->GetValue(x,y,z,3);
		  vector2.y=m_ProbLabels5->GetValue(x,y,z,4);
		  vector2.z=m_ProbLabels5->GetValue(x,y,z,5);
		  //ODF_mag=1;

		 
		  ODF_mag=m_ProbLabels5->GetMaximumValue(x,y,z,max_index,6,m_ProbLabels0->GetNSize());
		  //vector2=m_ProbLabels5->GetDirections()->GetDirection(max_index);
		}
	      else 
		{
		  ODF_mag=m_ProbLabels5->GetValue(x,y,z,rand_ODF_dir);
		  vector2=m_ProbLabels5->GetDirections()->GetDirection(rand_ODF_dir);
		}

	      break;
	    }
	    
	    break;
	  case fan:
	     
	    
	    //todo: scaling needs to come from occurrence rate of fanning; see curve inference output
	    m_scale=1.0;
	
	      if (m_it<4 && m_ProbLabels1->GetValue(x,y,z,3)+m_ProbLabels1->GetValue(x,y,z,4)+m_ProbLabels1->GetValue(x,y,z,5)>0)
		//if (0)
		{
		  
		  
		  vector2.x=m_ProbLabels1->GetValue(x,y,z,3);
		  vector2.y=m_ProbLabels1->GetValue(x,y,z,4);
		  vector2.z=m_ProbLabels1->GetValue(x,y,z,5);
		  //ODF_mag=1;

		 
		  ODF_mag=m_ProbLabels1->GetMaximumValue(x,y,z,max_index,6,m_ProbLabels0->GetNSize());
		  //vector2=m_ProbLabels1->GetDirections()->GetDirection(max_index);
		}
	      else 
		{
		  //cout << "m_scale " << m_scale << endl;
		  ODF_mag=m_fullpdf->GetValue(x,y,z,rand_ODF_dir-6);
		  vector2=m_fullpdf->GetDirections()->GetDirection(rand_ODF_dir);
		}
	    break;
	  }
	  
	    }
	}
      if (m_it>=4)
	    {
	      //cout << "random ODF dir " << rand_ODF_dir << " ODF mag " << ODF_mag << endl;
	    }
      //cout << ODF_mag << endl;
      //test:
      //cout << vector2.x << " " << vector2.y << " " << vector2.z << endl;
      

    } else {
      m_scale=1.0;

      if (m_sigma_theta_const>FLT_EPSILON) {
	cone_theta=m_sigma_theta_const;
      } else {
	cone_theta=GetSigmaTheta(x,y,z).sigma_theta[closest_index];
      }
    }


    
  



  


  if (m_ProbLabels0!=NULL) 
    {
      //note old data wasn't normalized to unit volume; need anther scaling factor if so
      m_OTT=m_scale*ODF_mag;
      //cout << "m_OTT " << m_OTT << endl;
      //directions->Release(); 
    }

  directions->Release(); 
  

  vector3=Jitter(vector2); //random within cone around this sampling dir
  return vector3;






}

//reworked for the 6-file histo case:
void FibreTracking::SetProbLabels(ImageData *prob_labels0, ImageData *prob_labels1, ImageData *prob_labels2, ImageData *prob_labels3, ImageData *prob_labels4, ImageData *prob_labels5)
{

  //cout << "SetProbLabels\n";
  //get sampling from first minc file:
  //assume all the same: up to user


  m_image_geom.x_size=prob_labels0->GetXSize();
  m_image_geom.y_size=prob_labels0->GetYSize();
  m_image_geom.z_size=prob_labels0->GetZSize();
  m_image_geom.x_start=prob_labels0->GetXStart();
  m_image_geom.y_start=prob_labels0->GetYStart();
  m_image_geom.z_start=prob_labels0->GetZStart();
  m_image_geom.x_step=prob_labels0->GetXStep();
  m_image_geom.y_step=prob_labels0->GetYStep();
  m_image_geom.z_step=prob_labels0->GetZStep();
  m_image_geom.x_direction_cosine=prob_labels0->GetXDirectionCosine();
  m_image_geom.y_direction_cosine=prob_labels0->GetYDirectionCosine();
  m_image_geom.z_direction_cosine=prob_labels0->GetZDirectionCosine();


  m_ProbLabels0=prob_labels0;
  m_ProbLabels1=prob_labels1;
  m_ProbLabels2=prob_labels2;
  m_ProbLabels3=prob_labels3;
  m_ProbLabels4=prob_labels4;
  m_ProbLabels5=prob_labels5;

  m_cone=true;


}

VECTOR3D FibreTracking::GetNextVectorFromODF(int x, int y, int z)
{

  //need to look at this piece of code: does it check whether the ODF value is zero in the direction selected, or just propagate anywhere?  what all goes on in this function?  I haven't been testing the -ODF option.

  VECTOR3D vector;
  int r;
  int n;
  int i,j;
  int samples;
  float xx,yy,zz;
  double rd;
  float *basis1=(float *) malloc(3*sizeof(float));
  float *basis2=(float *) malloc(3*sizeof(float));
  float orig_vector[3];
  float radius;

  //DO: change this so guaranteed to work with up to 500 samples, if not possible memorywise, will have to recalculate each time reach voxel
  float sample_factor=100/sqrt(m_ODF_data->GetNSize());




  if (m_fibre_tract_set->GetBinaryConnectivityMap()->GetValue(x,y,z)<FLT_EPSILON) {
    m_ODF_sampler[z*m_image_geom.x_size*m_image_geom.y_size+y*m_image_geom.x_size+x]=new Directions(0);
    //calculate the sampled directions:

    //the circle to sample on the plane:
    //radius=sqrt(M_PI/m_ODF_data->GetNSize());
    //check: (this is approx of course)
    radius=sqrt(2.0/m_ODF_data->GetNSize());

    for (n=0;n<m_ODF_data->GetNSize();n++) {
      //determine the number of samples to take in this region.  BIG assumption: that can approximate the region by a cartesian plane.
      //DO: play with factor
      //can make differences btwn high and low diffusion greater by making nonlinear or subtracting min nonzero or increasing factor or using Stretch or something..

      //samples=(int) (m_ODF_data->GetValue(x,y,z,n)*sample_factor);

      //best so far:
      //DO: need to create a LUT that is reasonable for the ODF in question
      //samples=(int) (pow(m_ODF_data->GetValue(x,y,z,n),4.0)*sample_factor);

      //test: linear mapping from ODF range (lowest nonzero to max) to the range 1-sample_factor

      samples=(int)((m_ODF_data->GetValue(x,y,z,n)-m_ODF_data->GetMinimumNonZeroValue(x,y,z))*(sample_factor-1)/(m_ODF_data->GetMaximumValue(x,y,z)-m_ODF_data->GetMinimumNonZeroValue(x,y,z))+1);



      //debug:

      if (m_ODF_data->GetValue(x,y,z,n)>FLT_EPSILON && samples < FLT_EPSILON) {

        //cout << "insufficient sampling: " << samples << endl;
        //samples=10;

      }



      //samples=10;

      //end debug

      //define the plane:
      orig_vector[0]=m_ODF_data->GetDirections()->GetXComponent(n);

      orig_vector[1]=m_ODF_data->GetDirections()->GetYComponent(n);

      orig_vector[2]=m_ODF_data->GetDirections()->GetZComponent(n);

      MakeTwoOrthogonalVectors(orig_vector,basis1,basis2);

      for (i=0;i<samples;i++) {
        for (j=0;j<samples;j++) {

          //define the point on plane
          xx=orig_vector[0]+(radius*(-1.0+1.0/(float)samples+i*2.0/samples))*basis1[0]+(radius*(-1.0+1.0/(float)samples+j*2.0/samples))*basis2[0];
          yy=orig_vector[1]+(radius*(-1.0+1.0/(float)samples+i*2.0/samples))*basis1[1]+(radius*(-1.0+1.0/(float)samples+j*2.0/samples))*basis2[1];
          zz=orig_vector[2]+(radius*(-1.0+1.0/(float)samples+i*2.0/samples))*basis1[2]+(radius*(-1.0+1.0/(float)samples+j*2.0/samples))*basis2[2];

          //if euclidian distance between the point and the vector lies within the circle:

          if (sqrt(pow(orig_vector[0]-xx,2.0)+pow(orig_vector[1]-yy,2.0)+pow(orig_vector[2]-zz,2.0))<radius) {

            m_ODF_sampler[z*m_image_geom.x_size*m_image_geom.y_size+y*m_image_geom.x_size+x]->IncreaseNumberOfDirections(1);
            m_ODF_sampler[z*m_image_geom.x_size*m_image_geom.y_size+y*m_image_geom.x_size+x]->SetXComponent(m_ODF_sampler[z*m_image_geom.x_size*m_image_geom.y_size+y*m_image_geom.x_size+x]->GetNumberOfDirections()-1, xx);
            m_ODF_sampler[z*m_image_geom.x_size*m_image_geom.y_size+y*m_image_geom.x_size+x]->SetYComponent(m_ODF_sampler[z*m_image_geom.x_size*m_image_geom.y_size+y*m_image_geom.x_size+x]->GetNumberOfDirections()-1, yy);
            m_ODF_sampler[z*m_image_geom.x_size*m_image_geom.y_size+y*m_image_geom.x_size+x]->SetZComponent(m_ODF_sampler[z*m_image_geom.x_size*m_image_geom.y_size+y*m_image_geom.x_size+x]->GetNumberOfDirections()-1, zz);
            //cout << "xx: " << xx << " yy: " << yy << " zz: " << zz<< endl;


          }

        }


      }


    }


    m_ODF_sampler[z*m_image_geom.x_size*m_image_geom.y_size+y*m_image_geom.x_size+x]->Normalize();

    //plot the samples

    //DisplayWindow *display_window=new DisplayWindow();
    //display_window->AddObject(m_ODF_sampler[z*m_image_geom.x_size*m_image_geom.y_size+y*m_image_geom.x_size+x]->GetActor());
    //display_window->Display();

  //and exit:

    //Error("just checking");


  }

  //generate a random number and draw the direction:
  rd=((double)rand() / ((double)(RAND_MAX)+(double)(1)));

  r=(int)(rd*(m_ODF_sampler[z*m_image_geom.x_size*m_image_geom.y_size+y*m_image_geom.x_size+x]->GetNumberOfDirections()-1));
  //test:

  //r=counter%(m_ODF_sampler[z*m_image_geom.x_size*m_image_geom.y_size+y*m_image_geom.x_size+x]->GetNumberOfDirections()-1);

  //counter+=1;



  //cout << r << endl;


  vector.x=m_ODF_sampler[z*m_image_geom.x_size*m_image_geom.y_size+y*m_image_geom.x_size+x]->GetXComponent(r);

  vector.y=m_ODF_sampler[z*m_image_geom.x_size*m_image_geom.y_size+y*m_image_geom.x_size+x]->GetYComponent(r);

  vector.z=m_ODF_sampler[z*m_image_geom.x_size*m_image_geom.y_size+y*m_image_geom.x_size+x]->GetZComponent(r);

  free(basis1);

  free(basis2);



  return vector;

}


VECTOR3D FibreTracking::GetNextVector(int x, int y, int z, bool swap, PropagationInfo propagation_info)
{
  VECTOR3D vectorv;

  switch (m_tracking_algorithm) {
  case PDD:
    vectorv=GetNextVectorFromPDDCone(x,y,z);
    break;

  case ODF_MC:
    vectorv= GetNextVectorFromODF(x,y,z);
    break;

  case MultiMax:

    //vectorv=GetNextVectorFromMultiMaxCone(x,y,z,propagation_info);
    vectorv=GetNextVectorFromMultiMax(x,y,z,propagation_info);//no more sigma_theta
    break;

  case Labels:
    vectorv=GetNextVectorForLabelTracking(x,y,z,propagation_info);
    //cout << vectorv.x << " " << vectorv.y << " " << vectorv.z << endl;
    break;
  }

  if (!propagation_info.start) {
    if (swap && (propagation_info.last_vector.x*vectorv.x+propagation_info.last_vector.y*vectorv.y+propagation_info.last_vector.z*vectorv.z)<0) {
      //cout << "swap dir: last: " << propagation_info.last_vector.x << " " << propagation_info.last_vector.y << " " << propagation_info.last_vector.z << "this: " << vectorv.x << " " << vectorv.y << " " << vectorv.z << endl;
      vectorv.x=-vectorv.x;
      vectorv.y=-vectorv.y;
      vectorv.z=-vectorv.z;
    }
  }

  //cout << "vector " << vectorv.x << " " << vectorv.y << " " << vectorv.z << " " << x << " " << y << " " << z << endl;

  return vectorv;

}

//this is only for Prob Labels so far:
VECTOR3D FibreTracking::Jitter(VECTOR3D v)
{
  float cone_theta=sqrt(2*M_PI/(m_ProbLabels0->GetNSize()-6));
  VECTOR3D vector2;

    //now generate a random vector within the cone:



    float cone_r=sin(cone_theta);
    //basis orthogonal to the direction:

    float *basis1=(float *) malloc(3*sizeof(float));
    float *basis2=(float *) malloc(3*sizeof(float));
    float orig_vector[3];


    orig_vector[0]=v.x;
    orig_vector[1]=v.y;
    orig_vector[2]=v.z;

    MakeTwoOrthogonalVectors(orig_vector,basis1,basis2);

    float r,theta;
    float xx,yy;

    r=(float)((double)rand() / ((double)(RAND_MAX)+(double)(1.0)))*cone_r;
    theta=(float)((double)rand() / ((double)(RAND_MAX)+(double)(1.0)))*2*M_PI;

    xx=r*cos(theta);
    yy=r*sin(theta);
    vector2.x=v.x+(xx*basis1[0]+yy*basis2[0]);
    vector2.y=v.y+(xx*basis1[1]+yy*basis2[1]);
    vector2.z=v.z+(xx*basis1[2]+yy*basis2[2]);
    float len;
    len=sqrt(vector2.x*(vector2.x)+vector2.y*(vector2.y)+vector2.z*(vector2.z));
    vector2.x=vector2.x/len;
    vector2.y=vector2.y/len;
    vector2.z=vector2.z/len;

    free(basis1);
    free(basis2);

   
    return vector2;

}

VECTOR3D FibreTracking::GetNextVectorFromPDDCone(int x, int y, int z)
{
  VECTOR3D vector,vector2;
  //SIGMA_THETA cone_theta;
  float cone_theta;

  if (m_sigma_theta_const>FLT_EPSILON) {
    cone_theta=m_sigma_theta_const;
  } else {
    cone_theta=GetSigmaTheta(x,y,z).sigma_theta[0];
  }

  if (m_it==0 || cone_theta<FLT_EPSILON) {
    vector.x=m_processed_DTI_data->Gete1x()->GetValue(x,y,z);
    vector.y=m_processed_DTI_data->Gete1y()->GetValue(x,y,z);
    vector.z=m_processed_DTI_data->Gete1z()->GetValue(x,y,z);

    if (cone_theta>FLT_EPSILON) {
      m_OTT=normal(0.0,cone_theta,M_PI/2);
      //m_OTT=1.0-0.0/(M_PI/2);
    }

    //no code for sigma_theta=0 here: don't expect it.

    return vector;
  }

  else {
    vector.x=m_processed_DTI_data->Gete1x()->GetValue(x,y,z);
    vector.y=m_processed_DTI_data->Gete1y()->GetValue(x,y,z);
    vector.z=m_processed_DTI_data->Gete1z()->GetValue(x,y,z);

    //now generate a random vector within the cone:



    float cone_r=sin(cone_theta);
    //basis orthogonal to the direction:

    float *basis1=(float *) malloc(3*sizeof(float));
    float *basis2=(float *) malloc(3*sizeof(float));
    float orig_vector[3];


    orig_vector[0]=vector.x;
    orig_vector[1]=vector.y;
    orig_vector[2]=vector.z;

    MakeTwoOrthogonalVectors(orig_vector,basis1,basis2);

    float r,theta;
    float xx,yy;

    r=(float)((double)rand() / ((double)(RAND_MAX)+(double)(1.0)))*cone_r;
    theta=(float)((double)rand() / ((double)(RAND_MAX)+(double)(1.0)))*2*M_PI;

    xx=r*cos(theta);
    yy=r*sin(theta);
    vector2.x=vector.x+(xx*basis1[0]+yy*basis2[0]);
    vector2.y=vector.y+(xx*basis1[1]+yy*basis2[1]);
    vector2.z=vector.z+(xx*basis1[2]+yy*basis2[2]);
    float len;
    len=sqrt(vector2.x*(vector2.x)+vector2.y*(vector2.y)+vector2.z*(vector2.z));
    vector2.x=vector2.x/len;
    vector2.y=vector2.y/len;
    vector2.z=vector2.z/len;

    free(basis1);
    free(basis2);

    m_OTT=normal(float_abs(acos(dot(vector,vector2))),cone_theta,M_PI/2);
    //m_OTT=1.0-float_abs(acos(dot(vector,vector2)))/(M_PI/2);

    /*
    if (m_OTT>0)
    {
    cout << "ott " << m_OTT << endl;
    cout << "vector " << vector.x << vector.y << vector.z << endl;
    cout << "vector2 " << vector2.x << vector2.y << vector2.z << endl;
    }
    */

    return vector2;

  }

}

VECTOR3D FibreTracking::GetNextVectorFromMultiMax(int x, int y, int z, PropagationInfo propagation_info)
{
  VECTOR3D vector;
  Directions *directions;
  float dot;
  float max_dot=0;
  int i;

  directions=m_ODF_data->GetODFMaxima(x,y,z);


  //note that this should never get called if on first step, and hence last vector should always be set.

  

  //find closest direction to the ingoing one:

  for (i=0;i<directions->GetNumberOfDirections();i++) {
    dot=float_abs(propagation_info.last_vector.x*directions->GetXComponent(i)+propagation_info.last_vector.y*directions->GetYComponent(i)+propagation_info.last_vector.z*directions->GetZComponent(i));

    if (dot>max_dot) {
      max_dot=dot;

      vector.x=directions->GetXComponent(i);
      vector.y=directions->GetYComponent(i);
      vector.z=directions->GetZComponent(i);
    }
  }


  //cout << "vector from multi max: " <<  vector.x << " " << vector.y << " " << vector.z << endl;
  //directions->Release();



  return vector;


}



void FibreTracking::SetPruneOn()
{
  //JC changed this so it gets set regardless of whether save tracts is on at this instant.

  m_prune=true;

}

void FibreTracking::SetMCIterations(int iterations)
{
  m_MC_iterations=iterations;

}


void FibreTracking::SetLabelsOn(char *labels_filename)
{
  m_labels_on=true;
  m_labels_exist=true;
  m_labels=new ImageData(labels_filename);

}

ImageData *FibreTracking::GetBinaryTrackingMask()
{
  return m_mask;
}



Label FibreTracking::GetLabel(int x, int y, int z, PropagationInfo propagation_info)
{
  
  //cout << "GetLabel\n";
  

  float confidence=0.0;
  Label label;
  float random;
  Label thisoldlabel,thisnewlabel;

  if (!m_labels_on && m_ProbLabels0==NULL) {
    Error("no labels");
  } else if (m_ProbLabels0!=NULL) {

    //a check for wrong voxels -> stop:
    if ((m_ProbLabels0->GetValue(x,y,z,0)<FLT_EPSILON || m_ProbLabels0->GetValue(x,y,z,0)>8) && (m_ProbLabels1->GetValue(x,y,z,0)<FLT_EPSILON || m_ProbLabels1->GetValue(x,y,z,0)>8) && (m_ProbLabels3->GetValue(x,y,z,0)<FLT_EPSILON || m_ProbLabels3->GetValue(x,y,z,0)>8))
      {
	//cout << x << " " << y << " " << z << endl;
	return none;
      
      }


    //will return one of single_curve, fan,doublecross,triplecross: one that has confidence >0


    //on first iteration, return most likely geom:
    
    if (m_it<3)
      {
	
	 if (m_ProbLabels0->GetValue(x,y,z,0)>=m_ProbLabels1->GetValue(x,y,z,0) && m_ProbLabels0->GetValue(x,y,z,0)>=m_ProbLabels3->GetValue(x,y,z,0))
	  {
	    label=single_curve;
	    //cout << "first 3 its label " << label << " " << m_ProbLabels0->GetValue(x,y,z,0) << " " << m_ProbLabels1->GetValue(x,y,z,0) << " " << m_ProbLabels3->GetValue(x,y,z,0) << " " << x << " " << y << " " << z << endl;
	  }
	else if (m_ProbLabels1->GetValue(x,y,z,0)>m_ProbLabels0->GetValue(x,y,z,0) && m_ProbLabels1->GetValue(x,y,z,0)>=m_ProbLabels3->GetValue(x,y,z,0))
	  {
	    label=doublecross;
	    //cout << "first 3 its label " << label << endl;
	  }
	else if (m_ProbLabels3->GetValue(x,y,z,0)>m_ProbLabels0->GetValue(x,y,z,0) && m_ProbLabels3->GetValue(x,y,z,0)>m_ProbLabels1->GetValue(x,y,z,0))
	  {
	    label=triplecross;
	    //cout << "first 3 its label " << label << endl;
	  }

      }
  
    else //not first three its
      {

	label=none;
    

      while (confidence<FLT_EPSILON){
        random=(float)((double)rand() / ((double)(RAND_MAX)+(double)(1.0)));

        if (random<0.33) {
         confidence=m_ProbLabels0->GetValue(x,y,z,0);

        if (confidence>FLT_EPSILON) {
	  label=single_curve;
        }
	}

	else if (random<0.66) {
        
        confidence=float_min(m_ProbLabels1->GetValue(x,y,z,0),m_ProbLabels2->GetValue(x,y,z,0));
	
        if (confidence>FLT_EPSILON) {

	  if (confidence<0.2)//if low confidence and no uncert was detected:
	    {
	      if (m_ProbLabels1->GetValue(x,y,z,2)<FLT_EPSILON || m_ProbLabels2->GetValue(x,y,z,2)<FLT_EPSILON)
		{
		  confidence=0;
		  //cout << "doublecross sigma_theta==0; confidence " << m_ProbLabels->GetValue(x,y,z,5) << " " << m_ProbLabels->GetValue(x,y,z,10) << endl;

		}
	      }
	
	  label=doublecross;
	}
	}
	else
	  {
	    //DO: retain the ones that have nonzero sigma_theta: relabel: requires a lot of checking elsewhere...
	    confidence=float_min3(m_ProbLabels3->GetValue(x,y,z,0),m_ProbLabels4->GetValue(x,y,z,0),m_ProbLabels5->GetValue(x,y,z,0));
	    if (confidence>FLT_EPSILON)
	      {
		if (confidence<0.2)//if low confidence and sigma_theta==0
		  {
			if ((m_ProbLabels3->GetValue(x,y,z,2)<FLT_EPSILON || m_ProbLabels4->GetValue(x,y,z,2)<FLT_EPSILON) || m_ProbLabels5->GetValue(x,y,z,2)<FLT_EPSILON)
			  {
			    confidence=0;
			    //cout << "triplecross sigma_theta==0; confidence " << m_ProbLabels->GetValue(x,y,z,15) << " " << m_ProbLabels->GetValue(x,y,z,20) << " " << m_ProbLabels->GetValue(x,y,z,25) << endl;

			  }
			}
		    label=triplecross;
		  }
	      }

      }
      }


    //if we have a nonzero polarity vector and dot product with polarity positive, we fan
    //todo: get fanning prob from curve inference output and use it to determine how often to fan; orig output had fan occurrence=100%!
    //however, we need to not fan sometimes, because the fan probability it wrong and I'm not using it, so fan 50% of the time:

    if (m_fanning)
      {
	
	if (propagation_info.this_vector.x*m_polarity_vector->GetValue(x,y,z,0)+propagation_info.this_vector.y*m_polarity_vector->GetValue(x,y,z,1)+propagation_info.this_vector.z*m_polarity_vector->GetValue(x,y,z,2)>0 && m_it>m_MC_iterations/2)
	  {
	    label=fan;
	  }
      }


    thisnewlabel=label;
    //cout << "new label " << thisnewlabel << endl;
    return thisnewlabel;
      
      
  }
  else if (0) 
    {
      // always need new label?
      //else if (m_labels_exist) {
      //cout << m_labels->GetValue(x,y,z) << endl;
      thisoldlabel=(Label)(int) m_labels->GetValue(x,y,z); //two casts because it started off as a float
      //cout << "old label " << thisoldlabel << endl;
      return thisoldlabel;
    }

  //else //branch everywhere: we don't do this anymore
  //{
  //  return 4;
  //}
}




//void FibreTracking::CopySetsToFinalSet(FibreTractSet *tract_set_array[2], FibreTractSet *(*final_set), bool prune_same_indices)
void FibreTracking::CopySetsToFinalSet(FibreTractSet *tract_set_array[2], bool prune_same_indices)
{
  FibreTractSet *fibre_tract_set_tmp;
  int forward;
  int t;
  FibreTract *fibre_tract;
  int i;
  INDEX3D index;


  for (forward=0;forward<2;forward++) {
    for (t=0;t<tract_set_array[forward]->GetNumberOfTracts();t++) {
      fibre_tract=tract_set_array[forward]->GetTract(t);


      if (fibre_tract->GetNumberOfPoints()>1) {
        //if ((*final_set)->SaveTractsOn())
        if (m_fibre_tract_set->SaveTractsOn()) {
          m_fibre_tract_set->AddTract(fibre_tract);
          //(*final_set)->AddTract(fibre_tract);
        }

        //assign scalar maps here.  The were never set for the individual forward/backward sets (pruning needed them together).

        for (i=0;i<fibre_tract->GetNumberOfPoints();i++) {
          index=fibre_tract->GetIndex(i);

          if (m_binary_connectivity_map_on) {
            //(*final_set)->GetBinaryConnectivityMap()->SetValue(index.x,index.y,index.z,1.0);
            m_fibre_tract_set->GetBinaryConnectivityMap()->SetValue(index.x,index.y,index.z,1.0);

          }

          /*
          if (m_min_DOTT_map_on)
            {
             (*final_set)->GetMinDOTTMap()->SetValue(index.x,index.y,index.z,float_max(tract_set_array[forward]->GetMinDOTTMap()->GetValue(index.x,index.y,index.z),(*final_set)->GetMinDOTTMap()->GetValue(index.x,index.y,index.z)));
            }
          */
          if (m_cum_min_OTT_map_on) {
            //we could look at the scalar map IF after the pruning because the pruning (same indices) didn't guarantee the most likely tract...
            //(*final_set)->GetCumMinOTTMap()->SetValue(index.x,index.y,index.z,float_max(tract_set_array[forward]->GetCumMinOTTMap()->GetValue(index.x,index.y,index.z),(*final_set)->GetCumMinOTTMap()->GetValue(index.x,index.y,index.z)));
            //this should already have been doen on the fly if not brute force:
            if (m_brute_force || m_reference_roi>0) {
              m_fibre_tract_set->GetCumMinOTTMap()->SetValue(index.x,index.y,index.z,float_max(tract_set_array[forward]->GetCumMinOTTMap()->GetValue(index.x,index.y,index.z),m_fibre_tract_set->GetCumMinOTTMap()->GetValue(index.x,index.y,index.z)));
            }

            //(*final_set)->GetCumMinOTTMap()->SetValue(index.x,index.y,index.z,float_max(fibre_tract->GetPointInfo(i).cumulative_min_OTT,(*final_set)->GetCumMinOTTMap()->GetValue(index.x,index.y,index.z)));
          }

        }

        //fibre_tract->Release(); //not nec because we didn't retain it?


        //we prune AFTER set scalar indices, because would miss stuff if not
        //note: it would be safe (and faster) to do the binary map after pruning, provided we knew whether it got kept or not...PruneSameIndices return boolean - not implemented.
        //if (prune_same_indices && (*final_set)->SaveTractsOn())
        if (prune_same_indices && m_fibre_tract_set->SaveTractsOn()) {
          fibre_tract_set_tmp=m_fibre_tract_set->PruneSameIndices(true);
          //fibre_tract_set_tmp=(*final_set)->PruneSameIndices(true);
          delete(m_fibre_tract_set);
          m_fibre_tract_set=fibre_tract_set_tmp;
        }
      }


      else {

        //fibre_tract->Release(); //not nec because we didn't retain it
      }

    }


  }


}



void FibreTracking::PruneBranchedSets(FibreTractSet *tract_set_array[2], ImageData *roi, bool exclude)
{


  //this assumes the roi in question is only on one side or the other of the seed voxel.  In most cases this is true, but a funny-shaped (concave) ROI would cause problems.

  bool roi_here[2];
  int forward;
  int i,k,t;
  FibreTractSet *fibre_tract_set_tmp;

  if (!exclude) {
    for (forward=0;forward<2;forward++) {
      //go through set and see whether it goes through roi:
      roi_here[forward]=false;

      for (i=0;i<tract_set_array[forward]->GetNumberOfTracts() && !roi_here[forward];i++) {
        for (k=0;k<tract_set_array[forward]->GetTract(i)->GetNumberOfPoints() && !roi_here[forward];k++) {

          if (roi->GetValue(tract_set_array[forward]->GetTract(i)->GetIndex(k).x,tract_set_array[forward]->GetTract(i)->GetIndex(k).y,tract_set_array[forward]->GetTract(i)->GetIndex(k).z)>FLT_EPSILON) {
            roi_here[forward]=true;
          }
        }
      }
    }

    if (roi_here[0] && !roi_here[1]) {
      fibre_tract_set_tmp=tract_set_array[0]->Prune(roi,exclude,false);
      delete(tract_set_array[0]);
      tract_set_array[0]=fibre_tract_set_tmp;
    } else if (roi_here[1] && !roi_here[0]) {
      fibre_tract_set_tmp=tract_set_array[1]->Prune(roi,exclude,false);
      delete(tract_set_array[1]);
      tract_set_array[1]=fibre_tract_set_tmp;
    } else if (!roi_here[0] && !roi_here[1]) {
      //test:
      //cout << "reinit forward tract\n";


      tract_set_array[0]->ReInitialize();

      //cout << "reinit backward tract; tracts " << tract_set_array[1]->GetNumberOfTracts() << endl;





      tract_set_array[1]->ReInitialize();



      //delete(tract_set_array[0]);
      //delete(tract_set_array[1]);
      //tract_set_array[0]=new FibreTractSet();
      //tract_set_array[1]=new FibreTractSet();
    }

  }

  else if (exclude) {

    for (forward=0;forward<2;forward++) {
      //go through set and see whether it goes through roi:
      roi_here[forward]=false;

      for (i=0;i<tract_set_array[forward]->GetNumberOfTracts();i++) {
        for (k=0;k<tract_set_array[forward]->GetTract(i)->GetNumberOfPoints() && !roi_here[forward];k++) {

          if (roi->GetValue(tract_set_array[forward]->GetTract(i)->GetIndex(k).x,tract_set_array[forward]->GetTract(i)->GetIndex(k).y,tract_set_array[forward]->GetTract(i)->GetIndex(k).z)>FLT_EPSILON) {
            roi_here[forward]=true;
          }
        }
      }
    }

    if (roi_here[0] && !roi_here[1]) {
      fibre_tract_set_tmp=tract_set_array[0]->Prune(roi,exclude,false);
      delete(tract_set_array[0]);
      tract_set_array[0]=fibre_tract_set_tmp;
      //test:
      tract_set_array[1]->ReInitialize();
      //delete(tract_set_array[1]);
      //tract_set_array[1]=new FibreTractSet();

    } else if (roi_here[1] && !roi_here[0]) {
      fibre_tract_set_tmp=tract_set_array[1]->Prune(roi,exclude,false);
      delete(tract_set_array[1]);
      tract_set_array[1]=fibre_tract_set_tmp;
      //test:
      tract_set_array[0]->ReInitialize();
      //delete(tract_set_array[0]);
      //tract_set_array[0]=new FibreTractSet();
    }

    //I guess it could go through on both sides:
    else if (roi_here[1] && roi_here[1]) {

      fibre_tract_set_tmp=tract_set_array[0]->Prune(roi,exclude,false);

      delete(tract_set_array[0]);
      tract_set_array[0]=fibre_tract_set_tmp;
      fibre_tract_set_tmp=tract_set_array[1]->Prune(roi,exclude,false);

      delete(tract_set_array[1]);
      tract_set_array[1]=fibre_tract_set_tmp;
      //tract_set_array[0]->ReInitialize();
      //tract_set_array[1]->ReInitialize();
    }

    else if (!roi_here[0] && !roi_here[1]) {
      //keep as is
    }




  }
}

void FibreTracking::AssignCumulativeMinOTTs(FibreTractSet *tract_set_array[2], ImageData *roi)
{
  //this is done here because must consider forward and backward tracts, which are in different sets.
  //these sets are potentially branching
  //this could be sped up: currently goes through all tracts in the set from start to finish, but a lot of that is repetitive because they are branched tracts.  This is only an issue if branching (for all Maxima), which is rarely used at this point
  //note that the roi is potentially only in one of the forward and backward paths in the case of brute_force or reference_roi>1(otherwise in both)

  bool blur_along_tract=true; 
  bool roi_found;
  bool roi_here;
  int forward;
  int i,k,t,f,l;
  float CMOTT;
  FibreTract *tract;

  roi_found=false;
  roi_here=false;

  int roi_index;
  float sigma;

  forward=0;
  CMOTT=FLT_MAX;
  bool roi_second_half=false;
  bool roi_found_this_tract;
  char write_file;


  //temporary non bf code:

  if (!m_brute_force && m_reference_roi==0) { //we know we start in the ROI:

    for (forward=0;forward<2;forward++) {

      for (i=0;i<tract_set_array[forward]->GetNumberOfTracts();i++) {
        roi_found=false;
        tract=tract_set_array[forward]->GetTract(i);
        //cout << "points " << tract->GetNumberOfPoints() << endl;

        for (k=0;k<tract->GetNumberOfPoints();k++) {
          if (roi->GetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z)>FLT_EPSILON) {
            tract->SetCumulativeMinOTT(-10.0,k);//for writing tract coords
            m_fibre_tract_set->GetCumMinOTTMap()->SetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z,-10.0);

            roi_found=true;
            CMOTT=FLT_MAX;
            /*
            if (m_it==40)
            {
            //cout << "cmott " << CMOTT << endl;
            }
            */
          } else {
            assert(roi_found);//true if not brute force

            //cout << "CMOTT before" << CMOTT << endl;

            //we are getting nan here sometimes ???? voxel lineup problem?? the check does nothing here: FLT_MAX will still get assigned.

	   if (blur_along_tract)
	     {
	       if (float_min(CMOTT,tract->GetBlurredOTT(k))<FLT_MAX)//detect nan? which shouldn't happen.
		 {
		   //cout << "CMOTT " << CMOTT << " blurred " << tract->GetBlurredOTT(k) << " not blurred " << tract->GetOTT(k) << endl;
	       CMOTT=float_min(CMOTT,tract->GetBlurredOTT(k));
		 }
	     }
	   else
	     {
	       if (float_min(CMOTT,tract->GetOTT(k))<FLT_MAX)//detect nan? which shouldn't happen.
		 {
	       CMOTT=float_min(CMOTT,tract->GetOTT(k));
		 }
	     }
	     
	    /* 
            else//DO: check why this happens
            {
            //cout << "cmott " << float_min(CMOTT,tract->GetOTT(k)) << endl;
            }
            */
            /*
            if (m_it==40)
            {
            //cout << "cmott " << CMOTT << endl;
            }
            */
            //cout << "CMOTT after" << CMOTT << endl;
            tract->SetCumulativeMinOTT(CMOTT,k);//for if write the tracts out...

            //cout << "map value before: " << m_fibre_tract_set->GetCumMinOTTMap()->GetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z) << endl;

            //what about just setting the master map?  we are guaranteed to go through the seed, so:
            m_fibre_tract_set->GetCumMinOTTMap()->SetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z,float_max(CMOTT,m_fibre_tract_set->GetCumMinOTTMap()->GetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z)));

            /*
            if (m_it==40)
            {
            //cout << "map cmott " << float_max(CMOTT,m_fibre_tract_set->GetCumMinOTTMap()->GetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z)) << " " << tract->GetIndex(k).x << " " << tract->GetIndex(k).y << " " << tract->GetIndex(k).z << endl;
            }
            */

            //cout << "map value after: " << m_fibre_tract_set->GetCumMinOTTMap()->GetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z) << endl;




          }
        }
      }

      /*
      cout << "write intermediate file? [y] [n]\n";
      cin >> write_file;
      if (write_file=='y')
        {
          m_fibre_tract_set->GetCumMinOTTMap()->WriteMincFile("test.mnc",true);
        }
      */

    }
    return;
  }

  //brute force and reference_roi>0:

  for (forward=0;forward<2;forward++) {

    //go through set and see whether it goes through roi.  if so, continue to end, assigning the cumulative min OTT, and back up from beginning of ROI, assigning as well (if greater than previous).
    //HERE go through and find out why we only get CMOTT from roi2 on: need to back up along its segment AND the other half.
    for (i=0;i<tract_set_array[forward]->GetNumberOfTracts();i++) {
      //if we've already found the roi in the other direction, set CMOTT to that of the first point on that tract (implies l>0)
      if (roi_found) {
        CMOTT=tract_set_array[forward]->GetCumMinOTTMap()->GetValue(tract->GetIndex(0).x,tract->GetIndex(0).y,tract->GetIndex(0).z);
      } else {
        CMOTT=FLT_MAX;
      }

      tract=tract_set_array[forward]->GetTract(i);

      roi_found_this_tract=false;
      roi_index=-1;

      for (k=0;k<tract->GetNumberOfPoints();k++) {
        if (roi->GetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z)>FLT_EPSILON) {

          if (roi_index<0) {
            roi_index=k;
            CMOTT=FLT_MAX;
          } else {
            roi_index=-1;//only back up once, at beginning
          }

          roi_found=true;

          roi_found_this_tract=true;

          if (forward==1) {
            roi_second_half=true;
          }

          //assigning -10 inside the ROI: then search for the maximum CMOTT in the volume afterward and put something larger but on the same scale in the ROI...algorithm/acq/sigma dependant...

          tract->SetCumulativeMinOTT(-10.0,k);

          tract_set_array[forward]->GetCumMinOTTMap()->SetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z,-10.0);

          //for the 1.0-approach:
          //tract->SetCumulativeMinOTT(1.0,k);
          //tract_set_array[forward]->GetCumMinOTTMap()->SetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z,1.0);
        }

        else {

          if (roi_found_this_tract) { //we have now passed out of the roi
  	  //cout << "CMOTT before" << CMOTT << endl;
	      if (float_min(CMOTT,tract->GetOTT(k))<FLT_MAX)//detect nan? which shouldn't happen.
		{
		  if (blur_along_tract)
		    {
		      CMOTT=float_min(CMOTT,tract->GetBlurredOTT(k));

		    }
		  else
		    {
		      CMOTT=float_min(CMOTT,tract->GetOTT(k));
		    }

		}

            //cout << "CMOTT after" << CMOTT << endl;
            tract->SetCumulativeMinOTT(CMOTT,k);

            //cout << "map value before: " << tract_set_array[forward]->GetCumMinOTTMap()->GetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z) << endl;
            tract_set_array[forward]->GetCumMinOTTMap()->SetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z,float_max(CMOTT,tract_set_array[forward]->GetCumMinOTTMap()->GetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z)));

            //cout << "map value after: " << tract_set_array[forward]->GetCumMinOTTMap()->GetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z) << endl;

          }
        }

      }

      //now go backward from roi if it was found on this branch: as back up, reassign all if greater
      if (roi_index>-1) {
        CMOTT=FLT_MAX;

        for (k=roi_index;k>=0;k--) {
          if (float_min(CMOTT,tract->GetOTT(k))<FLT_MAX) { //detect nan
 		  if (float_min(CMOTT,tract->GetOTT(k))>tract->GetCumulativeMinOTT(k)){
		      if (blur_along_tract){
			  CMOTT=float_min(CMOTT,tract->GetBlurredOTT(k));
			}
		      else{
			  CMOTT=float_min(CMOTT,tract->GetOTT(k));
			}

		    }
		  else//we stop going backward
		    {
		      k=-1;
		    }
	    }

          if (k<roi_index && k>-1) { //if hit roi again, stop
            if (roi->GetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z)>FLT_EPSILON)

            {
              k=-1;
            }
          }

          if (k>=0) {
            tract->SetCumulativeMinOTT(CMOTT,k);
            tract_set_array[forward]->GetCumMinOTTMap()->SetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z,float_max(CMOTT,tract_set_array[forward]->GetCumMinOTTMap()->GetValue(tract->GetIndex(k).x,tract->GetIndex(k).y,tract->GetIndex(k).z)));
          }
        }
      }


    }//each branch (starting at seed)


    if (forward==1 && roi_second_half) {
      forward=-1;//will become 0 again for third round
    }

    if (forward==0 && roi_second_half) {
      forward=2;//force stop.
    }
  }//forward and backward

  //so if roi never found, nothing happens

}






void FibreTracking::SetSigmaTheta()
{
  //this assumes already have m_ODF_data set.  sigma_theta_const will be based on AR of ODF.
  m_sigma_theta_const=sqrt(2*M_PI/m_ODF_data->GetNSize())/2;
  //now set iterations if they aren't set:


  m_binary_connectivity_map_on=true;

  if (m_MC_iterations==1) {
    m_MC_iterations=10000;//user needs to take care to change this afterward if less/more desired
  }

  m_cone=true;

}

void FibreTracking::SetSigmaTheta(float sigma_theta_const)
{
  //this assumes already have m_ODF_data set.  sigma_theta_const will be based on AR of ODF.
  m_sigma_theta_const=sigma_theta_const;
  //now set iterations if they aren't set:


  m_binary_connectivity_map_on=true;

  if (m_MC_iterations==1 && m_sigma_theta_const> FLT_EPSILON) {
    m_MC_iterations=10000;//user needs to take care to change this afterward if less/more desired
  }

  m_cone=true;
}

/*
void FibreTracking::SetSigmaTheta(ImageData *sigma_theta_image)
{
  //this assumes already have ODF set.  sigma_theta_const will be based on AR of ODF.
  m_sigma_theta_image=sigma_theta_image;
  //now set iterations if they aren't set:


  m_binary_connectivity_map_on=true;

  if (m_MC_iterations==1)
    {
      m_MC_iterations=10000;
    }
  m_cone=true;
}
*/

SIGMA_THETA FibreTracking::GetSigmaTheta(int x, int y, int z)
{

  SIGMA_THETA sigma_theta;
  int n;



  /*
  if (m_sigma_theta_image!=NULL)
    {
      return m_sigma_theta_image->GetValue(x,y,z);
    }
  */




  if (m_tracking_algorithm==PDD) {

    //exists and calculated: return
    if (m_sigma_theta_image!=NULL) {
      if (m_sigma_theta_image->GetSigmaTheta(x,y,z).calculated) {
        return m_sigma_theta_image->GetSigmaTheta(x,y,z);
      }
    }
  } else { //not PDD

    if (m_ODF_data->SigmaThetaAvailable()) {
      if (m_ODF_data->GetSigmaTheta(x,y,z).calculated) {
        return m_ODF_data->GetSigmaTheta(x,y,z);
      }
    }
  }


  //if we have to assign, we assign m_sigma_theta_const: shouldn't really ever call this method if this needs to happen:

  sigma_theta.calculated =true;

  if (m_tracking_algorithm==PDD) {
    sigma_theta.number=1;
  } else { //not PDD
    sigma_theta.number=m_ODF_data->GetODFMaxima(x,y,z)->GetNumberOfDirections();
  }


  sigma_theta.sigma_theta=(float *) malloc(sigma_theta.number*sizeof(float));


  for (n=0;n<sigma_theta.number;n++) {
    sigma_theta.sigma_theta[n]=m_sigma_theta_const;

  }

  return sigma_theta;




}


void FibreTracking::SetReferenceROI(int reference_roi)
{
  m_reference_roi=reference_roi;
}

bool FibreTracking::CheckPolarity(VECTOR3D vector,PropagationInfo propagation_info)
{
  bool polarity_choice;

  if (propagation_info.start==true) { //this means there is no last_vector to consult, so randomize
    if ((float)((double)rand() / ((double)(RAND_MAX)+(double)(1.0)))<0.5) {
      polarity_choice=true;
    } else {
      polarity_choice=false;
    }
  } else if (dot(propagation_info.last_vector,vector)>=0) {
    polarity_choice=true;
  } else {
    polarity_choice=false;
  }

  return polarity_choice;

}


