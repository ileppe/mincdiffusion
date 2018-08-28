/************************************************************************

   File: FibreTrack.c++

   Author: Jennifer Campbell

   Created:

   Revisions:

   Description: Fibretracking

   Copyright (c) 2005 by Jennifer Campbell, McConnell Brain Imaging
   Center and Centre for Intelligent Machines, McGill University,
   Montreal, QC.

***********************************************************************/
#include "jc_tools.h"
#include <string>
//#include "otherlibs_include.h"

/*
the tract scalars are currently output with -10 for the seed, inf for inf, and otherwise correct values.  minc3Dvis requires 0-1 to color-code them properly (left to user right now)


sigma_theta input as a minc file for PDD option may not be complete and has certanily never been tested

most of OTT/DOTT are not fully implemented.  DOTT I am leaving out for now (can add later and implement it given a *diffusion* ODF input for Maxima algorithm)

//resampling: use_input_sampling means the fibre tracking has the same step size as the input.  unless we resampling the diffusion data in a "smart" way, this very much blurs it, so might not be a good idea.
resample just means get all the inputs in the same voxel space, where they most probably aren't.  that can be done with -nearest for things that would otherwise be corrupted.  that's what's done with maxima (if input), eigenvectors...  masks and ROIs are currently resampled trilinear but watch out for missed small ones (?)

//note: the labels algo expects sigma theta

sigma_theta inputs are exclusive and there is currently no checking for this

for labels tracking, resampling isn't handled right now, as maxima must be left alone, and labels themselves better off untouched.  i.e., three "maxima" for fan must remain as such. must manually resample for labels algorithm first to match diffusion data

can amalgamate a bunch of these switches.

the labels algorithm uses the ODFMaxima (e.g., not done for tensor, but could by calculating the tensor ODF and maxima (which won't work that well because GetMaxima is not that great)).  Need to do more than one it to fan, but 1 it will select the ODF maximum only.  so one it is the same as Maxima algo.

for Cone and Labels, there is no uncertainty: will be added.

eventually, perhaps have an input ODF that incorporates uncert properly: so constant for the cone and then falloff.  Not sure: would require big input file (high N?)

right now, we only fan or single or cross in Labels algorithm: the rest can be done though

you can branch in  the Maxima algorithm by giving a file with 4 where you want to branch (nothing else matters)

many inputs required together: not checked for at this point

brute force if there are rois: I personally think multiple seeds per seed voxel is better than one seed in each voxel of the volume and checking whether they get back to the seed, however multiple seeds in each voxel of the volume is better still...at this point it is terrible slow

minDOTT: require the rest of the DT: e2,e3, lambdas;


for DOTT ODF (volume normalized), threshold at 0.6204 (gt)


added -resample option: if you are happy with resampling, don't redo.  must input this request else resamples even if they are already matching (it doesn't check - DO: could) note: that copying ODF files is time consuming

DO: add resampling of output minc files back to base (define how?) res afterward.  ie res of essential input of the algorithm (e1x, ODF, etc...)

this program assumes that diffusion data has the same direction cosines (and can therefore be resampled -like other input files without confounding the directions....)

DO: change -tmpdir input so user gives the base and the actual directory name is automatic: currently, user gives the whole thing


DO:change curvature_constraint* to curvature_constraint_radius here and in FibreTRacking, where not angle

OTT would in every case imaginable be fibre ODF

really need a new file format for the geometries and uncertainties that makes no assumption about the number of maxima

 */

typedef enum algorithm {FACTalgo, ODF, Maxima, Labelalgo} ALGORITHM;
//NOTE: I am not supporting the ODF algorithm any more, so have removed it from the input description, although the supporting code remains.  Probabilistic is really a new name in the help for the Labels option, the old non probabilistic labels being gone at this point

void Usage()
{

   std::cout << "Usage: FibreTrack [options]\n\n";

  std::cout << "algorithm options:\n\n-algorithm: PDD, Maxima.  The following summarizes these algorithm choices:\n\n";
  //Probabilistic.
  //Probabilistic: the tracking incorporates the output of ProbDeconvolve, which labels voxels as containing single fibres or double or triple crossings.  This is done in a probabilistic framework.  (Notes: (1) This option is non-probabilistic if the input file has no uncertainty in it.  (2) The code supports the output of a program called ProbLabels, which adds fanning to the possible voxel geometries, but this program is not yet part of the mincdiffusion tools.)\n\n"
  std::cout << "PDD: tracking follows the principal eigenvector (\"principal diffusion direction\") of the diffusion tensor - no quantification of uncertainty\n\nMaxima: tracking follows one of multiple fibre directions that have been obtained using a high angular resolution approach.  The fibre direction closest to the incoming fibre direction is selected.\n\n";
  std::cout << "PDD:\n";
  std::cout << "-e1x <e1x.mnc> -e1y <e1y.mnc> -e1z <e1z.mnc>\n";
  std::cout << "or: -namebase <namebase>\n\n";



  /*
  std::cout << "\nrequired for DOTT options (below):\n";
  std::cout << "-e2x <e2x.mnc> -e2y <e2y.mnc> -e2z <e2z.mnc>\n-e3x <e3x.mnc> -e3y <e3y.mnc> -e3z <e3z.mnc>\n-lambda1 <lambda1.mnc> -lambda2 <lambda2.mnc> -lambda3 <lambda3.mnc>";
  std::cout << "\nor: -namebase <namebase>\n";
  */
  //std::cout << "\nif ODF:\n";
  //std::cout << "-ODF <ODF_filename.mnc>\n\n";

  std::cout << "Maxima:\n";
  std::cout << "-ODF <ODF_filename.mnc>\n";
  //std::cout << "-ODF_maxima <maxima_filename>: optionally input these if already calculated: can only input them if use -noresample option\n";
  //std::cout << "-sigma_theta <sigma_theta_filename> (optional: input these if already calculated: can only input them if use -noresample.)\n";

  //the labels algorithm requires a maxima file that is maxima except in fanning case, there are three vectors: the first and third give the range of the fan and the second the merge direction.

  //std::cout << "\n-labels [<labels.mnc>]:  For the Maxima algorithm, this is a 3D minc file constructed as follows: integer value 4: branch, otherwise will follow nearest to incoming direction. This branching option with the Maxima algorithm is not recommended.
  /*
  std::cout << "\nProbabilistic:\n";

  std::cout << "\n-labels <labels>: the (6) files must be of the format output by ProbDeconvolve: see ProbDeconvolve -help for details.  Input only the base name\n";
  */

  //DO:  think about sigma_theta stuff: probably remove it entirely

  std::cout << "\n\noutput files and options:\n";

  std::cout << "\n-tracts <output_tracts.vtk>.  This option is discouraged for the PDD and Maxima algorithms if doing extensive brute force tractography (see below for the brute_force option).\n";
  //for the Probabilistic algorithm, and 
  std::cout << "-conn <binary_connectivity_map.mnc>: 1 where tracking went, 0 where it didn't go\n";
  //std::cout << "-conn_index <connectivity_index_filename.mnc>: (for Probabilistic algorithm): at each voxel tract passes through, write the connectivity index for the tract segment between the reference ROI (see above) and that voxel.\n";

  std::cout << "\n-binary: save vtk file in binary format (default ascii)\n";
  std::cout << "-obj: output tracts in MINC .obj format rather than vtk format\n";
  std::cout << "-nobbox: save full-size output minc files [default: smallest bounding box]\n";
  std::cout << "-clobber: overwrite existing output file\n";



  std::cout << "\nother options for all algorithms:\n\n";
  std::cout << "-brute_force: start at all voxels and check whether went through ROIs\n";
  std::cout << "-mask <mask.mnc> -(gt,lt,ge,le) <threshold> (only one comparison allowed per mask)\n";
  std::cout << "-roi <roi.mnc> (all values greater than 0 are treated is being inside the roi).  multiple -roi inputs may be given: you must type -roi before all input rois, one at a time; tracts that pass through ALL ROIs are retained\n";
  //std::cout << "-reference_roi <integer>: integer corresponding to the -roi in the input list that connectivity index values (cumulative minimum ODF tangent to tracts; see below) are relative to.  default: 1.  1 is the first. no 0\n";
  std::cout << "-exclude <exclusion_mask.mnc>: remove all tracts that pass through this ROI\n";
  //std::cout << "-curvature_constraint_constant: (radius of curvature - mm)\n";
  //std::cout << "-curvature_constraint <curv_constr.mnc> (radius of curvature - mm)\n";
  std::cout << "-curvature_constraint_angle <angle>: angle in degrees, applied for each step regardless of -step input\n";
  std::cout << "-roi_t <roi.txt> -(v,w) : text file listing coordinates of roi (x1 y1 z1\\n ... xn yn zn\\n) in either voxel(v) or world (w) coordinates\n";
  std::cout << "-max_length <length>: maximum tract length in mm\n";

  //for DOTT, need to implement this feature for Maxima (given *diffusion* ODF) and check Pruning, FibreTracking, etc., to make sure it all works. applications not clear, and could perhaps just do AI instead of DOTT.
  /*
  std::cout << "-min_DOTT <min_DOTT.mnc> (requires full diffusion ODF; implemented only for -FACT)\n";
  std::cout << "-cumulative_min_DOTT: write cumulative minimum diffusion ODF value tangent to tract at each point of all tracts (requires full diffusion ODF; implemented only for -FACT)\n";
  */

  //std::cout << "-OTT <OTT.mnc> (for ODF algorithm, and for Maxima and Labels algorithmms IF sigma_theta (see below) > 0): at each voxel tract passes through, write ODF value tangent to the tract (OTT).  If multiple tracts pass through a voxel, the maximal OTT is assigned. If the -tracts option is given, the tracts will have point scalars for OTT\n";
  //std::cout << "-min_OTT <min_OTT.mnc> (for ODF algorithm, and for Maxima and Labels algorithmms IF sigma_theta (see below) > 0): at each voxel tract passes through, write minimum ODF value tangent to the tract (minOTT).  If multiple tracts pass through a voxel, the maximal minOTT is assigned. If the -tracts option is given, the tracts will have scalars for minOTT\n";


  //for cumulative_min_OTT, if the -tracts option is also given, the tracts will have point scalars for cumulative minOTT: it makes no sense to do this because the tracts are pruned (DO: could not prune in this case, but be careful about memory: find limits!!  system specific), and it doesn't make sense to eliminate any.  (if it is used, there is no error checking: the tracts output will be the longest with arbitrary cum_min_OTT profile along tract.)
  std::cout << "-step <step>: grid size for tracking.  default is the voxel sampling, and this is generally recommended for speed when using brute_force option\n";
  std::cout << "-tmpdir: a temporary directory name (to be created by the program).  Default is /tmp/tmpXXXXXXXXX\n";
  std::cout << "-save_temp_dir: save the temporary directory and resampled files (so can run without resampling in the future if desired)\n";
  std::cout << "-noresample: don't resample any input files: assume in same space\n";
  //std::cout << "-sigma_theta <file.mnc>: input file with sigma_theta (in radians) for ODF maximum/maxima at each voxel\n";//not doing this because need unique sigma theta for each maximum: see above; use txt file
  //std::cout << "-sigma_theta_AR: use angular resolution of acquisition as only source of uncertainty\n";
  //std::cout << "-sigma_theta_const <const>: input a constant sigma_theta, in radians\n";

  //std::cout << "-iterations: iterations of algorithm (default: 10000) (Probabilistic algorithm)\n";


  std::cout << "-seedfreq <freq>: number of starts per voxel per dimension (default: 3)\n ";

  exit(0);
}

int main(int argc, char *argv[])
{

  int i;

  bool clobber=false;
  std::string output_tracts_name;

  std::string e1x_filename;
  std::string e1y_filename;
  std::string e1z_filename;

  std::string e2x_filename;
  std::string e2y_filename;
  std::string e2z_filename;
  std::string e3x_filename;
  std::string e3y_filename;
  std::string e3z_filename;
  std::string lambda1_filename;
  std::string lambda2_filename;
  std::string lambda3_filename;

  std::string namebase;

  std::string mask_filename_array[100];
  int number_of_masks=0;
  Comparison comparison_array[100];
  float threshold_array[100];


  std::string roi_filename_array[100];
  int number_of_rois=0;
  int reference_roi=1;

  std::string roi_t_filename_array[100];//text roi array

  int number_of_rois_t=0;
  bool use_text_roi = false;
  bool voxel=true; /*default input roi in voxel coordinates*/
  std::string input_tracts_name;

  bool max_tract_length_exists=false; //if the user wants to
  float max_tract_length; //length in mm

  std::string exclusion_mask_filename_array[100];
  int number_of_exclusion_masks=0;


  bool curvature_constraint_exists=false; //if there is one
  bool variable_curvature_constraint; //if it is variable, else it is const float
  std::string curvature_constraint_filename;
  float curvature_constraint;
  float curvature_constraint_angle;
  bool curvature_constraint_angle_exists=false;

  ALGORITHM algorithm;

  bool writeBCM=false;
  std::string BCM_filename;
  //bool write_min_DOTT=false;
  //std::string min_DOTT_filename;
  //bool write_cumulative_min_DOTT=false;

  ProcessedDTIData *processed_DTI_data;

  ImageData *ODF_data;
  std::string ODF_filename;
  //bool write_min_OTT=false;
  std::string cum_min_OTT_filename;
  bool write_cumulative_min_OTT=false;

  int iterations=10000;//this doesn't get set in FibreTracking unless iterations would be called for.  So for FACT & Multimax with no uncetainty it should never get set although will be 10000 here

  bool save_tracts=false;

  bool brute_force=false;

  bool ASCII=true; //IL default output VTK format is ascii
  bool OBJ_FORMAT=false; //VF by default save in VTK format

  std::string ODF_maxima_filename;


  bool use_input_sampling=true;
  float step;
  //std::string templatefilename;
  //std::string templatefilename=(std::string ) malloc(200*sizeof(char));
  std::string templatefilename;
  
  bool templatefilename_set=false;
  std::string templatefilename_for_roi;

  //char tempdirname[200];
  std::string tempdirname;
  char command[1024];



  std::string testfilename;
  std::string testfilename2;

  bool labels=false;
  bool labels_exist=false;
  std::string labels_filename;

  ImageData *prob_labels0, *prob_labels1, *prob_labels2, *prob_labels3, *prob_labels4, *prob_labels5;

  bool save_tmp_dir=false;

  FibreTracking *fibre_tracking;

  bool resample=true;

  bool bbox=true; // IL by default save output minc files as smallest possbile volume

  bool cone_AR=false;
  bool cone_const=false;
  bool cone=false;

  float sigma_theta_const;
  std::string sigma_theta_name;


  int seedfreq=3;

  //if type FibreTrack with no args, give help

  if (argc == 1) {
    Usage();
  }

  try { //VF

    //set up temporary directory for resampled files:  if tmpdir is input, this pointer gets overwritten
    char tmp[1024];
    sprintf(tmp,"/tmp/tmp%i",(int) time(NULL)); //VF: bad idea
    tempdirname=tmp;

    //get the inputs
    //note that some of these inputs are not mentioned in the help, but are here for historical/development/experimental purposes

    for (i = 1; i < argc; i++) {
      if (!strcmp(argv[i], "-help")) {
        Usage();
      } else if (!strcmp(argv[i], "-clobber")) {
        clobber = true;
      } else if (!strcmp(argv[i], "-save_temp_dir")) {
        save_tmp_dir = true;
      } else if (!strcmp(argv[i], "-noresample")) {
        resample = false;
      }

      // IL option to save smallest bounding box
      else if (!strcmp(argv[i], "-nobbox")) {
        bbox = false;
      } else if (!strcmp(argv[i], "-brute_force")) {
        brute_force = true;
      }

      // IL binary option
      else if (!strcmp(argv[i], "-binary")) {
        ASCII = false;
      } else if (!strcmp(argv[i],"-obj")) {
        OBJ_FORMAT = true;  //VF
      } else if (!strcmp(argv[i], "-labels")) {
        labels = true;

        if (i+1<argc) {
          if (argv[i+1][0] != '-') {
            labels_filename=argv[++i];
            labels_exist=true;
          }
        }
      } else if (!strcmp(argv[i], "-mask")) {
        if (++i >= argc) Usage();

        mask_filename_array[number_of_masks] = argv[i++];

        if (!strcmp(argv[i], "-gt")) {
          comparison_array[number_of_masks]=gt;
        }

        else if (!strcmp(argv[i], "-ge")) {
          comparison_array[number_of_masks]=ge;
        } else if (!strcmp(argv[i], "-le")) {
          comparison_array[number_of_masks]=le;
        } else if (!strcmp(argv[i], "-lt")) {
          comparison_array[number_of_masks]=lt;
        }

        if (++i >= argc) Usage();

        sscanf(argv[i], "%f", &threshold_array[number_of_masks]);

        number_of_masks+=1;

      }

      else if (!strcmp(argv[i], "-roi")) 
      {
        if (++i >= argc) Usage();

        roi_filename_array[number_of_rois]=argv[i];

        number_of_rois+=1;
      }

      else if (!strcmp(argv[i], "-roi_t")) 
      {
        if (++i >= argc) Usage();

        roi_t_filename_array[number_of_rois_t]=argv[i++];

        if (!strcmp(argv[i], "-w")) {
          voxel=false;
        }

        else if (!strcmp(argv[i], "-v")) {
          voxel=true;
        }

        number_of_rois_t+=1;

        use_text_roi = true;

      } else if (!strcmp(argv[i], "-reference_roi")) {
        if (++i >= argc) Usage();

        sscanf(argv[i], "%i", &reference_roi);
      } else if (!strcmp(argv[i], "-exclude")) {
        if (++i >= argc) Usage();

        exclusion_mask_filename_array[number_of_exclusion_masks]=argv[i];

        number_of_exclusion_masks+=1;
      } else if (!strcmp(argv[i], "-curvature_constraint_constant")) {
        if (++i >= argc) Usage();

        curvature_constraint_exists=true;

        variable_curvature_constraint=false;

        sscanf(argv[i], "%f", &curvature_constraint);

      } else if (!strcmp(argv[i], "-curvature_constraint_angle")) {
        if (++i >= argc) Usage();

        curvature_constraint_angle_exists=true;

        variable_curvature_constraint=false;

        sscanf(argv[i], "%f", &curvature_constraint_angle);

      } else if (!strcmp(argv[i], "-iterations")) {
        if (++i >= argc) Usage();

        sscanf(argv[i], "%i", &iterations);

      }


      else if (!strcmp(argv[i], "-seedfreq")) {
        if (++i >= argc) Usage();

        sscanf(argv[i], "%i", &seedfreq);

      } else if (!strcmp(argv[i], "-curvature_constraint")) {
        if (++i >= argc) Usage();

        curvature_constraint_exists=true;

        variable_curvature_constraint=true;

        curvature_constraint_filename=argv[i];
      }

      else if (!strcmp(argv[i], "-algorithm")) {
        if (++i >= argc) Usage();

        if (!strcmp(argv[i], "FACT")) {
          algorithm=FACTalgo;
	  //alt name:
	} else if (!strcmp(argv[i], "PDD")) {
	  algorithm=FACTalgo;
        } else if (!strcmp(argv[i], "ODF")) {
          algorithm=ODF;
        } else if (!strcmp(argv[i], "Maxima")) {
          algorithm=Maxima;
        } else if (!strcmp(argv[i], "Labels")) {
          algorithm=Labelalgo;
          labels=true;
        }
        //same as Labels, different name:
        else if (!strcmp(argv[i], "Probabilistic")) {
          algorithm=Labelalgo;
          labels=true;
  	}
      }

      // IL segmentation fault when usind this option due to lack of mem allocation. Done here for all 'e's and 'lambda's (could maybe be smaller than 200 chars?)
      else if (!strcmp(argv[i], "-e1x")) {
        if (++i >= argc) Usage();


        e1x_filename=argv[i];

        //e1x_filename = argv[i];
      } else if (!strcmp(argv[i], "-e1y")) {
        if (++i >= argc) Usage();

        e1y_filename=argv[i];

        //e1y_filename = argv[i];
      } else if (!strcmp(argv[i], "-e1z")) {
        if (++i >= argc) Usage();

        e1z_filename=argv[i];

        //e1z_filename = argv[i];
      } else if (!strcmp(argv[i], "-e2x")) {
        if (++i >= argc) Usage();

        e2x_filename=argv[i];

        //e2x_filename = argv[i];
      } else if (!strcmp(argv[i], "-e2y")) {
        if (++i >= argc) Usage();

        e2y_filename=argv[i];

        //e2y_filename = argv[i];
      } else if (!strcmp(argv[i], "-e2z")) {
        if (++i >= argc) Usage();

        e2z_filename=argv[i];

        //e2z_filename = argv[i];
      } else if (!strcmp(argv[i], "-e3x")) {
        if (++i >= argc) Usage();

        e3x_filename=argv[i];

        //e3x_filename = argv[i];
      } else if (!strcmp(argv[i], "-e3y")) {
        if (++i >= argc) Usage();

        e3y_filename=argv[i];

        //e3y_filename = argv[i];
      } else if (!strcmp(argv[i], "-e3z")) {
        if (++i >= argc) Usage();

        e3z_filename=argv[i];

        //e3z_filename = argv[i];
      } else if (!strcmp(argv[i], "-lambda1")) {
        if (++i >= argc) Usage();

        lambda1_filename=argv[i];

        //lambda1_filename = argv[i];
      } else if (!strcmp(argv[i], "-lambda2")) {
        if (++i >= argc) Usage();

        lambda2_filename=argv[i];

        //lambda2_filename = argv[i];
      } else if (!strcmp(argv[i], "-lambda3")) {
        if (++i >= argc) Usage();

        lambda3_filename=argv[i];

        //lambda3_filename = argv[i];
      } else if (!strcmp(argv[i], "-namebase")) {
        if (++i >= argc) Usage();

        namebase = argv[i];

        //and create the other strings:
        e1x_filename=namebase;

        e1x_filename+="_e1x.mnc";

        e1y_filename=namebase;

        e1y_filename+="_e1y.mnc";

        e1z_filename=namebase;

        e1z_filename+="_e1z.mnc";

        e2x_filename=namebase;

        e2x_filename+="_e2x.mnc";

        e2y_filename=namebase;

        e2y_filename+="_e2y.mnc";

        e2z_filename=namebase;

        e2z_filename+="_e2z.mnc";

        e3x_filename=namebase;

        e3x_filename+="_e3x.mnc";

        e3y_filename=namebase;

        e3y_filename+="_e3y.mnc";

        e3z_filename=namebase;

        e3z_filename+="_e3z.mnc";

        lambda1_filename=namebase;

        lambda1_filename+="_lambda1.mnc";

        lambda2_filename=namebase;

        lambda2_filename+="_lambda2.mnc";

        lambda3_filename=namebase;

        lambda3_filename+="_lambda3.mnc";

      }

      else if (!strcmp(argv[i], "-tracts")) {
        if (++i >= argc) Usage();

        output_tracts_name = argv[i];

        save_tracts=true;

      } else if (!strcmp(argv[i], "-max_length")) { /*IL*/
        if (++i >= argc) Usage();

        max_tract_length_exists=true;

        sscanf(argv[i], "%f", &max_tract_length);

      } else if (!strcmp(argv[i], "-conn")) {
        if (++i >= argc) Usage();

        BCM_filename = argv[i];

        writeBCM=true;
      }

      /*
      else if (!strcmp(argv[i], "-min_DOTT"))
	{
	  if(++i >= argc) Usage();
	  min_DOTT_filename = argv[i];
	  write_min_DOTT=true;

	}
      else if(!strcmp(argv[i], "-cumulative_min_DOTT"))
	{
	  write_cumulative_min_DOTT = true;
	}
      else if (!strcmp(argv[i], "-min_OTT"))
	{
	  if(++i >= argc) Usage();
	  min_OTT_filename = argv[i];
	  write_min_OTT=true;

	}
      */
      else if(!strcmp(argv[i], "-cumulative_min_OTT"))
	{
	  if(++i >= argc) Usage();
	  cum_min_OTT_filename = argv[i];
	  write_cumulative_min_OTT = true;
	}
      else if(!strcmp(argv[i], "-conn_index"))
	{
	  if(++i >= argc) Usage();
	  cum_min_OTT_filename = argv[i];
	  write_cumulative_min_OTT = true;
	}
      else if (!strcmp(argv[i], "-ODF"))
	{
	  if(++i >= argc) Usage();
	  ODF_filename = argv[i];

	}
      else if (!strcmp(argv[i], "-ODF_maxima"))
	{
	  if(++i >= argc) Usage();
	  ODF_maxima_filename = argv[i];

	}
      else if (!strcmp(argv[i], "-step"))
	{
	  if(++i >= argc) Usage();
	  sscanf(argv[i], "%f", &step);
	  use_input_sampling = false;

	}
      else if (!strcmp(argv[i], "-tmpdir"))
	{
	  if(++i >= argc) Usage();
	  tempdirname = argv[i];

	}
      else if (!strcmp(argv[i], "-sigma_theta"))
	{
	  if(++i >= argc) Usage();
	  sigma_theta_name = argv[i];
	  cone=true;

	}
      else if (!strcmp(argv[i], "-sigma_theta_AR"))
	{
	  cone_AR=true;
	  cone=true;
	}
      else if (!strcmp(argv[i], "-sigma_theta_const"))
	{
	  if(++i >= argc) Usage();
	  sscanf(argv[i], "%f", &sigma_theta_const);
	  cone_const=true;
	  cone=true;
	}

      else if(argv[i][0] == '-')
	{
	  Usage();
	}

    }


    sprintf(command,"mkdir %s\n",tempdirname.c_str());

    //cout << command;
    system(command);

    fibre_tracking= new FibreTracking();


    //DO: check for required inputs here.
    switch (algorithm) {
    case FACTalgo:

      if (resample) {
        // IL: segmentation fault when not using 'namebase' option: can't realloc without malloc first, did this in arg parsing above
        Resample(e1x_filename,e1x_filename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,true);
        Resample(e1y_filename,e1y_filename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,true);
        Resample(e1z_filename,e1z_filename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,true);
      }

      /*not doing this right now: scalar file needs to be cumulative and relative to something and it currently isn't
      if (write_cumulative_min_DOTT || write_min_DOTT)
      {

      if (resample)
      {

       e2x_filename=(std::string ) realloc(e2x_filename,200*sizeof(char));
       Resample(e2x_filename,e2x_filename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,true);
       e2y_filename=(std::string ) realloc(e2y_filename,200*sizeof(char));
       Resample(e2y_filename,e2y_filename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,true);
       e2z_filename=(std::string ) realloc(e2z_filename,200*sizeof(char));
       Resample(e2z_filename,e2z_filename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,true);
       e3x_filename=(std::string ) realloc(e3x_filename,200*sizeof(char));
       Resample(e3x_filename,e3x_filename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,true);
       e3y_filename=(std::string ) realloc(e3y_filename,200*sizeof(char));
       Resample(e3y_filename,e3y_filename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,true);
       e3z_filename=(std::string ) realloc(e3z_filename,200*sizeof(char));
       Resample(e3z_filename,e3z_filename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,true);
       lambda1_filename=(std::string ) realloc(lambda1_filename,200*sizeof(char));
       Resample(lambda1_filename,lambda1_filename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,true);
       lambda2_filename=(std::string ) realloc(lambda2_filename,200*sizeof(char));
       Resample(lambda2_filename,lambda2_filename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,true);
       lambda3_filename=(std::string ) realloc(lambda3_filename,200*sizeof(char));
       Resample(lambda3_filename,lambda3_filename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,true);
      }


      processed_DTI_data = new ProcessedDTIData(e1x_filename,e1y_filename,e1z_filename, e2x_filename,e2y_filename,e2z_filename, e3x_filename,e3y_filename,e3z_filename, lambda1_filename, lambda2_filename, lambda3_filename);
      }
      else
      */
      //{

      processed_DTI_data = new ProcessedDTIData(e1x_filename.c_str(),e1y_filename.c_str(),e1z_filename.c_str());

      if (use_text_roi) {
        templatefilename_for_roi =  e1x_filename;
      }

      //}

      break;

    case ODF:
      if (resample) {
        Resample(ODF_filename,testfilename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,true);
      } else {
        testfilename=ODF_filename;
      }

      ODF_data=new ImageData(testfilename.c_str());

      if (use_text_roi) {
        templatefilename_for_roi =  testfilename;
      }

      break;

    case Maxima:

      if (resample) {
        Resample(ODF_filename,testfilename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,true);
      } else {
        testfilename=ODF_filename;
      }

      ODF_data=new ImageData(testfilename.c_str());

      if (use_text_roi) {
        templatefilename_for_roi =  testfilename;
      }

      break;

    case Labelalgo:

      if (resample) {
        Error("all files must be sampled the same way at input: use -noresample option.");
      }

      //no longer requires ODF: just ProbLabels minc file, input later (below)
      //note no resampling of ProbLabels minc file as code is now.  it doesn't make sense to do so

      /*
      //should not resample (as our resampling is now) and should always input maxima (in the case of fanning, these will define limits of fan, so incorrect to recalculate).  this is because the labels and the maxima HAVE to match, and resmapling might mess this up
      if (ODF_maxima_filename==NULL)
      {
       //hardcoded for problabels: this is not used
       //Error("ODF Maxima are required.");
      }
      if (resample)
	{
	  cout << "NOTE: should never resample in this case\n";
	  Resample(ODF_filename,testfilename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,false);

	}
      else
      {
      strcpy(testfilename,ODF_filename);
      }
      ODF_data=new ImageData(testfilename);
      */
      if (use_text_roi) {
        templatefilename_for_roi =  testfilename;
      }

      break;
    }
    char * _history=CreateHistoryString(argc, argv);
    std::string history=_history;
    free(_history);
    //cout << "history " << history << endl;

    if (!ODF_maxima_filename.empty()) {
      if (!resample) {
        //ODF maxima will only be used if no resampling is being done.  If it is, we will recalculate
        ODF_data->SetODFMaxima(ODF_maxima_filename.c_str());
      }
    }



    //if sigma_theta exists, give it to the ODF
    if (!sigma_theta_name.empty()) {
      if (!resample) {
        ODF_data->SetSigmaTheta(sigma_theta_name.c_str());
      } else {
        Error("cannot resample if inputting sigma_theta file");
      }
    }

    //set up the FibreTracking:

    if ((algorithm==Maxima || algorithm==Labelalgo)&& labels) {
      if (labels_exist) {
        if (resample) {
          Resample(labels_filename,testfilename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,false);
        } else {
          testfilename=labels_filename;
        }


      }

    }

    //first the type:
    switch (algorithm) {
    case FACTalgo:
      fibre_tracking->SetTrackingAlgorithm(PDD);
      fibre_tracking->SetProcessedDTIData(processed_DTI_data);
      break;
    case ODF:
      fibre_tracking->SetTrackingAlgorithm(ODF_MC);
      fibre_tracking->SetODFData(ODF_data);
      break;
    case Labelalgo:
      //probabilistic label tracking: no option for non-probabilistic right now

      if (1) {
        fibre_tracking->SetTrackingAlgorithm(Labels);
      
	//create the full filenames:
	for (i=0; i<6; i++)
	{
	  char inttostring[5];
	  char OutFile[200];
		memset(OutFile, '\0', 200);
		strcpy(OutFile, labels_filename.c_str());
		sprintf(inttostring,"%i",i);
		strcat(OutFile, inttostring);
		strcat(OutFile, ".mnc");
		switch (i)
		  {
		  case 0:
		    prob_labels0=new ImageData(OutFile);
		    break;
		  case 1:
		    prob_labels1=new ImageData(OutFile);
		    break;
		  case 2:
		    prob_labels2=new ImageData(OutFile);
		    break;
		  case 3:
		    prob_labels3=new ImageData(OutFile);
		    break;
		  case 4:
		    prob_labels4=new ImageData(OutFile);
		    break;
		  case 5:
		    prob_labels5=new ImageData(OutFile);
		    break;
		  }
	}

        fibre_tracking->SetProbLabels(prob_labels0, prob_labels1, prob_labels2, prob_labels3, prob_labels4, prob_labels5);

      } else { //not used
        fibre_tracking->SetTrackingAlgorithm(Labels);
        fibre_tracking->SetODFData(ODF_data);
      }

      break;

    case Maxima:
      fibre_tracking->SetTrackingAlgorithm(MultiMax);
      fibre_tracking->SetODFData(ODF_data);
      break;
    }



    //now some general options:

    fibre_tracking->SetIntegrationType(FACT);//only option right now

    if (brute_force) {
      fibre_tracking->SetBruteForceOn();
    }

    if (number_of_rois<1 && number_of_rois_t<1) { /*add text roi option*/
      Error("must input at least one roi");

      //DO: this only if requested: gives every tract in the masks!
      //fibre_tracking->SetBruceForceOn();
    }


    fibre_tracking->SetSeedFrequency(seedfreq);

    if (cone || algorithm ==ODF || algorithm ==Labelalgo) {

      fibre_tracking->SetMCIterations(iterations);
    }


    for (i=0;i<number_of_masks;i++) {
      if (resample) {
        Resample(mask_filename_array[i],testfilename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,false);
      } else {
        testfilename=mask_filename_array[i];
      }
      fibre_tracking->SetMask(testfilename.c_str(),comparison_array[i],threshold_array[i]);
    }

    //VF foolproofing 
    if (number_of_masks>0)//JC foolproofing the foolproofing
      {
    if(fibre_tracking->GetBinaryTrackingMask()->GetNumberOfVoxelsAbove(0.5)==0.0)
    {
      Error("Your ROI doesn't have any voxels set, check your masks and -gt -lt options !");
    }
      }

    for (i=0;i<number_of_rois;i++) {

      if (resample) {
        Resample(roi_filename_array[i],testfilename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,false);
      } else {
        testfilename=roi_filename_array[i];
      }

      fibre_tracking->SetROI(testfilename.c_str());
    }

    /*rois entered as text*/
    for (i=0;i<number_of_rois_t;i++) {
      /*fibre_tracking->ReadROI(roi_t_filename_array[i].c_str(), 
                              testfilename2.c_str(), tempdirname.c_str(), 
                              templatefilename_for_roi.c_str(), voxel);     */
        fibre_tracking->ReadROI(roi_t_filename_array[i], 
                            testfilename2, tempdirname, 
                            templatefilename_for_roi, voxel); 
      /*read roi txt file into mnc file*/

      if (resample) {
        Resample(testfilename2,testfilename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,false);
      } else {
        testfilename=testfilename2;
      }

      fibre_tracking->SetROI(testfilename.c_str());
    }

    fibre_tracking->SetReferenceROI(reference_roi-1);

    for (i=0;i<number_of_exclusion_masks;i++) {
      if (resample) {
        Resample(exclusion_mask_filename_array[i],testfilename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,false);
      } else {
        testfilename=exclusion_mask_filename_array[i];
      }

      fibre_tracking->SetExclusionMask(testfilename.c_str());
    }




    if (writeBCM) {
      fibre_tracking->BinaryConnectivityMapOn();

    }

    //if (write_min_DOTT || write_cumulative_min_DOTT)
    //{
    //  fibre_tracking->MinDOTTMapOn();
    //}

    if (write_cumulative_min_OTT) {
      fibre_tracking->CumMinOTTMapOn();
    }



    if (curvature_constraint_exists || curvature_constraint_angle_exists) {
      if (curvature_constraint_angle_exists) {
        fibre_tracking->SetCurvatureConstraintAngle(curvature_constraint_angle);
      } else if (!variable_curvature_constraint && curvature_constraint_exists) {
        fibre_tracking->SetCurvatureConstraint(curvature_constraint);

      } else {

        if (resample) {
          Resample(curvature_constraint_filename,testfilename,tempdirname,step,use_input_sampling,templatefilename,templatefilename_set,false);
        } else {
          testfilename=curvature_constraint_filename;
        }
        fibre_tracking->SetCurvatureConstraint(testfilename.c_str());
      }

    }

    if (max_tract_length_exists) { /*IL*/
      fibre_tracking->SetMaxTractLength(max_tract_length);
    }





    if (save_tracts) {

      // IL uncomment this to test pruning
      //Pruning OK now for FACT
      //this isn't default because theoretically one might want all of the tracts even if same indices: i.e., different scalars.

      fibre_tracking->SetPruneOn();

      if (iterations>10 && (cone || algorithm ==ODF || algorithm ==Labelalgo)) { //DO: check where memory issues begin (with bf and no)
        cout << "warning: memory overflow possible\n";
      }
    } else {
      fibre_tracking->GetFibreTractSet()->TractsOff();

    }
    if (cone_AR) {
      fibre_tracking->SetSigmaTheta();
    }

    if (cone_const) {
      fibre_tracking->SetSigmaTheta(sigma_theta_const);
    }


    fibre_tracking->Track();

    //DO: for all output files, resample -like original diffusion data: right now, if !use_input_sampling, output will be with the new step

    //there is always a FibreTractSet. Scalar maps are in it even if !save_tracts
    FibreTractSet *fibre_tract_set=fibre_tracking->GetFibreTractSet();


    //writing history here because seems to be failing in ImageData:
    char mycommand[1000];

    //cout << history.c_str() << endl;

    if (writeBCM) {

      fibre_tract_set->GetBinaryConnectivityMap()->WriteMincFile(BCM_filename.c_str(),clobber,history.c_str(),bbox);// IL changed this to include bbox option
     
      sprintf(mycommand,"minc_modify_header -sinsert :history=\"%s\" %s\n",history.c_str(),BCM_filename.c_str());
      //cout << mycommand << endl;
      system(mycommand);// happens for PDD, Probabilistic, Maxima. 

    }

    if (write_cumulative_min_OTT) {

      fibre_tract_set->SetPointScalarOutputMode(cumulative_min_OTT);
      fibre_tract_set->GetCumMinOTTMap()->WriteMincFile(cum_min_OTT_filename.c_str(),clobber,history.c_str(),bbox);

      sprintf(mycommand,"minc_modify_header -sinsert :history=\"%s\" %s\n",history.c_str(),cum_min_OTT_filename.c_str());

      system(mycommand);
     

    }

    if (save_tracts) {
      FibreTractSet *fibre_tract_set=fibre_tracking->GetFibreTractSet();


      /* got rid of these options because they all have theoretical problems
      if (write_min_DOTT)
      {

      fibre_tract_set->GetMinDOTTMap()->WriteMincFile(min_DOTT_filename,clobber,history,bbox);// IL changed this to include bbox option
      }

      if (write_cumulative_min_DOTT)+
      {
      fibre_tract_set->SetPointScalarOutputMode(cumulative_min_DOTT);
      fibre_tract_set->GetCumMinDOTTMap()->WriteMincFile(min_DOTT_filename,clobber,history,bbox);

      }

      if (write_min_OTT)
      {
      fibre_tract_set->GetMinOTTMap()->WriteMincFile(min_OTT_filename,clobber,history,bbox);// IL changed this to include bbox option
      }
      */





      if (ASCII) { // IL added for ASCII option

        fibre_tract_set->OutputVTKAscii();
      }

      if (OBJ_FORMAT)
        fibre_tract_set->OutputOBJFile(output_tracts_name.c_str(),clobber); //VF
      else
        fibre_tract_set->OutputVTKFile(output_tracts_name.c_str(),clobber);

    }



    }
  
  catch (...) { //VF
    printf("Program unexpectedly crashed!\n");
  }
  

  if (!save_tmp_dir) {
    char command[1024];
    sprintf(command,"rm -f %s/tmp*\n",tempdirname.c_str());
    //cout << command << endl;
    system(command);
    sprintf(command,"rmdir %s\n",tempdirname.c_str());
    //cout << command << endl;
    system(command);
  }


  return(0);

}
