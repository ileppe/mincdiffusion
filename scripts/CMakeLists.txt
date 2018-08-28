# copy the scripts to the INSTALL directory
INSTALL(PROGRAMS
 ascii_binary_vtk
 concat-diff-header.pl
 diff_preprocess.pl
 dti.m
 dti_octave.m
 duplicate_frame
 extract-diff-frames.pl
 flip_vector.m
 grad2FSL.pl
 header2FSL.pl 
 #mincFibreTrack.pl
 bvecs_bvals2mnc.pl 
 minc_modify_header-diff.pl
 mincbet
 mincresample_vector.pl
 minctensor.pl
 rgbvector
 DESTINATION bin )

