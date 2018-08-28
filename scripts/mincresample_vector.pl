#! /usr/bin/env perl

### Applies a transformation to vector data
### 
##########################################################
##########################################################
## usage: mincresample_vector.pl -xfm <transform> -vector <x y z>
##
##  Inputs: - xfm : transformation file
##    	    - vector : input files (should be the 3 seprarate components of a vector)
##	    - target : name of target volume
##  Output: resampled transformed volumes with  "norm_res_tx" appended before filename
##			NOTE: automatically clobbers
###########################################
#### Created by Ilana Leppert
####  Jan 2007
####

require 5.001;
use Getopt::Tabular;
use MNI::Startup qw(nocputimes);
use MNI::Spawn;
use MNI::FileUtilities qw(check_output_dirs);
use File::Basename;
#use ileppe_lib;

if($0 =~ /[\/A-Za-z_0-9]+\/([A-Za-z0-9_]+)\.pl/) {$program = $1;}	#program name without path
elsif($0 =~ /([A-Za-z0-9_-]+\.pl)/) {$program = $1;}
$Usage = <<USAGE;

Usage: $program -xfm <transform.xfm> -target <target.mnc> -vector <x.mnc y.mnc z.mnc> [-mask <mask.mnc>]
-help for options

USAGE

@args_table = (["-xfm","string",1,\$xfm,"transformation file","<tx.xfm>"],
	       ["-target","string",1,\$target,"target volume name","<target.mnc>"],
	       ["-vector","string",3,\@vector,"3 vector components","<x.mnc> <y.mnc> <z.mnc>"],
	       ["-mask","string",1,\$maskname,"mask file: tensor will not be calculated where mask file value is 0","<mask.mnc>"],
);

$maskname="none";

Getopt::Tabular::SetHelp ($Usage, '');

GetOptions(\@args_table, \@ARGV) || exit 1;


## Vector components
($x,$dir,$ext) = fileparse($vector[0],'\..*');
($y,$dir,$ext) = fileparse($vector[1],'\..*');
($z,$dir,$ext) = fileparse($vector[2],'\..*');

if ($maskname ne "none") 
{
    system("mincmath -clobber -nocheck_dimensions -copy_header -mult $vector[0] $maskname $x-mask.mnc");
    system("mincmath -clobber -nocheck_dimensions -copy_header -mult $vector[1] $maskname $y-mask.mnc");
    system("mincmath -clobber -nocheck_dimensions -copy_header -mult $vector[2] $maskname $z-mask.mnc");
    @input=($x.'-mask.mnc',$y.'-mask.mnc',$z.'-mask.mnc');
}
else {
@input=@vector;
}

#print "$x \t $y \t $z\n\n";
## Get params for rotation from xfm file
## NOTE: only rotation will modify the relative intensity of elements of the vector
$param = `xfm2param $xfm`;
#print "param: $param\n";
if($param =~ /rotation\s+([-]?\d+\.\d+)\s+([-]?\d+\.\d+)\s+([-]?\d+\.\d+)\s*/){$xrot = $1; $yrot = $2; $zrot = $3;}else{die "Can't find rotation from xfm file\n";}

#print "$xrot\t$yrot\t$zrot\n";

## create xfm file with rotation only
`param2xfm -clobber -rotations $xrot $yrot $zrot tmprot.xfm`;

open(ROT,'tmprot.xfm');
@rot = <ROT>;
$rot = join("",@rot);
#print "rot:$rot\n";
## get transformation matrix from new xfm file (3x3 matrix + 1 column for translation)
$float = '[-]?\d+\.?\d*[Ee]?[+-]?\d*';
$i=0;
foreach $line (@rot){ # each line in xfm file
 if($line =~ /\s+($float)\s+($float)\s+($float)\s+$float;?\s*\n/){$m[$i]=$1; $m[$i+1]=$2; $m[$i+2]=$3; $i=$i+3;}
}
#print "m:@m\n";

## prepare for vector interpolation, make sure all vectors are in top(arbitrary) quadrants of 3-D space
## first remove other files or matlab crashes
`rm -f flip*`;
system("matlab_cmdline.pl \"flip_vector('$input[0]','$input[1]','$input[2]')\"");

## apply rotation to vectors
`minccalc -clobber -expression $m[0]*A[0]+$m[1]*A[1]+$m[2]*A[2] flip_x.mnc flip_y.mnc flip_z.mnc tx_$x$ext`;
`minccalc -clobber -expression $m[3]*A[0]+$m[4]*A[1]+$m[5]*A[2] flip_x.mnc flip_y.mnc flip_z.mnc tx_$y$ext`;
`minccalc -clobber -expression $m[6]*A[0]+$m[7]*A[1]+$m[8]*A[2] flip_x.mnc flip_y.mnc flip_z.mnc tx_$z$ext`;

##`minccalc -clobber -expression $m[0]*A[0]+$m[1]*A[1]+$m[2]*A[2] $vector[0] $vector[1] $vector[2] tx_$x$ext`;
##`minccalc -clobber -expression $m[3]*A[0]+$m[4]*A[1]+$m[5]*A[2] $vector[0] $vector[1] $vector[2] tx_$y$ext`;
##`minccalc -clobber -expression $m[6]*A[0]+$m[7]*A[1]+$m[8]*A[2] $vector[0] $vector[1] $vector[2] tx_$z$ext`;

##apply regular transformation
## now that we've taken care of the ambiguous directionality, should be able to just interpolate
#`mincresample -clobber new_x.mnc -like $vector[0] -transformation $xfm res_tx_$x$ext`;
#`mincresample -clobber new_y.mnc -like $vector[0] -transformation $xfm res_tx_$y$ext`;
#`mincresample -clobber new_z.mnc -like $vector[0] -transformation $xfm res_tx_$z$ext`;

## now that we've taken care of the ambiguous directionality, should be able to just interpolate
`mincresample -clobber -tricubic -use_input tx_$x$ext -transformation $xfm res_tx_$x$ext`;
`mincresample -clobber -tricubic -use_input tx_$y$ext -transformation $xfm res_tx_$y$ext`;
`mincresample -clobber -tricubic -use_input tx_$z$ext -transformation $xfm res_tx_$z$ext`;
## DO: should clean these up

## have to re-normalize the vectors
`minccalc -clobber -expression A[0]^2+A[1]^2+A[2]^2 res_tx_$x$ext res_tx_$y$ext res_tx_$z$ext norm$ext`;
`mincmath -clobber -sqrt norm$ext norm2$ext`;
`mv norm2$ext norm$ext`;
`mincmath -clobber -div res_tx_$x$ext norm$ext norm_res_tx_$x$ext`;
`mincmath -clobber -div res_tx_$y$ext norm$ext norm_res_tx_$y$ext`;
`mincmath -clobber -div res_tx_$z$ext norm$ext norm_res_tx_$z$ext`;
## clamp volume
`mincmath -clobber -clamp -const2 -1 1 norm_res_tx_$x$ext norm2_res_tx_$x$ext`;
`mincmath -clobber -clamp -const2 -1 1 norm_res_tx_$y$ext norm2_res_tx_$y$ext`;
`mincmath -clobber -clamp -const2 -1 1 norm_res_tx_$z$ext norm2_res_tx_$z$ext`;
`mv norm2_res_tx_$x$ext norm_res_tx_$x$ext`;
`mv norm2_res_tx_$y$ext norm_res_tx_$y$ext`;
`mv norm2_res_tx_$z$ext norm_res_tx_$z$ext`;
## NOTE: the files should END with _e1{x,y,z}.mnc to be recognized by the -namebase option of Fibretrack
## clean up
print "Removing temporary files...\n";
`rm -f res_tx_*.mnc tx_*.mnc`;
`rm -f rot.xfm`;
`rm -f norm.mnc`;
`rm -f flip*`;

## write nicer output file name
($x2,$dir,$ext) = fileparse($input[0],'\..*');
($y2,$dir,$ext) = fileparse($input[1],'\..*');
($z2,$dir,$ext) = fileparse($input[2],'\..*');

`mv norm_res_tx_$x$ext $x2-res$ext`;
`mv norm_res_tx_$y$ext $y2-res$ext`;
`mv norm_res_tx_$z$ext $z2-res$ext`;

exit(1);


#### Using matlab instead of minccalc####
## apply vector transformation on resampled volumes
$matrix = "\[$m[0] $m[1] $m[2];$m[3] $m[4] $m[5];$m[6] $m[7] $m[8] \]";
system("matlab_cmdline.pl \"transform_vector('resamp_x$ext','resamp_y$ext','resamp_z$ext', $matrix)\"");



#############################################
#rmext: remove last extension if there is one
#############################################


sub rmext {
my($in) = @_;
if (index($in,".")<0) {
return $in;
}
else {
my($flag) = 0;
while (chop($in) ne ".") {}
return $in
}
}

