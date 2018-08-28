#! /usr/bin/env perl

###################
### This function will write diffusion directions from bvecs bvals to a mincfile header 
###
###  When using dcm2nii, .bvec and .bval files are ouputted, giving the correct number of elements for all the b=0
###      diffusion images, we can the write this directly to the minc header
###
### HOWEVER: it appears that the sign of the y and z directions are flipped !!
###
### Writtent by Ilana R. Leppert
###################

require 5.001;
use Getopt::Tabular;
use MNI::Startup qw(nocputimes);
use MNI::Spawn;
use MNI::FileUtilities qw(check_output_dirs);
use File::Basename;

if($0 =~ /[\/A-Za-z_0-9]+\/([A-Za-z0-9_]+)\.pl/) {$program = $1;}	#program name without path
$Usage = <<USAGE;

Usage: $program.pl -bvec <file.bvec> -bval <file.bval> <file.mnc>
-help for options


Write bvec and bval files (as used by .nii diffusion files) to a diffusion minc file.

USAGE

@args_table = (["-bvec","string",1,\$bvec,"b vector file"],
	       ["-bval","string",1,\$bval,"b values files"]
);

Getopt::Tabular::SetHelp ($Usage, '');

GetOptions(\@args_table, \@ARGV,\@args)|| exit 1;
$input=pop(@args);

## open vector and bvalue files
open(BVEC,"< $bvec") or die "can't open $bvec: $!";
@bvecs = <BVEC>; #slurp all files, each line at an index
open(BVAL,"< $bval") or die "can't open $bval: $!";
@bval = <BVAL>; #slurp all files, each line at an index

@b_values = split(" ", $bval[0]);
$b_values = join(",", @b_values);
@x_direction = split(" ", $bvecs[0]);
$x_directions = join(",", @x_direction);
### Y and Z polarities have to be changed!!!!
@y_direction = split(" ", $bvecs[1]);
@z_direction = split(" ", $bvecs[2]);
for ($i=0; $i<scalar(@y_direction); $i++){
	push(@y,-1*$y_direction[$i]);
	push(@z,-1*$z_direction[$i]);
}
$y_directions = join(",", @y);
$z_directions = join(",", @z);

print "bvalues: $b_values\n\n xdir: $x_directions\n\n ydir: $y_directions\n\n zdir: $z_directions\n\n";


### add to mincheader
`minc_modify_header -dinsert acquisition:bvalues=$b_values $input`;
`minc_modify_header -dinsert acquisition:direction_x=$x_directions $input`;
`minc_modify_header -dinsert acquisition:direction_y=$y_directions $input`;
`minc_modify_header -dinsert acquisition:direction_z=$z_directions $input`;
