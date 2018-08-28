#! /usr/bin/perl -w

###################
### This function will modify the minc header to include info about b-values and gradient directions
###
###	usage: minc_modify_header-diff.pl -diffdirs <direction_file> -bval <b_value> -num_b0 <num> -coord_system <xxx>
###
###  direction_file: list of diffusion gradient directions (formatted like: /data/mril/mril5/ilana/DTI_docs/MGH/gradient_mgh_dti)
###				{x,y,z}
###		note: number of lines in this file should be number of directions!!
###  input_file: input file for which to change header
###################
### dicom y = -minc y 
### DO: check Phase-Read-Slice (PRS) when from bic_ep2d_diff or jc_diffNdirs_400_4 - should be A-P acquisition


require 5.001;
use Getopt::Tabular;
use MNI::Startup qw(nocputimes);
use MNI::Spawn;
use MNI::FileUtilities qw(check_output_dirs);
use File::Basename;
#use ileppe_lib;

if($0 =~ /[\/A-Za-z_0-9]+\/([A-Za-z0-9_-]+\.pl)/) {$program = $1;}#program name without path
elsif($0 =~ /([A-Za-z0-9_-]+\.pl)/) {$program = $1;}
	
$Usage = <<USAGE;

Usage: $program [options] <input.mnc>
-help for options

USAGE

@args_table = (["-diffdirs","string",1,\$dir_file,"list of diffusion gradient directions (formatted like ~ilana/DTI_docs/MGH/gradient_mgh_dti/gradient_mgh_dti*.gdt))"],
	       ["-bval","string",1,\$b_value,"b value"],
	       ["-num_b0","string",1,\$num_b0,"number of baseline images"],
	       ["-b_baseline","integer",1,\$b_baseline,"b value of baseline images"],
	       ["-coord_system","string",1,\$coord_system,"coordinate system of vectors (xyz) or (prs)"]
);

Getopt::Tabular::SetHelp ($Usage, '');
$coord_system='xyz';
$b_baseline= 0; # default b value for baseline images is 0

GetOptions(\@args_table, \@ARGV,\@args) || exit 1;

$input=$args[0];

##get coordinate system
if($coord_system eq 'xyz'){$type=1;}
elsif($coord_system eq 'prs'){$type=2;}
else {die "Incorrect coordinate system $coord_system\n";}

## get num of directions
open(FILE,"< $dir_file") or die "can't open $dir_file: $!";
@dirs = <FILE>; #slurp all files, each line at an index
#$n_dir = scalar(@dirs); # num directions
#print "n_dir:$n_dir\n\n";

### DO: need to ask bval and number of b0 images from user because many protocols do not indicate in header
## get num of b=0 and b n.e. 0 scans
#$n_scans =  sprintf("%.0f", `mincinfo -attvalue acquisition:num_dyn_scans $input`); # total num of scans
#$n_b_0 = $n_scans - $n_dir; #number of scans with b=0
##print "n_b_0:$n_b_0\n";

#$b_value = sprintf("%.0f", `mincinfo -attvalue acquisition:b_value $input`);
#print "b:@b\n\n";

$n_b_0 = $num_b0;

## fill b_values array (should have $n_scans elements with $n_b_0 of them = to 0 and the rest = acquisition:b_value
## fill direction arrays (should have $n_scans elements with $n_b_0 of them = to 0 and the rest fomr input direction file

for ($i=0; $i < $n_b_0; $i++){
  push(@b, $b_baseline);
  push(@x, 0);
  push(@y, 0);
  push(@z, 0);
}



$d=0; # count number of directions
## make an array with info for each direction
if($type == 1){ ## xyz coordinate system
 ### NOTE DICOM y = -MINC Y
 foreach $dir (@dirs){
   if($dir =~ /\s*([0-9.\-]+),*\s*([0-9.\-]+),*\s*([0-9.\-]+)/){push(@x,$1);push(@y,-1*$2);push(@z,$3); } else {die "can't get directions from $dir_file at line $dir\n"}  
   #print "x:$x\t y:$y\t z:$z\n";
   $d++;
   push(@b,$b_value); ## add values to b_values array
 }
}
elsif($type == 2){ ## prs coordinate system
 ### NOTE if A-P acquisition: P=-Y  R=X S=Z
 foreach $dir (@dirs){
   if($dir =~ /\s*([0-9.\-]+),*\s*([0-9.\-]+),*\s*([0-9.\-]+)/){push(@y,-1*$1);push(@x,$2);push(@z,$3); } else {die "can't get directions from $dir_file at line $dir\n"}  
   #print "x:$x\t y:$y\t z:$z\n";
   $d++;
   push(@b,$b_value); ## add values to b_values array
 }
}
print "There are $d directions...\n";

#print "x: @x size : $#x\n\n";
#print "b:@b\ size : $#b\n\n";
## make one string with elements separated by commas
$b_values = join(",", @b);
$x_direction = join(",", @x);
$y_direction = join(",", @y);
$z_direction = join(",", @z);

#print "z_direction:$z_direction\n\n";

### add to mincheader
`minc_modify_header -dinsert acquisition:bvalues=$b_values $input`;
`minc_modify_header -dinsert acquisition:direction_x=$x_direction $input`;
`minc_modify_header -dinsert acquisition:direction_y=$y_direction $input`;
`minc_modify_header -dinsert acquisition:direction_z=$z_direction $input`;
