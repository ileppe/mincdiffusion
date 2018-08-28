#! /usr/bin/perl -w

###################
### This function will create a gradient direction file compatible with FSL from a test file
###
###	usage: mghgrad2FSL.pl -dirs <direction text file> -b <bvalue> -nb0 <num_b0> <bvals> <bvecs>
###
###  grad dirs file: list of diffusion gradient directions (1 direction per line, separated by " " or ",")
###  b : bvalue
###  num_b0 : number of b=0 images
###  output bvals : text file of bvals used by FSL
###  output bvecs: text file with grad directions used by FSL
###			check out:	http://www.fmrib.ox.ac.uk/fsl/fdt/index.html
###################


require 5.001;
use Getopt::Tabular;
use MNI::Startup qw(nocputimes);
use MNI::Spawn;
use MNI::FileUtilities qw(check_output_dirs);
use File::Basename;
$num_b0=1;

if($0 =~ /[\/A-Za-z_0-9]+\/([A-Za-z0-9_]+)\.pl/) {$program = $1;}	#program name without path
$Usage = <<USAGE;

Usage: $program.pl  -dirs <direction text file> -b <bvalue> -nb0 <num_b0> <bvals> <bvecs>
-help for options

USAGE

@args_table = (["-dirs","string",1,\$dirs,"file of directions(format: 1 direction/line, separated by \" \" or \",\")"],
 	       ["-b","string",1,\$bvalue,"b value"],
	       ["-nb0","string",1,\$num_b0,"number of b=0 images"]
);

Getopt::Tabular::SetHelp ($Usage, '');

GetOptions(\@args_table, \@ARGV,\@args)|| exit 1;
$bvecs=pop(@args);
$bvals=pop(@args);

## open vector and bvalue files
open(BVEC,"> $bvecs") or die "can't open $bvecs: $!";
open(BVAL,"> $bvals") or die "can't open $bvals: $!";
#print "n_dir:$n_dir\n\n";


## get num of directions
open(FILE,"< $dirs") or die "can't open $dirs: $!";
@dirs = <FILE>; #slurp all files, each line at an index
#print "n_dir:$n_dir\n\n";

## fill b_values array (should have $n_scans elements with $n_b_0 of them = to 0 and the rest = acquisition:b_value
## fill direction arrays (should have $n_scans elements with $n_b_0 of them = to 0 and the rest from input direction file

for ($i=0; $i < $num_b0; $i++){
  push(@b, 0);
  push(@x, 0);
  push(@y, 0);
  push(@z, 0);
}

## make an array with info for each direction
foreach $dir (@dirs){
  if($dir =~ /\s*([0-9.\-]+),*\s*([0-9.\-]+),*\s*([0-9.\-]+)/){push(@x,$1);push(@y,$2);push(@z,-1*$3); } else {die "can't get directions from $dirs at line $dir\n"}  
  #print "x:$x\t y:$y\t z:$z\n";
  push(@b,$bvalue); ## add values to b_values array
}
#print "x: @x size : $#x\n\n";
#print "b:@b\ size : $#b\n\n";


## Make text files as expected by FSL
################NOTE: BVECS should be normalized to 1 !!!


## bvals text file
# | an example bvals would look as such:
# |   # bvalues  
## b1 b2 b3 b4 ...bn

## bvecs text file
# | an example bvecs would look as such:
# |   gradient directions
## x1 x2 x3 ...xn
## y1 y2 y3 ...yn
## z1 z2 z3 ...zn

print BVAL "@b\n";
print BVEC "@x\n";
print BVEC "@y\n";
print BVEC "@z";

