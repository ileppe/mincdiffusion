#! /usr/bin/env perl


###################
### This function will create a gradient direction file compatible with FSL from directions in the mincheader
###
###	usage: header2mincdti.pl <dti-series.mnc> <bvec> <bval>
###
###  dti-series.mnc: input file (to get direction and bvalues)
###  output_file: output file with b and grad directions suitable for mncdti
###################

require 5.001;
use Getopt::Tabular;
use MNI::Startup qw(nocputimes);
use MNI::Spawn;
use MNI::FileUtilities qw(check_output_dirs);
use File::Basename;

if($0 =~ /[\/A-Za-z_0-9]+\/([A-Za-z0-9_]+)\.pl/) {$program = $1;}	#program name without path
$Usage = <<USAGE;

Usage: $program <dti-series.mnc> <bvecs> <bvals>

This function will create a gradient direction file (bvecs) and a bvals file compatible with FSL from directions in the mincheader

USAGE
@args_table = ();

Getopt::Tabular::SetHelp ($Usage, '');

GetOptions(\@args_table, \@ARGV, \@args) || exit 1;
(@args > 0) || die $Usage;

$inputfilename=$args[0];
if(!(-e $inputfilename)){die "ERROR Input file $inputfilename does not exist\n";}
$fbvec=$args[1];
$fbval=$args[2];
open(FBVEC,"> $fbvec") or die "can't open $fbvec: $!";
open(FBVAL,"> $fbval") or die "can't open $fbval: $!";

## get directions
$tmp = `mincinfo -attvalue acquisition:direction_x $inputfilename`;
chomp($tmp);
@x = split(/ /,$tmp);
$tmp = `mincinfo -attvalue acquisition:direction_y $inputfilename`;
chomp($tmp);
@y = split(/ /,$tmp);
$tmp = `mincinfo -attvalue acquisition:direction_z $inputfilename`;
chomp($tmp);
@z = split(/ /,$tmp);

## get bvalues
$tmp = `mincinfo -attvalue acquisition:bvalues $inputfilename`;
chomp($tmp);
@b = split(/ /,$tmp);

@yneg='';
$tmp = scalar(@y);
## have to change the sign of the y direction!! NO! Seems like the newer FSL does no require this
for($i=0;$i<$tmp;$i++){
  $yneg[$i]=1*$y[$i];
 }

@zneg='';
$tmp = scalar(@z);
## have to change the sign of the y direction!!NO! Seems like the newer FSL does no require this
for($i=0;$i<$tmp;$i++){
  $zneg[$i]=1*$z[$i];
 }

## bvec file should look like
#x1 x2 x3 ... xn
#y1 y2 y3 ... yn
#z1 z2 z3 ... zn
##
#####################################
## FSLVIEW is not displaying the directions correctly, correct dyads does not mean correct tracking!!!!!
##############################
#print FBVEC "@x\n@yneg\n@zneg";
print FBVEC "@x\n@yneg\n@zneg"; 
## bval file should look like
#b1 b2 b3 b4 ... bn
## CARRIAGE RETURN!!
##
print FBVAL "@b\n";
