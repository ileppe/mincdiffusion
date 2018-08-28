#! /usr/bin/env perl

###################
### This function reshape minc files (using mincreshape) while maintaining diffusion info in header
###
###	usage: extract_diff_frames.pl [<options>] <infile> <outfile> 
###
###  can have any number of input files but each has to have diffusion info already specified in header
###################

require 5.001;
use Getopt::Tabular;
use MNI::Startup qw(nocputimes);
use MNI::Spawn;
use MNI::FileUtilities qw(check_output_dirs);
use File::Basename;
#use ileppe_lib;

if($0 =~ /[\/A-Za-z_0-9]+\/([A-Za-z0-9_-]+\.pl)/) {$program = $1;}#program name without path
elsif($0 =~ /([A-Za-z0-9_-]+\.pl)/) {$program = $1;}
	

$clobber=0;

$Usage = <<USAGE;

Cuts a slab in the time dimension out of a minc file while handling diffusion header info.

Usage: $program [<options> ...] <infile> <outfile>
       $program [-help]
USAGE

@args_table = (["-clobber","boolean",0,\$clobber,"clobber existing file"],
	       ["-start","integer",1,\$start,"Specifies corner of hyperslab (C conventions for indices)","<start>"],
	       ["-count","integer",1,\$count,"Specifies number of frames to read","<count>"]
);

Getopt::Tabular::SetHelp ($Usage, '');

GetOptions(\@args_table, \@ARGV, \@args) || exit 1;


### Get input and output file names from non-options arguments
$input = $args[0];
$output = $args[1];
if($start+$count> `mincinfo -dimlength time $input`){print "ERROR: Exceeded maximum time dimension\n"; die}
#print "input:$input\toutput:$output\tstart:$start\tcount:$count\n";

### Extract diffusion info
$bval = `mincinfo -attvalue acquisition:bvalues $input`;
chomp($bval);
@bvals = split(' ',$bval),
@bvals=@bvals[$start..$start+$count-1];
$bval=join(',',@bvals);

$xdir = `mincinfo -attvalue acquisition:direction_x $input`;
chomp($xdir);
@xdirs = split(' ',$xdir),
@xdirs=@xdirs[$start..$start+$count-1];
$xdir=join(',',@xdirs);

$ydir = `mincinfo -attvalue acquisition:direction_y $input`;
chomp($ydir);
@ydirs = split(' ',$ydir),
@ydirs=@ydirs[$start..$start+$count-1];
$ydir=join(',',@ydirs);

$zdir = `mincinfo -attvalue acquisition:direction_z $input`;
chomp($zdir);
@zdirs = split(' ',$zdir),
@zdirs=@zdirs[$start..$start+$count-1];
$zdir=join(',',@zdirs);

##Extract volume
if($clobber==0){$option='';}else{$option='-clobber';}
`mincreshape $option -dimrange time=$start,$count $input $output`;
#print "mincreshape $option -dimrange time=$start,$count $input $output\n";

##Write concatenated info to header
`minc_modify_header -dinsert acquisition:bvalues=$bval $output`;
`minc_modify_header -dinsert acquisition:direction_x=$xdir $output`;
`minc_modify_header -dinsert acquisition:direction_y=$ydir $output`;
`minc_modify_header -dinsert acquisition:direction_z=$zdir $output`;
