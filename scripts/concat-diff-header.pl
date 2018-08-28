#! /usr/bin/env perl

###################
### This function will concatenate minc files, as well as the diffusion directions in the header
###
###	usage: mincconcat-diff.pl [<options>] {<infile1> ...] <outfile> 
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

Concatenate minc files along the time dimension as well as diffusion information in header.


Usage: $program [<infile1> ...] <outfile>
       $program -help 	 Print summary of command-line options and abort

USAGE

@args_table = (["-clobber","boolean",0,\$clobber,"clobber existing file"]
);

Getopt::Tabular::SetHelp ($Usage, '');

GetOptions(\@args_table, \@ARGV, \@args) || exit 1;



### Get input and output file names from ARGV
$output_index=scalar @args;
$output = $args[$output_index-1];
@inputs = @args[0..$output_index-2];


### Create new array with combined diffusion info
$bval=$xdir=$ydir=$zdir='';


foreach $file (@inputs){
    
	$tmp = `mincinfo -attvalue acquisition:bvalues $file`;
	## replace last carriage return with space
	$bval=$bval.$tmp;
	$bval=~s/\s{1}/,/g;
	##have to remove last comma
	chop($bval);
	
	$tmp = `mincinfo -attvalue acquisition:direction_x $file`;
	$xdir=$xdir.$tmp;
	$xdir=~s/\s{1}/,/g;
	#print "xdir:$xdir\n";
	##have to remove last comma
	chop($xdir);
	
	$tmp = `mincinfo -attvalue acquisition:direction_y $file`;
	$ydir=$ydir.$tmp;
	$ydir=~s/\s{1}/,/g;
	##have to remove last comma
	chop($ydir);
	
	$tmp = `mincinfo -attvalue acquisition:direction_z $file`;
	$zdir=$zdir.$tmp;
	$zdir=~s/\s{1}/,/g;
	##have to remove last comma
	chop($zdir);
}
## have to remove last comma
chop($bval);
chop($xdir);
chop($ydir);
chop($zdir);

#print "bval:$bval\n";
##Concatenate volume
if($clobber==0){$option='';}else{$option='-clobber';}
`mincconcat $option -sequential -concat_dimension time @inputs $output`;


##Write concatenated info to header
`minc_modify_header -dinsert acquisition:bvalues=$bval $output`;
`minc_modify_header -dinsert acquisition:direction_x=$xdir $output`;
`minc_modify_header -dinsert acquisition:direction_y=$ydir $output`;
`minc_modify_header -dinsert acquisition:direction_z=$zdir $output`;
