#! /usr/bin/env perl

##################################################################
## This script can be used to preprocess diffusion data, steps can
## be added or removed as seem fit by user
##
##Indludes:
##	- writting diffusion info to the header
##	- non-uniformity correction of the anatomical
##	- averaging of diffusion series
##	- registration to the anatomical
##	- creation of a mask
##################################################################
##
##  Usage: diff_preprocess.pl -dirs <direction_file> -b <b_value> -b0 <num> -cod <xxx> -anat <t1_anatomical.mnc> <input1.mnc ...> 
##
##  -dirs : list of directions formatted like ~ilana/DTI_docs/MGH/MGHdiffgrad/diffgrad*.txt or gradient_mgh_dti/gradient_mgh_dti*.gdt
##  -b : b value
##  -b0 : number of b0 images
##  -cod : coordinate system (xyz) or (prs)
##  -anat : anatomical image
##  -input : series of DTI images
##################################################################

require 5.001;
use Getopt::Tabular;
use MNI::Startup qw(nocputimes);
use MNI::Spawn;
use MNI::FileUtilities qw(check_output_dirs);
use File::Basename;

if($0 =~ /[\/A-Za-z_0-9]+\/([A-Za-z0-9_]+\.pl)/) {$program = $1;}	#program name without path
elsif($0 =~ /([A-Za-z0-9_-]+\.pl)/) {$program = $1;}

$Usage = <<USAGE;

Usage: $program [options] [-dirs <dirs.txt> -b <bvalue> -b0 <num>] [-b_baseline <bvalue>] [-cod <xyz or prs>] 
-anat <anatomical.mnc> <input1.mnc ...> <output.mnc> 
Note: by default, multiple series will be concatenated, not averaged, the SNR gain is the same
Intermediate output files are stored in the same directory as the output file.
-help for options

USAGE

$coord_system='xyz';
$clobber=0;
$b_baseline = 0;# default baseline b value is 0
$dir_file = 'already in mincheader'; ## if no input file given, assume info is already in header
$b_value = 'already in mincheader';
$num_b0 = 'already in mincheader';
$keep=0;

@args_table = (["-clobber","boolean",undef,\$clobber,"clobber all output files"],
	       ["-dirs","string",1,\$dir_file,"list of diffusion gradient directions (formatted like ~ilana/links/DTI_docs/MGH/gradient_mgh_dti/gradient_mgh_dti*.gdt)"],
	       ["-b","string",1,\$b_value,"b value"],
	       ["-b0","string",1,\$num_b0,"number of baseline images"],
	       ["-b_baseline","integer",1,\$b_baseline,"b value of baseline images"],
	       ["-cod","string",1,\$coord_system,"coordinate system of vectors (xyz) or (prs)"],
	       ["-anat","string",1,\$anatfullname,"anatomical image"],
	       ["-keep","boolean",undef,\$keep,"keep all intermediate files (stored in tmp dir)"]
);

Getopt::Tabular::SetHelp ($Usage, '');

GetOptions(\@args_table, \@ARGV,\@args) || exit 1;

## parse inputs and ouput
$output=pop(@args);

## clobber option
if($clobber == 1){$clob = '-clobber';}else{$clob='';}

if(-e($output) & ($clobber ==0) ){die "Output file $output already exixts and noclobber\n"};

## output directory
($tmp, $outdir, $ext) = fileparse($output,'\..*');

## temp directory
chomp($unique=`date +%m%w%H%M%S`);  
$tmpdirroot="tmp$unique";  
system("mkdir /tmp/$tmpdirroot");

foreach $dtifullname (@args){

   ## check file permissions
    $ls=`\\ls -l $dtifullname`;
    if($ls =~ /^([a-z-]{9})/){$perm=$1;} else {die "Can't get file permissions of $dtifullname\n";}
    if($perm =~ /w/){}else{print "Permission denied for $dtifullname......Change permissions with chmod or make a copy\n"; die}
}

($anatname, $dir, $ext) = fileparse($anatfullname,'\..*');

## Apply non-uniformity correction twice to the anatomical
if(!(-e "$anatname-n3.mnc")){ # only if not done
  `nu_correct $anatfullname $outdir/$anatname-n3.mnc`;
  `nu_correct $outdir/$anatname-n3.mnc $outdir/$anatname-n3-2.mnc`;
  `mv -f $outdir/$anatname-n3-2.mnc $outdir/$anatname-n3.mnc`;
}
## flip images to get pos step value for spatial axes
`mincreshape +direction $outdir/$anatname-n3.mnc $outdir/$anatname-n3+dirs.mnc`;
`mv -f $outdir/$anatname-n3+dirs.mnc $outdir/$anatname-n3.mnc`;

### Some processing for each input dti series
#@dtinames=@dirs=@exts=();
foreach $dtifullname (@args){

  ##Get path info for each input file
  ($dtiname, $dir, $ext) = fileparse($dtifullname,'\..*');
  #push(@dtinames,$dtiname);
  #push(@dirs,$dir);
  #push(@exts,$ext);
  
  ### write directions and bvalues to header  of dti-series if not already there
  if ($dir_file ne 'already in mincheader'){
    print "Writing directions to header of $dtifullname\n";
    system("minc_modify_header-diff.pl -diffdirs $dir_file -bval $b_value -num_b0 $num_b0 -b_baseline $b_baseline -coord_system $coord_system $dtifullname")==0 or die;
  }else{ print "Using directions already in header of $dtifullname\n";}  

  ## get first frame of DTI series (b=0 image)
  print "Extracting first frame from $dtifullname...\n\n";
  `mincreshape -clobber $dtifullname +direction -dimrange time=0 $outdir/$dtiname-frame0.mnc`;
  # keep the b=0 image from the first series as the target
  if($dtifullname eq $args[0]){
  	$target="$outdir/$dtiname-frame0.mnc";
	$series[0]=$dtifullname;
  }
  ##register all series to the first
  else{ # don't do this for the first series
  print "Registration within session...\n";
  `mritoself $clob -close $outdir/$dtiname-frame0.mnc $target $outdir/$dtiname-to-self.xfm`;
  `mincresample $clob -use_input_sampling -transformation $outdir/$dtiname-to-self.xfm $dtifullname $outdir/$dtiname-reg-with-self.mnc`;
  # tranform diffusion directions
  `transform_diff_dirs -xfm $outdir/$dtiname-to-self.xfm $outdir/$dtiname-reg-with-self.mnc`;
  
  ### Add registered series 
  push(@series,"$outdir/$dtiname-reg-with-self.mnc");
  
  }
}
if(scalar(@args) > 1){ # more than 1 input
  ##if($avg == 1){
  	##print "Averaging multiple runs...\n";
	##`mincaverage -copy_header @series /tmp/dti-avg.mnc`;
  ##}else{
  	##concat all series together (same SNR gain as averaging)
  	print "Concatenating all runs...\n";
	print "concat-diff-header.pl @series /tmp/$tmpdirroot/dti-avg.mnc\n";
  	`concat-diff-header.pl @series /tmp/$tmpdirroot/dti-avg.mnc`;
  ##}
  $final_dti="/tmp/$tmpdirroot/dti-avg.mnc";
}else {$final_dti=$args[0];}

## registration  !!!!!IL  removed center input
$reg_out = "$outdir/dti-to-t1.xfm";
if(!(-e $reg_out)){
  print "Registration to the anatomical $anatfullname...\n";
  `minctracc $clob -identity -mi $target $outdir/$anatname-n3.mnc -lsq6 -debug -threshold 30 5 $outdir/dti-to-t1.xfm -simplex 1 -clobber`;
  `minctracc $clob -mi $target $outdir/$anatname-n3.mnc -lsq6 -debug -threshold 30 5 -transformation $outdir/dti-to-t1.xfm $outdir/dti-to-t1b.xfm -simplex 1 -step 2 2 2`;
  `mv -f $outdir/dti-to-t1b.xfm $outdir/dti-to-t1.xfm`;
}

`mincresample $clob -use_input_sampling -transformation $outdir/dti-to-t1.xfm $final_dti $output`;
#unless -e "$dtiname-reg-with-t1.mnc";

##optionally (is transformation is smalll; should put a threshold) transform the diffusion encoding directions and change minc header:
`transform_diff_dirs -xfm $outdir/dti-to-t1.xfm $output`;


## create brain mask
print "Making brain mask...\n"; ## using executable, migt need to recompile on other systems
#`/data/aces/aces1/quarantines/Linux-i686/Dec-20-2006/bin/mincbet $anatname-n3.mnc $anatname-n3-bet -m -n`;
`mincbet $outdir/$anatname-n3.mnc $outdir/$anatname-n3-bet -h 1.15 -m -n`;
#`mv $anatname-n3-bet_mask $anatname-n3-bet_mask.mnc`;



## Clean up some of the processing files (you can chose to keep them all)
if (!$keep){
  `rm -rf /tmp/$tmpdirroot`;
}
#`rm -f *-frame0.mnc`;
#`rm -f *-reg-with-self.mnc`;

