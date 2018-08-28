#! /usr/bin/env perl
#
#perl script to run matlab script dti or dti_octave.m
#
#
#only works for transverse images right now (so as to simplify this code)
#
#DO: remove dti call entirely soon and just use niak always through dti_octave


require 5.001;
use Getopt::Tabular;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use strict;

my $keep_tmp=0;
my $maskname='';
my $sigma=0;
my $outputdir='.';
my $me=basename($0);
my $mydir=dirname($0);
my $octave=0;#use octave
my $clobber=0;
my $fake=0;
my $verbose=0;
my $niak=1; #always use by default
my $niak_dir='';
my $split=0;

my $Usage = <<USAGE;

Script to run dti.m or dti_octave.m: calculates mean diffusivity, FA, eigenvectors and eigenvalues of diffusion tensor

Usage: $me [options] <input.mnc>

input:  minc file with b values in time dimension
output: by default, output filenames are based on the input filename and written to the current directory.
        Use '-outputdir' option for arbitrary output directory. 


-help for options

will resample everything to transverse images!!!

to use this script, you will need: mincinfo, mincreshape, mincconcat, matlab or octave,  

USAGE

my @args_table = (["-clobber","boolean",undef,\$clobber,"clobber all output files"],
                  ["-sigma","string",1,\$sigma,"standard deviation of the noise in the image (sigma=mean intensity in noise * sqrt(2/pi)).  If this option is not given or is 0, no chisquare map will be produced.",],
                  ["-mask","string",1,\$maskname,"mask file: tensor will not be calculated where mask file value is 0"],
                  ["-outputdir","string",1,\$outputdir,"output directory, will use the input filename base."],
                  ["-octave","boolean",undef,\$octave,"Use octave"],
                  ["-verbose","boolean",undef,\$verbose,"Be verbose"],
                  #["-niak","boolean",undef,\$niak,"Run NIAK version (always used for octave; not necessary if specify -niakdir input)"],
                  ["-niakdir","string",1,\$niak_dir,"Specify NIAK directory (if not already on your path)"],
                  ["-split","boolean",undef,\$split,"Split input into subsets to avoid memory problems in matlab"]
);

Getopt::Tabular::SetHelp ($Usage, '');

my @args;
GetOptions(\@args_table, \@ARGV, \@args) || exit 1;
die $Usage unless $#ARGV >=0;

my $inputfilename=$args[0];

if( $maskname && !-e $maskname ){die "ERROR Mask $maskname does not exist\n";}
if(!(-e $inputfilename)){die "ERROR Input file $inputfilename does not exist\n";}

my $tmpdirroot=&tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => !$keep_tmp );

my $compress=$ENV{MINC_COMPRESS};
delete $ENV{MINC_COMPRESS} if $compress;

my $namebase;
my $masknamenodir;
my $maskbase;
my $info;
my $count;

$namebase=rmext(rmext(basename($inputfilename))); #called twice in case zipped

if ($outputdir){
  if(!(-e $outputdir)){do_cmd("mkdir -p $outputdir");} # only create if doesn't exist
}else{
  $outputdir='./';
}


$namebase="$outputdir/$namebase";

if($maskname) 
{
    $masknamenodir=rmext(basename($maskname));
    $maskbase="$outputdir/$masknamenodir";
}


unless($clobber) {
  check_file("$namebase\_MD.mnc");
  check_file("$namebase\_e1x.mnc");
  check_file("$namebase\_e1y.mnc");
  check_file("$namebase\_e1z.mnc");
  check_file("$namebase\_e2x.mnc");
  check_file("$namebase\_e2y.mnc");
  check_file("$namebase\_e2z.mnc");
  check_file("$namebase\_lambda1.mnc");
  check_file("$namebase\_lambda2.mnc");
  check_file("$namebase\_lambda3.mnc");
  check_file("$namebase\_FA.mnc");
  check_file("$namebase\_e3x.mnc");
  check_file("$namebase\_e3y.mnc");
  check_file("$namebase\_e3z.mnc");
}

if ($maskname) 
{
  do_cmd("mincresample -nearest -like $inputfilename $maskname $maskbase-diffspace.mnc");
  duplicate_frames("$maskbase-diffspace.mnc",$inputfilename,"$tmpdirroot/diffspace-nframes.mnc");
  do_cmd("mincmath -float -nocheck_dimensions -copy_header -mult $inputfilename $tmpdirroot/diffspace-nframes.mnc $tmpdirroot/mask.mnc -clob");
  ## resample to transverse if is not already
  $info=`mincinfo $inputfilename`;
  if($info !~/time zspace yspace xspace/){
     print "Resampling image to transverse...\n";
     do_cmd("mincresample -transverse $tmpdirroot/mask.mnc $tmpdirroot/mask-trans.mnc");
     $inputfilename="$tmpdirroot/mask-trans.mnc";
  } else { 
    $inputfilename="$tmpdirroot/mask.mnc";
  }
  
} else {
 ## resample to transverse
 $info=`mincinfo $inputfilename`;
 if($info !~/time zspace yspace xspace/){
   print "Resampling image to transverse...\n";
   do_cmd("mincresample -transverse $inputfilename $tmpdirroot/trans.mnc");
   $inputfilename="$tmpdirroot/trans.mnc"; 
 }
}
 
##########################################################################################
# input subsets of data at a time: 
##########################################################################################

#get number of slices in input file:

my ($xdim,$ydim,$zdim,$time);

chomp($zdim=`mincinfo -dimlength zspace $inputfilename`); #this works only for transverse images; zspace is slice
chomp($ydim=`mincinfo -dimlength yspace $inputfilename`);
chomp($xdim=`mincinfo -dimlength xspace $inputfilename`);
chomp($time=`mincinfo -dimlength time   $inputfilename`);

my $steps=0;

if ($zdim/20>1 && $split) {
  my $i;
  for ($i=0; $i<=$zdim-20; $i+=20) {
    $steps+=1;
    do_cmd("mincreshape -dimrange time=0,$time -dimrange zspace=$i,20 -dimrange yspace=0,$ydim  -dimrange xspace=0,$xdim $inputfilename $tmpdirroot/file$steps.mnc");
  }


  if ($zdim-$i != 0) {

    $steps+=1;

    $count=$zdim-$i;

    do_cmd("mincreshape -dimrange time=0,$time -dimrange zspace=$i,$count -dimrange yspace=0,$ydim  -dimrange xspace=0,$xdim $inputfilename $tmpdirroot/file$steps.mnc");

  }

} else {
  $steps=1;
  do_cmd("cp $inputfilename $tmpdirroot/file$steps.mnc");
}

##########################################################################################
# input these to matlab dti.m or dti_octave.m (runs niak):
##########################################################################################

my $i;
for ($i=1; $i<=$steps; $i++) {

my $command='dti_octave';

if ($sigma>0)
  {
    if($octave) 
    {
      run_octave("$command('$tmpdirroot/file$i.mnc',$sigma)");
    } else {
      run_matlab("$command('$tmpdirroot/file$i.mnc',$sigma)");
    }
  }
else
  {
    if($octave) 
    {
      run_octave("$command('$tmpdirroot/file$i.mnc')");
    } else {
      run_matlab("$command('$tmpdirroot/file$i.mnc')");
    }
  }
}


$ENV{MINC_COMPRESS}=$compress if $compress;

my @files=qw/ MD e1x e1y e1z e2x e2y e2z e3x e3y e3z lambda1 lambda2 lambda3 FA /;
push @files,"chi2" if $sigma>0;


foreach $i(@files) {
  if($steps>1) {
      
    do_cmd("mincconcat -clobber -concat_dimension zspace $tmpdirroot/*${i}* ${namebase}_${i}.mnc");
  } else {
    do_cmd("cp $tmpdirroot/*${i}* ${namebase}_${i}.mnc");
  }
}

rgb_vector("${namebase}_e1x.mnc","${namebase}_e1y.mnc","${namebase}_e1z.mnc","${namebase}_FA.mnc","${namebase}_rgb.mnc");


#################################################################
#functions
#################################################################



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


sub do_cmd {
    print STDOUT "@_\n" if $verbose;
    if(!$fake) {
        system(@_) == 0 or die "DIED: @_\n";
    }
}

sub check_file {
  die("${_[0]} exists!\n") if -e $_[0];
}


sub run_matlab {
  my ($command)=@_;
  #system("(echo \"$command\"| matlab -nosplash $logging)");
  open MATLAB,"|matlab -nosplash -nojvm -nodisplay " or die "Can't start matlab!\n";
  print MATLAB "addpath('$mydir')\n";
  print MATLAB "addpath(genpath('$niak_dir'))\n" if $niak_dir;
  print MATLAB $command,"\n";
  close MATLAB;
  die "Error in MATLAB!\n" if $?;
}

sub run_octave {
  my ($command)=@_;
  #system("(echo \"$command\"| matlab -nosplash $logging)");
  open OCTAVE, "|octave  " or die "Can't start octave!\n";
  print OCTAVE "addpath('$mydir')\n";
  print OCTAVE "addpath(genpath('$niak_dir'))\n" if $niak_dir;
  print OCTAVE $command,"\n";
  close OCTAVE;
  die "Error in MATLAB!\n" if $?;
}

sub duplicate_frames {
  my ($input,$like,$output)=@_;
  my $samples=`mincinfo -dimlength time $like`;
  chomp($samples);
  my $i;
  my @dummy;
  for($i=0;$i<$samples;$i+=1)
  {
    push @dummy,$input;
  }
  do_cmd('mincconcat','-concat_dimension','time',@dummy,$output,'-clobber');
}

sub rgb_vector {
  my ($x,$y,$z,$fa,$out)=@_;
  do_cmd('minccalc','-express','A[1]*abs(A[0])',$x,$fa,"$tmpdirroot/x.mnc");
  do_cmd('minccalc','-express','A[1]*abs(A[0])',$y,$fa,"$tmpdirroot/y.mnc");
  do_cmd('minccalc','-express','A[1]*abs(A[0])',$z,$fa,"$tmpdirroot/z.mnc");
  do_cmd('mincmakevector',"$tmpdirroot/x.mnc","$tmpdirroot/y.mnc","$tmpdirroot/z.mnc",$out,'-clobber');
}
