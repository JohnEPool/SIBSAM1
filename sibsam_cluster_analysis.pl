#!/usr/bin/perl -w
use strict;
use cmdline_args;

if (!$ARGV[0] || $ARGV[0] eq "help"){
    print <<EOF;

This script analyzes output from sibsam_1locus.pl

It requires the following arguments:

-f [file name stem]
-x [cluster number]

There must be a window file with the above stem followed by "_windows.txt"
There must be one or more simulation files from sibsam_cluster.pl with the file stem followed by "_clusterXsims*" (where X is the cluster number and * is any characters)
Output file will have this file stem appended by "_clusterXresults.txt"

EOF
exit 1;
}

my $FileStem = "";
my $ClusterNumber = "";
my %args;
$args{'-f'} = [ \$FileStem, 1];
$args{'-x'} = [ \$ClusterNumber, 1];

cmdline_args::get_options(\%args, \@ARGV);

my $SimFile = $FileStem . '_cluster' . $ClusterNumber . 'sims';
my $WinFile = $FileStem . '_windows.txt';
my $ClusterFile = $FileStem . '_cluster' . $ClusterNumber . '.txt';
my $OutputFile = $FileStem . '_cluster' . $ClusterNumber . 'results.txt';

my $tolerance = 0.05; #height of sim peaks must be within this range of real peak height to accept

#Get cM and depth information from window files
my $CMsX = 0;
my $CMsA = 0;
my $depth = 0;
my @line = ();
my @DepthsXHigh = ();
my @DepthsAHigh = ();
my @DepthsXLow = ();
my @DepthsALow = ();
my @CMStopsX = ();
my @CMStopsA = ();
my @AutoCMBreaks = ();
my @WinChrs = ();
my @AncDiffs = ();

open W, "<$WinFile" or die;
scalar (<W>);
while (<W>){
  chomp;
  last if m/^$/;
  @line = split;
  if ($line[0] =~ m/X/){
    push @CMStopsX, $line[3];
    $depth = int($line[4]);
    push @DepthsXHigh, $depth;
    $depth = int($line[5]);
    push @DepthsXLow, $depth;
  }
  else{
    push @CMStopsA, $line[3];
    $depth = int($line[4]);
    push @DepthsAHigh, $depth;
    $depth = int($line[5]);
    push @DepthsALow, $depth;
  }
  push @WinChrs, $line[0];
  push @AncDiffs, $line[-1];
  $WinChrs[-1] =~ s/L//;
  $WinChrs[-1] =~ s/R//;
}
close W;
$CMsX = $CMStopsX[-1];

my $TotalWindows = @CMStopsA + @CMStopsX;
my $AWindows = @CMStopsA;
my $XWindows = @CMStopsX;


#Get real peak data
my $i = 0;
my @RealPeakWindows = ();
my @RealPeakHeights = ();
my @RealPeakValleys = ();
open C, "<$ClusterFile" or die;
scalar (<C>);
while (<C>){
  chomp;
  last if m/^$/;
  @line = split;
  push @RealPeakValleys, $line[0];
  push @RealPeakWindows, $line[1];
  push @RealPeakHeights, $line[2];
}
close C;
$i = @RealPeakHeights;
print "Found data for $i peaks\n";

#get simulation results
my $j = 0;
my @SimPropVars = ();
my @SimPropVarAoA = ();
my @SimPeakHeightAoA = ();
my @SimPeakWindowAoA = ();
my @SimPeakHeights = ();
my @SimPeakWindows = ();
my @AllFiles = ();
my @SimPeakFiles = ();
opendir DIR, "." or die "couldn't open directory\n";
@AllFiles = readdir(DIR);
closedir DIR;
for ($i = 0; $i < @AllFiles; $i++){
  if (($AllFiles[$i] =~ m/$SimFile/) && ($AllFiles[$i] =~ m/txt/)){
    unless ($AllFiles[$i] =~ m/zip/){
      push @SimPeakFiles, $AllFiles[$i];
    }
  }
}
for ($i = 0; $i < @SimPeakFiles; $i++){
  open S, "<$SimPeakFiles[$i]" or die;
  scalar (<S>);
  while (<S>){
    chomp;
    last if m/^$/;
    @line = split;
    @SimPropVars = ();
    @SimPeakHeights = ();
    @SimPeakWindows = ();
    for ($j = 0; $j < @line; $j = $j + 3){
      push @SimPropVars, $line[$j];
      push @SimPeakHeights, $line[$j+1];
      push @SimPeakWindows, $line[$j+2];
    }
    push @SimPeakHeightAoA, [ @SimPeakHeights ];
    push @SimPeakWindowAoA, [ @SimPeakWindows ];
    push @SimPropVarAoA, [ @SimPropVars ];
  }
  close S;
}
$i = @SimPeakHeightAoA;
print "Found data for $i simulated replicates\n";
die if ($i < 100);

if (@RealPeakHeights != @{$SimPeakHeightAoA[0]}){
  $j = @{$SimPeakHeightAoA[0]};
  print "$i != $j\n";
  die;
}

#Retain sim peak data only in cases where all cluster peaks are within tolerance of real peak heights
my @GoodPropVars = ();
my @GoodPeakWindows = ();
my @GoodSecDevs = ();
my @TotalSimsFit = ();
my @PropVar5 = ();
my @PropVar95 = ();
my @SecDev95 = ();
my @PropVar50 = ();
my @AdjustedSecValleyHeights = ();
my @AdjustedSecValleyWindows = ();
my @SecPeakPValues = ();
my @SecPeakDeviations = ();
my @LeftBounds = ();
my @RightBounds = ();
my $k = 0;
my $l = 0;
my $ClusterStart = 0;
my $ClusterStop = 0;
my $SecDev = 0;
my $AdjustedValleyHeight = 0;
my $AdjustedValleyWindow = 0;
my $LowestValleyAwayHeight = 0;
my $LowestValleyAwayWindow = 0;
my $HighThreshRank = 0;
my $LowThreshRank = 0;
my $PValue = 0;
my $half = 0;
my $MedianValue = 0;

for ($i = 0; $i < @SimPeakHeightAoA; $i++){
  for ($j = 0; $j < @RealPeakHeights; $j++){
    if ((abs($RealPeakHeights[$j] - $SimPeakHeightAoA[$i][$j])) > $tolerance){
      splice @SimPeakHeightAoA, $i, 1;
      splice @SimPeakWindowAoA, $i, 1;
      splice @SimPropVarAoA, $i, 1;
      $i--;
      last;
    }
  }
}

$i = @SimPeakHeightAoA;
print "Kept data for $i simulated replicates\n";


for ($i = 0; $i < @RealPeakHeights; $i++){
  @GoodPropVars = ();
  @GoodPeakWindows = ();
  for ($j = 0; $j < @SimPeakHeightAoA; $j++){
    push @GoodPropVars, $SimPropVarAoA[$j][$i];
    push @GoodPeakWindows, $SimPeakWindowAoA[$j][$i];
  }
  $j = @GoodPropVars;
  push @TotalSimsFit, $j;

#Obtain desired critical values for each peak  
  @GoodPropVars = sort {$a <=> $b} @GoodPropVars;
  @GoodPeakWindows = sort {$a <=> $b} @GoodPeakWindows;
  $HighThreshRank = int(@GoodPropVars * 0.95) + 1;
  $LowThreshRank = int(@GoodPropVars * 0.05);
  push @PropVar5, $GoodPropVars[$LowThreshRank];
  push @PropVar95, $GoodPropVars[$HighThreshRank];
  push @LeftBounds, $GoodPeakWindows[$LowThreshRank]; 
  push @RightBounds, $GoodPeakWindows[$HighThreshRank];
  if ((@GoodPropVars % 2) == 0){
    $half = @GoodPropVars / 2;
    $MedianValue = ($GoodPropVars[$half] + $GoodPropVars[$half-1]) / 2;
    push @PropVar50, $MedianValue;
  }
  else{
    $half = int(@GoodPropVars / 2);
    $MedianValue = $GoodPropVars[$half];
    push @PropVar50, $MedianValue;
  }
  print "Evaluated thresholds for peak $i\n";
}

#zip up the sim files to get them out of the way
my $cmd = 'zip ' . $SimFile . '.zip';
for ($i = 0; $i < @SimPeakFiles; $i++){
  $cmd = $cmd . ' ' . $SimPeakFiles[$i];
}
system($cmd);
for ($i = 0; $i < @SimPeakFiles; $i++){
  $cmd = 'rm ' . $SimPeakFiles[$i];
  system($cmd);
}

#Output
open O, ">$OutputFile";
print O "RealPeakWindow\tRealPeakHeight\tMatchingSims\tPropVar50th\tPropVar5th\tPropVar95th\tLeftWinBound\tRightWinBound\n";
for ($i = 0; $i < @RealPeakHeights; $i++){
  print O "$RealPeakWindows[$i]\t$RealPeakHeights[$i]\t$TotalSimsFit[$i]\t$PropVar50[$i]\t$PropVar5[$i]\t$PropVar95[$i]\t$LeftBounds[$i]\t$RightBounds[$i]\n";
}
close O;
