#!/usr/bin/perl -w
use strict;
use cmdline_args;

if (!$ARGV[0] || $ARGV[0] eq "help"){
    print <<EOF;

This script reads in primary and secondary peak data, then summarizes QTL clusters for subsequent simulation.

It requires the following arguments:

-f [file name stem]

There must be a primary peak file from sibsam_1locus_analysis with this stem followed by "_primpeaks_fullresults.txt"
There must be a secondary peak file from sibsam_1locus_analysis with this stem followed by "_secpeaks_pvalues.txt"
Output files will have this stem appended by "_cluster1.txt" and so on.

EOF
exit 1;
}

my $FileStem = "";
my %args;
$args{'-f'} = [ \$FileStem, 1];

cmdline_args::get_options(\%args, \@ARGV);

my $PThresh = 0.05;  #only include secondary peaks at or below this P value in clusters

my $PrimPeakFile = $FileStem . '_primpeaks_fullresults.txt';
my $SecPeakFile = $FileStem . '_secpeaks_pvalues.txt';

#Get real peak data
my @line = ();
my @PrimPeakAoA = ();
open P, "<$PrimPeakFile" or die;
scalar (<P>);
while (<P>){
  chomp;
  last if m/^$/;
  @line = split;
  push @PrimPeakAoA, [ @line ];
}
close P;

my @SecPeakAoA = ();
open C, "<$SecPeakFile" or die;
scalar (<C>);
while (<C>){
  chomp;
  last if m/^$/;
  @line = split;
  push @SecPeakAoA, [ @line ];
}
close C;

#declarations
my $i = 0;
my $j = 0;
my $k = 0;
my $OutputFile = "";
my @ClusterPeakWindows = ();
my @ClusterPeakHeights = ();
my @ClusterValleyWindows = ();

#remove nonsignificant secondary peaks
for ($i = 0; $i < @SecPeakAoA; $i++){
  if ($SecPeakAoA[$i][-1] > $PThresh){
    splice @SecPeakAoA, $i, 1;
    $i--;
  }
}

#build some arrays for secondary peaks associated with a single primary peak
for ($i = 0; $i < @SecPeakAoA; $i++){
  if ($SecPeakAoA[$i][1] eq 'L'){
    unshift @ClusterPeakWindows, $SecPeakAoA[$i][2];
    unshift @ClusterPeakHeights, $SecPeakAoA[$i][3];
    unshift @ClusterValleyWindows, $SecPeakAoA[$i][4];
  }
  else{
    push @ClusterPeakWindows, $SecPeakAoA[$i][2];
    push @ClusterPeakHeights, $SecPeakAoA[$i][3];
    push @ClusterValleyWindows, $SecPeakAoA[$i][4];
  }
  
#if this is the last secondary peak in the cluster, add primary peak data and output to file
  if (($i == (@SecPeakAoA - 1)) || ($SecPeakAoA[$i+1][0] != $SecPeakAoA[$i][0])){
    for ($j = 0; $j < @PrimPeakAoA; $j++){
      next if ($PrimPeakAoA[$j][0] < $SecPeakAoA[$i][0]);
      for ($k = 0; $k < @ClusterPeakWindows; $k++){
	if ($PrimPeakAoA[$j][1] < $ClusterPeakWindows[$k]){
	  splice @ClusterPeakWindows, $k, 0, $PrimPeakAoA[$j][1];
	  splice @ClusterPeakHeights, $k, 0, $PrimPeakAoA[$j][2];
	  unshift @ClusterValleyWindows, 0;
	  last;
	}
	elsif ($k == (@ClusterPeakWindows - 1)){
	  push @ClusterPeakWindows, $PrimPeakAoA[$j][1];
	  push @ClusterPeakHeights, $PrimPeakAoA[$j][2];
	  unshift @ClusterValleyWindows, 0;
	  last;
	}
      }
      last;
    }
    $OutputFile = $FileStem . '_cluster' . $SecPeakAoA[$i][0] . '.txt';
    open O, ">$OutputFile";
    print O "PriorValleyWindow\tPeakWindow\tPeakHeight\n";
    for ($k = 0; $k < @ClusterPeakWindows; $k++){
      print O "$ClusterValleyWindows[$k]\t$ClusterPeakWindows[$k]\t$ClusterPeakHeights[$k]\n";
    }
    close O;
    print "Made a cluster file for primary peak ID $SecPeakAoA[$i][0]\n";
    @ClusterValleyWindows = ();
    @ClusterPeakWindows = ();
    @ClusterPeakHeights = ();
  }

}
	





