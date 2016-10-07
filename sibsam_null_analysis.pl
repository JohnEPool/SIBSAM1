#!/usr/bin/perl -w
use strict;
use cmdline_args;


if (!$ARGV[0] || $ARGV[0] eq "help"){
    print <<EOF;

For each primary peak in the empirical data, calculate true positive probability based on frequency of occurrence of that height or greater in real vs. null simulations.

It requires the following arguments:

-f [file name stem]

There must be one or more simulation files from sibsam_null.pl with the file stem followed by "_nullsims*" (where * is any characters)
There must be a window file with the above stem followed by "_windows.txt"
Output files will have the file stem followed by "_primpeaks_pvalues.txt" and "_secpeaks.txt"

EOF
exit 1;
}

#Get command line arguments
my $FileStem = "";
my %args;
$args{'-f'} = [ \$FileStem, 1];

cmdline_args::get_options(\%args, \@ARGV);

#my $RealPeakFile = $FileStem . '_primpeaks.txt';
my $SimFile = $FileStem . '_nullsims';
my $OutputFile = $FileStem . '_primpeaks_pvalues.txt';

my $PThresh = 0.05;  #focus on peaks below this P value
my $PrimThresh = 0.1;  #test all primary peaks at least this large
my $SecThresh = 0.1;  #identify secondary peaks with a deviation at least this large
my $SmoothEachSide = 4;  #include this many windows to each side of focal window in the weighted average that replaces raw window values of ancestry difference





#get real data primary peaks from input file
#open R, "<$RealPeakFile" or die;
#scalar (<R>);
#while (<R>){
#  chomp;
#  last if m/^$/;
#  @line = split;
#  push @RealPeakWindows, $line[1];
#  push @RealPeakHeights, $line[2];
#}
#close R;

#read in ancestry differences from window file
my $InputFile = $FileStem . '_windows.txt';
my $i = 0;
my @AncDiffs = ();
my @WinChrs = ();
my @line = ();
open I, "<$InputFile" or die;
scalar (<I>);
while (<I>){
  chomp;
  last if m/^$/;
  @line = split;
  $i = $line[-1] + 0;
  push @AncDiffs, $line[-1];
  push @WinChrs, $line[0];
  $WinChrs[-1] =~ s/L//;
  $WinChrs[-1] =~ s/R//;
}
close I;

#Smooth ancestry difference array
my $j = 0;
my $pos = 0;
my $weight = 0;
my $SmoothNum = 0;
my $SmoothDenom = 0;
my @SmoothDiffs = ();
for ($i = 0; $i < @AncDiffs; $i++){
  $SmoothNum = $AncDiffs[$i] * ($SmoothEachSide + 1);
  $SmoothDenom = $SmoothEachSide + 1;
  for ($j = 1; $j <= $SmoothEachSide + 0.5; $j++){
    $weight = ( ($SmoothEachSide + 1) - $j);
    $pos = $i - $j;
    if ((($i - $j) >= 0) && ($WinChrs[$i] eq $WinChrs[$pos])){
      $SmoothDenom += $weight;
      $SmoothNum += ($AncDiffs[$pos] * $weight);
    }
    $pos = $i + $j;
    if ((($i + $j) < @AncDiffs) && ($WinChrs[$i] eq $WinChrs[$pos])){
      $SmoothDenom += $weight;
      $SmoothNum += ($AncDiffs[$pos] * $weight);
    }
  }
  $SmoothNum = $SmoothNum / $SmoothDenom;
  push @SmoothDiffs, $SmoothNum
}
@AncDiffs = @SmoothDiffs;

#Identify all primary peaks with anc diff above an arbitrary threshold, separated by local valleys that reach down to 0 or below
my $PeakNow = 0;
my $MaxDiff = 0;
my $MaxWin = -1;
my @RealPeakHeights = ();
my @RealPeakWindows = ();
for ($i = 0; $i < @AncDiffs; $i++){
  if ($PeakNow == 0){
    if ($AncDiffs[$i] > $PrimThresh){
      $PeakNow = 1;
      $MaxDiff = $AncDiffs[$i];
      $MaxWin = $i;
    }
    next;
  }
  else{
    if ($AncDiffs[$i] > $MaxDiff){
      $MaxDiff = $AncDiffs[$i];
      $MaxWin = $i;
    }
    if (($AncDiffs[$i] < 0) || ($i == (@AncDiffs - 1)) || ($WinChrs[$i] ne $WinChrs[$i+1])){
      push @RealPeakHeights, $MaxDiff;
      push @RealPeakWindows, $MaxWin;
      $MaxDiff = 0;
      $MaxWin = -1;
      $PeakNow = 0;
    }
  }
}
if (($PeakNow == 1) && ($MaxDiff > $PrimThresh)){
  push @RealPeakHeights, $MaxDiff;
  push @RealPeakWindows, $MaxWin;
}
$i = @RealPeakHeights;
print "Found $i primary peaks above initial threshold for analysis\n";
print @RealPeakHeights;

#get sim peaks from input file(s)
my @SimPeakAoA = ();
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
    shift @line;
    shift @line;
    push @SimPeakAoA, [ @line ];
  }
  close S;
}

#for each real peak, test how often peaks that size or greater occur in real vs sim data
my $k = 0;
my $RealPeaksThisHigh = 0;
my $SimPeaksThisHigh = 0;
my $FalsePosProb = 0;
my $LowestSignifPeakHeight = 1;
my @FalsePosProbs = ();
for ($i = 0; $i < @RealPeakHeights; $i++){
  $RealPeaksThisHigh = 0;
  for ($j = 0; $j < @RealPeakHeights; $j++){
    if (($i==$j) || ($RealPeakHeights[$j] >= $RealPeakHeights[$i])){
      $RealPeaksThisHigh++;
    }
  }
  $SimPeaksThisHigh = 0;
  for ($j = 0; $j < @SimPeakAoA; $j++){
    for ($k = 0; $k < @{$SimPeakAoA[$j]}; $k++){
      if ($SimPeakAoA[$j][$k] >= $RealPeakHeights[$i]){
	$SimPeaksThisHigh++;
      }
    }
  }
  $SimPeaksThisHigh = $SimPeaksThisHigh / @SimPeakAoA;
  if ($SimPeaksThisHigh == 0){
    $FalsePosProb = 0;
  }
  else{
    $FalsePosProb = 1 - (($RealPeaksThisHigh / $SimPeaksThisHigh) - 1) / ($RealPeaksThisHigh / $SimPeaksThisHigh);
    if ($FalsePosProb >= 1){
      $FalsePosProb = 1;
    }
  }
  push @FalsePosProbs, $FalsePosProb;
  if (($FalsePosProb < $PThresh) && ($RealPeakHeights[$i] < $LowestSignifPeakHeight)){
    $LowestSignifPeakHeight = $RealPeakHeights[$i];
  } 
}

#identify the critical value for primary peak height (used later to exclude shorter secondary peaks)
$RealPeaksThisHigh = 0;
my @SortedPeakHeights = sort { $b <=> $a } @RealPeakHeights;
for ($i = 0; $i < @SortedPeakHeights; $i++){
  if ($SortedPeakHeights[$i] >= $LowestSignifPeakHeight){
    $RealPeaksThisHigh++;
    shift @SortedPeakHeights;
  }
  else{
    last;
  }
}

my $HeightThresh = 0;
if ($RealPeaksThisHigh == 0){
  print "No primary peaks achieved significance\n";
}  
else{
$HeightThresh = (int(($LowestSignifPeakHeight + 0.0005) * 1000)) / 1000;
while ($HeightThresh > 0){
  $SimPeaksThisHigh = 0;
  for ($j = 0; $j < @SimPeakAoA; $j++){
    for ($k = 0; $k < @{$SimPeakAoA[$j]}; $k++){
      if ($SimPeakAoA[$j][$k] >= $HeightThresh){
	$SimPeaksThisHigh++;
      }
    }
  }
  $SimPeaksThisHigh = $SimPeaksThisHigh / @SimPeakAoA;
  if ($SimPeaksThisHigh == 0){
    $FalsePosProb = 0;
  }
  else{
    $FalsePosProb = 1 - (($RealPeaksThisHigh / $SimPeaksThisHigh) - 1) / ($RealPeaksThisHigh / $SimPeaksThisHigh);
    if ($FalsePosProb >= 1){
      $FalsePosProb = 1;
    }
  }
  if ($FalsePosProb > 0.05){
    $HeightThresh = $HeightThresh + 0.001;
    last;
  }
  $HeightThresh = $HeightThresh - 0.001;
}
}

#zip up the null sim files to get them out of the way
my $cmd = 'zip ' . $SimFile . '.zip';
for ($i = 0; $i < @SimPeakFiles; $i++){
  $cmd = $cmd . ' ' . $SimPeakFiles[$i];
}
system($cmd);
for ($i = 0; $i < @SimPeakFiles; $i++){
  $cmd = 'rm ' . $SimPeakFiles[$i];
  system($cmd);
}

#prim peak p value output
open O, ">$OutputFile";
print O "PeakID\tPrimPeakWin\tPrimPeakHeight\tPValue\n";
for ($i = 0; $i < @RealPeakWindows; $i++){
  print O "$i\t$RealPeakWindows[$i]\t$RealPeakHeights[$i]\t$FalsePosProbs[$i]\n";
}
close O;

#SECONDARY PEAK IDENTIFICATION
my $MinWin = 0;
my $MinSincePeak = 0;
my $PeakHeight = 0;
my @SecPeakHeights = ();
my @SecPeakWindows = ();
my @SecPeakLeftRights = ();
my @SecValleyHeights = ();
my @SecValleyWindows = ();
my @SecPeakPrimIDs = ();
my @AllSecPeakHeights = ();
my @AllSecPeakWindows = ();
my @AllSecPeakLeftRights = ();
my @AllSecValleyHeights = ();
my @AllSecValleyWindows = ();
my @AllSecPeakPrimIDs = ();

#For each SIGNIFICANT primary peak, identify any secondary peaks with secondary amplitude above a preliminary threshold.  First moving leftward...
for ($i = 0; $i < @RealPeakWindows; $i++){
  next if ($FalsePosProbs[$i] > $PThresh);
  $PeakHeight = $RealPeakHeights[$i];
  $MinSincePeak = $PeakHeight;
  $MaxDiff = 0;
  $MaxWin = -1;
  $MinWin = -1;
  $PeakNow = 0;
  for ($j = ($RealPeakWindows[$i] - 1); $j >= 0; $j--){
    if (($j == 0) || ($WinChrs[$j] ne $WinChrs[$j-1]) || ($AncDiffs[$j] <= 0)){
      if ($PeakNow == 1){
	if ($AncDiffs[$j] > $MaxDiff){
	  $MaxWin = $j;
	  $MaxDiff = $AncDiffs[$j];
	}
	push @SecPeakHeights, $MaxDiff;
	push @SecPeakWindows, $MaxWin;
	push @SecValleyHeights, $MinSincePeak;
	push @SecValleyWindows, $MinWin;
	push @SecPeakLeftRights, 'L';
	push @SecPeakPrimIDs, $i;
      }
      last;
    }
    elsif ($PeakNow == 0){
      if ($AncDiffs[$j] < $MinSincePeak){
	$MinSincePeak = $AncDiffs[$j];
	$MinWin = $j;
	next;
      }
      elsif (($AncDiffs[$j] >= ($MinSincePeak + $SecThresh)) && ($AncDiffs[$j] > $HeightThresh)){
	$PeakNow = 1;
	$MaxWin = $j;
	$MaxDiff = $AncDiffs[$j];
      }
    }
    else{
      if ($AncDiffs[$j] > $MaxDiff){
	$MaxWin = $j;
	$MaxDiff = $AncDiffs[$j];
      }
      elsif($AncDiffs[$j] < ($MaxDiff - $SecThresh)){
	push @SecPeakWindows, $MaxWin;
	push @SecPeakHeights, $MaxDiff;
	push @SecValleyHeights, $MinSincePeak;
	push @SecValleyWindows, $MinWin;
	push @SecPeakLeftRights, 'L';
	push @SecPeakPrimIDs, $i;
	$MinSincePeak = $AncDiffs[$j];
	$MaxDiff = 0;
	$MaxWin = -1;
	$MinWin = $j;
	$PeakNow = 0;
      }
    }
  }
  @SecPeakWindows = reverse @SecPeakWindows;
  @SecPeakHeights = reverse @SecPeakHeights;
  @SecValleyHeights = reverse @SecValleyHeights;
  @SecValleyWindows = reverse @SecValleyWindows;

#then moving rightward.
  $MinSincePeak = $PeakHeight;
  $MaxDiff = 0;
  $MaxWin = -1;
  $MinWin = -1;
  $PeakNow = 0;
  for ($j = ($RealPeakWindows[$i] + 1); $j < @AncDiffs; $j++){
    if (($j == (@AncDiffs - 1)) || ($WinChrs[$j] ne $WinChrs[$j+1]) || ($AncDiffs[$j] <= 0)){
      if ($PeakNow == 1){
	if ($AncDiffs[$j] > $MaxDiff){
	  $MaxWin = $j;
	  $MaxDiff = $AncDiffs[$j];
	}
	push @SecPeakHeights, $MaxDiff;
	push @SecPeakWindows, $MaxWin;
	push @SecValleyHeights, $MinSincePeak;
	push @SecValleyWindows, $MinWin;
	push @SecPeakLeftRights, 'R';
	push @SecPeakPrimIDs, $i;
      }
      last;
    }
    elsif ($PeakNow == 0){
      if ($AncDiffs[$j] < $MinSincePeak){
	$MinSincePeak = $AncDiffs[$j];
	$MinWin = $j;
	next;
      }
      elsif (($AncDiffs[$j] >= ($MinSincePeak + $SecThresh))  && ($AncDiffs[$j] > $HeightThresh)){
	$PeakNow = 1;
	$MaxWin = $j;
	$MaxDiff = $AncDiffs[$j];
      }
    }
    else{
      if ($AncDiffs[$j] > $MaxDiff){
	$MaxWin = $j;
	$MaxDiff = $AncDiffs[$j];
      }
      elsif($AncDiffs[$j] < ($MaxDiff - $SecThresh)){
	push @SecPeakWindows, $MaxWin;
	push @SecPeakHeights, $MaxDiff;
	push @SecValleyHeights, $MinSincePeak;
	push @SecValleyWindows, $MinWin;
	push @SecPeakLeftRights, 'R';
	push @SecPeakPrimIDs, $i;
	$MinSincePeak = $AncDiffs[$j];
	$MaxDiff = 0;
	$MaxWin = -1;
	$MinWin = $j;
	$PeakNow = 0;
      }
    }
  }
  push @AllSecPeakWindows, @SecPeakWindows;
  push @AllSecPeakHeights, @SecPeakHeights;
  push @AllSecValleyHeights, @SecValleyHeights;
  push @AllSecValleyWindows, @SecValleyWindows;
  push @AllSecPeakLeftRights, @SecPeakLeftRights;
  push @AllSecPeakPrimIDs, @SecPeakPrimIDs;
  @SecValleyHeights = ();
  @SecValleyWindows = ();
  @SecPeakHeights = ();
  @SecPeakWindows = ();
  @SecPeakLeftRights = ();
  @SecPeakPrimIDs = ();
}

my $file = $FileStem . '_secpeaks.txt';
open S, ">$file";
print S "PrimPeakID\tSecPeakLeftRight\tSecPeakWin\tSecPeakHeight\tSecValleyWin\tSecValleyHeight\tHeightThresh:\t$HeightThresh\n";
for ($i = 0; $i < @AllSecPeakPrimIDs; $i++){
  print S "$AllSecPeakPrimIDs[$i]\t$AllSecPeakLeftRights[$i]\t$AllSecPeakWindows[$i]\t$AllSecPeakHeights[$i]\t$AllSecValleyWindows[$i]\t$AllSecValleyHeights[$i]\n";
}
close S;
