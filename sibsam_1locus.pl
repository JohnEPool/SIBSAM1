#!/usr/bin/perl -w
use strict;
use cmdline_args;

if (!$ARGV[0] || $ARGV[0] eq "help"){
    print <<EOF;

This is a forward simulation program designed to test performance of bulk segregant analysis in the Drosophila genome.

It simulates with 1 randomly located causative locus, and returns primary and secondary peak data.

It requires the following arguments:
-g [number of generations of crosses before sampling (g>=2)]
-n [number of individuals (females plus males) in the generations before sampling]
-l [number of individuals used for phenotypic selection in the last generation]
-s [proportion of individuals in each phenotypic tail selected for sequencing (0 < s <= 0.5)]
-c [number of independent crosses to sum results across]   <DISABLED IN THIS VERSION>
-r [number of replicates to simulate]
-f [file name stem]
-b [batch number, optional to allow parallel processing]

There must be a window file with the above stem followed by "_windows.txt".
There must also be a file from sibsam_null_analysis.pl with this stem followed by "_primpeaks_pvalues.txt".
Output files will have this stem appended by "_1locsims.txt", for example.

EOF
exit 1;
}

#Get command line arguments
my $gens = "";
my $F2N = "";
my $LastNf = "";
my $SelectionProp = "";
my $BatchNum = "";
my $NCrosses = 1;
my $TrueReps = "";
my $FileStem = "";

my %args;
$args{'-g'} = [ \$gens, 1];
$args{'-n'} = [ \$F2N, 1];
$args{'-l'} = [ \$LastNf, 1];
$args{'-s'} = [ \$SelectionProp, 1];
$args{'-b'} = [ \$BatchNum, 1];
$args{'-c'} = [ \$NCrosses, 1];
$args{'-r'} = [ \$TrueReps, 1];
$args{'-f'} = [ \$FileStem, 1];

cmdline_args::get_options(\%args, \@ARGV);

my $WinFile = $FileStem . '_windows.txt';
my $PeakFile = $FileStem . '_primpeaks_pvalues.txt';
my $OutputFile = $FileStem . '_1locsims' . $BatchNum . '.txt';

if ($NCrosses > 1){
  die "Use sibsam_1locus_multicross.pl instead\n";
}

#my $QTLsPerReplicate = 5;  #number of loci simulated for separate/independent pool selections from each cage simulation.  number each for X and autosomes.
#if ($NCrosses > 1){        #can only do 1 locus per cage if we're combining multiple cages into 1 analysis
#  $QTLsPerReplicate = 1;
#}
my $MaxLocusProp = 1;   
my $MinLocusProp = 0;  #determines how much of phenotypic variance a locus's additive effect can explain
my $SmoothEachSide = 4;  #include this many windows on each side of the focal window in a weighted moving average (e.g. if = 4, weights are 1,2,3,4,5,4,3,2,1)
my $PThresh = 0.05;

my $F1N = $F2N;  #Not currently varying the number of F1's independently
my $LastN = $LastNf * 2;
my $replicates = $TrueReps * $NCrosses; #how many replicates to actually simulate
#my $multiple = $NCrosses;  #combine signals from this many independent crosses


my $i = 0;
my $pos = 0;
my $depth = 0;

#Get cM and depth information from window files
my $CMsX = 0;
my $CMsA = 0;
my @line = ();
my @DepthsXHigh = ();
my @DepthsAHigh = ();
my @DepthsXLow = ();
my @DepthsALow = ();
my @CMStopsX = ();
my @CMStopsA = ();
my @AutoCMBreaks = ();
my @WinChrsA = ();
my @WinChrsX = ();

open W, "<$WinFile" or die;
scalar (<W>);
while (<W>){
  chomp;
  last if m/^$/;
  @line = split;
  if ($line[0] =~ m/X/){
    push @WinChrsX, $line[0];
    push @CMStopsX, $line[3];
    $depth = int($line[4]);
    push @DepthsXHigh, $depth;
    $depth = int($line[5]);
    push @DepthsXLow, $depth;
  }
  else{
    push @WinChrsA, $line[0];
    push @CMStopsA, $line[3];
    $depth = int($line[4]);
    push @DepthsAHigh, $depth;
    $depth = int($line[5]);
    push @DepthsALow, $depth;
  }
}
close W;
$CMsX = $CMStopsX[-1];

my $TotalWindows = @CMStopsA + @CMStopsX;
my $AWindows = @CMStopsA;
my $XWindows = @CMStopsX;

#Look for chromosome breaks in autosomal window data (cases where this window's CM stop position is smaller than the last window's)
for ($i = 0; $i < (@CMStopsA - 1); $i++){
  if ($CMStopsA[$i+1] < $CMStopsA[$i]){
    $CMStopsA[$i] += $CMsA;
    $CMsA = $CMStopsA[$i];
    push @AutoCMBreaks, $CMsA;
  }
  else{
    $CMStopsA[$i] += $CMsA;
  }
}
$CMStopsA[-1] += $CMsA;
$CMsA = $CMStopsA[-1];
$i = @AutoCMBreaks + 1;
print "Found window data for $XWindows X chromosome windows ($CMsX cMs) and $AWindows windows on $i autosomes ($CMsA total cMs)\n";

#set up arrays
my @PositionsX = ();
my @AncDiffsX = ();
my @PositionsA = ();
my @AncDiffsA = ();
for ($i = 0; $i < @CMStopsX; $i++){
  $pos = $CMStopsX[$i] / $CMsX;
  push @PositionsX, $pos;
  push @AncDiffsX, 0;
}
for ($i = 0; $i < @CMStopsA; $i++){
  $pos = $CMStopsA[$i] / $CMsA;
  push @PositionsA, $pos;
  push @AncDiffsA, 0;
}

#open peak file and get peak locations to simulate
#Set up locus positions and effects
my @LocusPosX = ();
my @LocusPosA = ();
my @AllQTLsX = ();
my @AllQTLsA = ();
open P, "<$PeakFile" or die;
scalar (<P>);
while (<P>){
  chomp;
  last if m/^$/;
  @line = split;
  if ($line[-1] < $PThresh){
    if ($line[1] < @PositionsX){
      push @AllQTLsX, $line[1];
    }
    else{
      $line[1] -= @PositionsX;
      push @AllQTLsA, $line[1];
    }
  }
}
close P; 

$i = @AllQTLsX;
print "Found $i significant QTLs on the X:  ";
for ($i = 0; $i < @AllQTLsX; $i++){
  print "$AllQTLsX[$i]  "
}
print "\n";
$i = @AllQTLsA;
print "Found $i significant QTLs on the autosomes:  ";
for ($i = 0; $i < @AllQTLsA; $i++){
  print "$AllQTLsA[$i]  "
}
print "\n";


#Factorial subroutine
my $A = 0;
my $B = 0;
my $FactIn = 0;
sub factorial{
    $A = 1;
    $B = 1;
    while ($A <= $FactIn) {
	$A *= $B;
	$B++;
    }
    return $A;
}

#Probability of a given number of recombination events (Poisson) for each chromosome
#(allowing up to 10 recombination events per chromosome per generation)
#Cumulative probabilities are given for each number, for easy use with random numbers later.
my $TotalLength = $CMsX;
my $Pr = 0;
my $PrSum = 0;
my @RecombProbsX = ();
my $ExpRecEvents = 0.01 * $TotalLength;
for ($i = 0; $i < 11; $i++){
  $FactIn = $i;
  $Pr = ((exp(-1 * $ExpRecEvents)) * ($ExpRecEvents ** $i)) / &factorial;
  $PrSum = $PrSum + $Pr;
  push @RecombProbsX, $PrSum;
}
my @RecombProbsA = ();
$TotalLength = $CMsA;
$ExpRecEvents = 0.01 * $TotalLength;
$PrSum = 0;
for ($i = 0; $i < 11; $i++){
  $FactIn = $i;
  $Pr = ((exp(-1 * $ExpRecEvents)) * ($ExpRecEvents ** $i)) / &factorial;
  $PrSum = $PrSum + $Pr;
  push @RecombProbsA, $PrSum;
}

for ($i = 0; $i < @AutoCMBreaks; $i++){
  $AutoCMBreaks[$i] = $AutoCMBreaks[$i] / $CMsA;
}

#build breakpoint matrix for F1 individuals
my @OldFemaleXAoA = ();
my @OldFemaleAAoA = ();
my @OldMaleXAoA = ();
my @OldMaleAAoA = ();
my @NewFemaleXAoA = ();
my @NewFemaleAAoA = ();
my @NewMaleXAoA = ();
my @NewMaleAAoA = ();

my @IndBreaks1 = (); #chromosomes from from 1st parental line
my @IndBreaks2 = (0); #chromosomes from from 2nd parental line

#normal subroutine
sub gaussian_rand {
    my ($u1, $u2);  # uniformly distributed random numbers
    my $w;          # variance, then a weight
    my ($g1, $g2);  # gaussian-distributed numbers

    do {
        $u1 = 2 * rand() - 1;
        $u2 = 2 * rand() - 1;
        $w = $u1*$u1 + $u2*$u2;
    } while ( $w >= 1 );

    $w = sqrt( (-2 * log($w))  / $w );
    $g1 = $u2 * $w;
    return $g1;
}

#more declarations
my $g = 0;
my $j = 0;
my $k = 0;
my $r = 0;
my $NEachSex = 0;
my $random = 0;
my $random2 = 0;
my $parents = 0;
my $PAllele1 = 0;
my $PAllele2 = 0;
my $AncBeforeBreak = 0;
my $AncAfterBreak = 0;
my $CurrentAllele = 0;
my $LastBreak = 0;
my $AncProp = 0;
my $EnvVar = 0;
my $anc = 0;
my $allele = 0;
my $break = 0;
my $LowerThresh = 0;
my $UpperThresh = 0;
my $AncDiff = 0;
my $MaxDiff = 0;
my $MaxDiffA = 0;
my $MaxDiffX = 0;
my $MaxPos = 0;
my $RegionStart = 0;
my $RegionStop = 0;
my $PowerCount = 0;
my $SmoothNum = 0;
my $SmoothDenom = 0;
my $weight = 0;
my $PeakNow = 0;
my $LocusProp = 0;
my $LocusStrength = 0;
my $TargetWin = 0;
my $TargetDiff = 0;
my $MaxWin = 0;
my $q = 0;
my $VertDev = 0;
my $MaxDev = 0;
my $PeakWin = 0;
my $MinSincePeak = 0;
my $l = 0;
my $AdjustedValley = 0;
my $ClusterStart = 0;
my $ClusterStop = 0;
my $SecDev = 0;
my @RecombPos = ();
my @AncBreaks = ();
my @A1Breaks = ();
my @A2Breaks = ();
my @ancestries = ();
my @BreakpointNumbers = ();
my @AncPropsX = ();
my @AncPropsA = ();
my @AncestryTracker = ();
my @BreakpointTracker = ();
my @AncXAoA = ();
my @AncAAoA = ();
my @blank = ();
my @phenotypes = ();
my @sorted = ();
my @array = ();
my @AncXHighAoA = ();
my @AncAHighAoA = ();
my @AncXLowAoA = ();
my @AncALowAoA = ();
my @AncPropsXHigh = ();
my @AncPropsXLow = ();
my @AncPropsAHigh = ();
my @AncPropsALow = ();
my @PeakDist = ();
my @PeakDistAoA = ();
my @AllPeakDists = ();
my @MaxDiffs = ();
my @MaxDiffAs = ();
my @MaxDiffXs = ();
my @AncDiffsAoA = ();
my @AllPeaks = ();
my @AllPeaksAoA = ();
my @PropVar = ();
my @PeakHeights = ();
my @PeakWindows = ();
my @SecDevs = ();
my @PeakHeightAoA = ();
my @PeakWindowAoA = ();
my @SecDevAoA = ();
my @SecPeakDeviations = ();
my @SmoothDiffs = ();
my @LocusStrengths = ();

#Replicate loop (not indented)
for ($r = 1; $r <= $replicates; $r++){

@OldFemaleXAoA = ();
@OldFemaleAAoA = ();
@OldMaleXAoA = ();
@OldMaleAAoA = ();

for ($i = 0; $i < ($F1N / 2); $i++){
  push @OldFemaleXAoA, [ @IndBreaks1 ];
  push @OldFemaleXAoA, [ @IndBreaks2 ];
  push @OldFemaleAAoA, [ @IndBreaks1 ];
  push @OldFemaleAAoA, [ @IndBreaks2 ];
  push @OldMaleAAoA, [ @IndBreaks1 ];
  push @OldMaleAAoA, [ @IndBreaks2 ];
}
for ($i = 0; $i < ($F1N / 4); $i++){
  push @OldMaleXAoA, [ @IndBreaks1 ];
  push @OldMaleXAoA, [ @IndBreaks2 ];
}
$NEachSex = $F2N / 2;

#For each individual in the next generation, choose parents, enact recombination from mother, transfer breakpoints
for ($g = 2; $g <= $gens; $g++){
  $parents = @OldFemaleXAoA / 2;
  if ($g == $gens){
    $NEachSex = $LastN / 2;
  }

#for females in the next generation, choose mother
  for ($i = 0; $i < $NEachSex; $i++){
    $random = int(rand($parents));
    $random2 = rand;
    if ($random2 < 0.5){
      $PAllele1 = $random * 2;
      $PAllele2 = $PAllele1 + 1;   
    }
    else{
      $PAllele2 = $random * 2;
      $PAllele1 = $PAllele2 + 1;   
    }

#for X chromosome, determine the number of recombination events (store that many random breakpoints) from maternal alleles
    $random = rand;
    @RecombPos = ();
    for ($j = 0; $j < @RecombProbsX; $j++){
      if ($random > $RecombProbsX[$j]){
	push @RecombPos, rand;
      }
      else{
	last;
      }
    }
    if (@RecombPos > 1){
      @RecombPos = sort @RecombPos;
    }
#in light of recombination positions, determine new ancestry along recombinant chromosome
    if (@RecombPos == 0){
      @AncBreaks = @{$OldFemaleXAoA[$PAllele1]};
      push @NewFemaleXAoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldFemaleXAoA[$PAllele1]};
      @A2Breaks = @{$OldFemaleXAoA[$PAllele2]};
      $CurrentAllele = 1;
      $LastBreak = -1;
      for ($j = 0; $j < @RecombPos; $j++){
	$AncBeforeBreak = 1;
	$AncAfterBreak = 1;
	if ($CurrentAllele % 2 == 1){
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	      $LastBreak = $A1Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	else{
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	      $LastBreak = $A2Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	if ( ($AncBeforeBreak % 2) != ($AncAfterBreak % 2) ){
	  push @AncBreaks, $RecombPos[$j];
	  $LastBreak = $RecombPos[$j];
	}
	$CurrentAllele++;
      }
      if ($CurrentAllele % 2 == 1){
	for ($k = 0; $k < @A1Breaks; $k++){
	  if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	    }
	}
      }
      else{
	for ($k = 0; $k < @A2Breaks; $k++){
	  if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	    }
	}
      }
      push @NewFemaleXAoA, [ @AncBreaks ];
    }

#for the merged autosome, determine the number of recombination events (store that many random breakpoints) from maternal alleles
    $random = rand;
    @RecombPos = ();
    for ($j = 0; $j < @RecombProbsA; $j++){
      if ($random > $RecombProbsA[$j]){
	push @RecombPos, rand;
      }
      else{
	last;
      }
    }
#each time there's a chromosome break, switch to the parent's other copy with 50% probability
    for ($j = 0; $j < @AutoCMBreaks; $j++){
      $random = rand;
      if ($random > 0.5){
	push @RecombPos, $AutoCMBreaks[$j];
      }
    }
    if (@RecombPos > 1){
      @RecombPos = sort @RecombPos;
    }
#in light of recombination positions, determine new ancestry along recombinant chromosome
    if (@RecombPos == 0){
      @AncBreaks = @{$OldFemaleAAoA[$PAllele1]};
      push @NewFemaleAAoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldFemaleAAoA[$PAllele1]};
      @A2Breaks = @{$OldFemaleAAoA[$PAllele2]};
      $CurrentAllele = 1;
      $LastBreak = -1;
      for ($j = 0; $j < @RecombPos; $j++){
	$AncBeforeBreak = 1;
	$AncAfterBreak = 1;
	if ($CurrentAllele % 2 == 1){
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	      $LastBreak = $A1Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	else{
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	      $LastBreak = $A2Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	if ( ($AncBeforeBreak % 2) != ($AncAfterBreak % 2) ){
	  push @AncBreaks, $RecombPos[$j];
	  $LastBreak = $RecombPos[$j];
	}
	$CurrentAllele++;
      }
      if ($CurrentAllele % 2 == 1){
	for ($k = 0; $k < @A1Breaks; $k++){
	  if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	    }
	}
      }
      else{
	for ($k = 0; $k < @A2Breaks; $k++){
	  if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	    }
	}
      }
      for ($k = 0; $k < (@AncBreaks-1); $k++){
	if ($AncBreaks[$k] == $AncBreaks[$k+1]){
	  splice (@AncBreaks, $k, 2);
	  $k--;
	}
      }
      push @NewFemaleAAoA, [ @AncBreaks ];
    }

#for females in the next generation, choose father
#transfer paternal chromosomes without recombination (except allow autosomal independent assortment)
    $PAllele1 = int(rand($parents));
    @AncBreaks = @{$OldMaleXAoA[$PAllele1]};
    push @NewFemaleXAoA, [ @AncBreaks ];
    $random = rand;
    if ($random > 0.5){
      $PAllele1 = $PAllele1 * 2;
      $PAllele2 = $PAllele1 + 1;
    }
    else{
      $PAllele2 = $PAllele1 * 2;
      $PAllele1 = $PAllele2 + 1;
    }
#each time there's a chromosome break, switch to the parent's other copy with 50% probability
    @RecombPos = ();
    for ($j = 0; $j < @AutoCMBreaks; $j++){
      $random = rand;
      if ($random > 0.5){
	push @RecombPos, $AutoCMBreaks[$j];
      }
    }
#in light of recombination positions, determine new ancestry along recombinant chromosome
    if (@RecombPos == 0){
      @AncBreaks = @{$OldMaleAAoA[$PAllele1]};
      push @NewFemaleAAoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldMaleAAoA[$PAllele1]};
      @A2Breaks = @{$OldMaleAAoA[$PAllele2]};
      $CurrentAllele = 1;
      $LastBreak = -1;
      for ($j = 0; $j < @RecombPos; $j++){
	$AncBeforeBreak = 1;
	$AncAfterBreak = 1;
	if ($CurrentAllele % 2 == 1){
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	      $LastBreak = $A1Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	else{
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	      $LastBreak = $A2Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	if ( ($AncBeforeBreak % 2) != ($AncAfterBreak % 2) ){
	  push @AncBreaks, $RecombPos[$j];
	  $LastBreak = $RecombPos[$j];
	}
	$CurrentAllele++;
      }
      if ($CurrentAllele % 2 == 1){
	for ($k = 0; $k < @A1Breaks; $k++){
	  if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	    }
	}
      }
      else{
	for ($k = 0; $k < @A2Breaks; $k++){
	  if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	    }
	}
      }
      for ($k = 0; $k < (@AncBreaks-1); $k++){
	if ($AncBreaks[$k] == $AncBreaks[$k+1]){
	  splice (@AncBreaks, $k, 2);
	  $k--;
	}
      }
      push @NewFemaleAAoA, [ @AncBreaks ];
    }
  }

#for males in the next generation, choose mother
  for ($i = 0; $i < $NEachSex; $i++){
    $random = int(rand($parents));
    $random2 = rand;
    if ($random2 < 0.5){
      $PAllele1 = $random * 2;
      $PAllele2 = $PAllele1 + 1;   
    }
    else{
      $PAllele2 = $random * 2;
      $PAllele1 = $PAllele2 + 1;   
    }

#for X chromosome, determine the number of recombination events (store that many random breakpoints) from maternal alleles
    $random = rand;
    @RecombPos = ();
    for ($j = 0; $j < @RecombProbsX; $j++){
      if ($random > $RecombProbsX[$j]){
	push @RecombPos, rand;
      }
      else{
	last;
      }
    }
    if (@RecombPos > 1){
      @RecombPos = sort @RecombPos;
    }
#in light of recombination positions, determine new ancestry along recombinant chromosome
    if (@RecombPos == 0){
      @AncBreaks = @{$OldFemaleXAoA[$PAllele1]};
      push @NewMaleXAoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldFemaleXAoA[$PAllele1]};
      @A2Breaks = @{$OldFemaleXAoA[$PAllele2]};
      $CurrentAllele = 1;
      $LastBreak = -1;
      for ($j = 0; $j < @RecombPos; $j++){
	$AncBeforeBreak = 1;
	$AncAfterBreak = 1;
	if ($CurrentAllele % 2 == 1){
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	      $LastBreak = $A1Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	else{
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	      $LastBreak = $A2Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	if ( ($AncBeforeBreak % 2) != ($AncAfterBreak % 2) ){
	  push @AncBreaks, $RecombPos[$j];
	  $LastBreak = $RecombPos[$j];
	}
	$CurrentAllele++;
      }
      if ($CurrentAllele % 2 == 1){
	for ($k = 0; $k < @A1Breaks; $k++){
	  if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	    }
	}
      }
      else{
	for ($k = 0; $k < @A2Breaks; $k++){
	  if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	    }
	}
      }
      push @NewMaleXAoA, [ @AncBreaks ];
    }

#for the merged autosome, determine the number of recombination events (store that many random breakpoints) from maternal alleles
    $random = rand;
    @RecombPos = ();
    for ($j = 0; $j < @RecombProbsA; $j++){
      if ($random > $RecombProbsA[$j]){
	push @RecombPos, rand;
      }
      else{
	last;
      }
    }
#each time there's a chromosome break, switch to the parent's other copy with 50% probability
    for ($j = 0; $j < @AutoCMBreaks; $j++){
      $random = rand;
      if ($random > 0.5){
	push @RecombPos, $AutoCMBreaks[$j];
      }
    }
    if (@RecombPos > 1){
      @RecombPos = sort @RecombPos;
    }
#in light of recombination positions, determine new ancestry along recombinant chromosome
    if (@RecombPos == 0){
      @AncBreaks = @{$OldFemaleAAoA[$PAllele1]};
      push @NewMaleAAoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldFemaleAAoA[$PAllele1]};
      @A2Breaks = @{$OldFemaleAAoA[$PAllele2]};
      $CurrentAllele = 1;
      $LastBreak = -1;
      for ($j = 0; $j < @RecombPos; $j++){
	$AncBeforeBreak = 1;
	$AncAfterBreak = 1;
	if ($CurrentAllele % 2 == 1){
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	      $LastBreak = $A1Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	else{
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	      $LastBreak = $A2Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	if ( ($AncBeforeBreak % 2) != ($AncAfterBreak % 2) ){
	  push @AncBreaks, $RecombPos[$j];
	  $LastBreak = $RecombPos[$j];
	}
	$CurrentAllele++;
      }
      if ($CurrentAllele % 2 == 1){
	for ($k = 0; $k < @A1Breaks; $k++){
	  if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	    }
	}
      }
      else{
	for ($k = 0; $k < @A2Breaks; $k++){
	  if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	    }
	}
      }
      for ($k = 0; $k < (@AncBreaks-1); $k++){
	if ($AncBreaks[$k] == $AncBreaks[$k+1]){
	  splice (@AncBreaks, $k, 2);
	  $k--;
	}
      }
      push @NewMaleAAoA, [ @AncBreaks ];
    }

#for males in the next generation, choose father
#transfer paternal chromosomes without recombination (except allow autosomal independent assortment)
    $PAllele1 = int(rand($parents));
    $random = rand;
    if ($random > 0.5){
      $PAllele1 = $PAllele1 * 2;
      $PAllele2 = $PAllele1 + 1;
    }
    else{
      $PAllele2 = $PAllele1 * 2;
      $PAllele1 = $PAllele2 + 1;
    }
#each time there's a chromosome break, switch to the parent's other copy with 50% probability
    @RecombPos = ();
    for ($j = 0; $j < @AutoCMBreaks; $j++){
      $random = rand;
      if ($random > 0.5){
	push @RecombPos, $AutoCMBreaks[$j];
      }
    }
#in light of recombination positions, determine new ancestry along recombinant chromosome
    if (@RecombPos == 0){
      @AncBreaks = @{$OldMaleAAoA[$PAllele1]};
      push @NewMaleAAoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldMaleAAoA[$PAllele1]};
      @A2Breaks = @{$OldMaleAAoA[$PAllele2]};
      $CurrentAllele = 1;
      $LastBreak = -1;
      for ($j = 0; $j < @RecombPos; $j++){
	$AncBeforeBreak = 1;
	$AncAfterBreak = 1;
	if ($CurrentAllele % 2 == 1){
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	      $LastBreak = $A1Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	else{
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	      $LastBreak = $A2Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	if ( ($AncBeforeBreak % 2) != ($AncAfterBreak % 2) ){
	  push @AncBreaks, $RecombPos[$j];
	  $LastBreak = $RecombPos[$j];
	}
	$CurrentAllele++;
      }
      if ($CurrentAllele % 2 == 1){
	for ($k = 0; $k < @A1Breaks; $k++){
	  if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	    }
	}
      }
      else{
	for ($k = 0; $k < @A2Breaks; $k++){
	  if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	    }
	}
      }
      for ($k = 0; $k < (@AncBreaks-1); $k++){
	if ($AncBreaks[$k] == $AncBreaks[$k+1]){
	  splice (@AncBreaks, $k, 2);
	  $k--;
	}
      }
      push @NewMaleAAoA, [ @AncBreaks ];
    }
  }

#transfer new ancestry breakpoint matrices to old matrices to prepare for next generation
  @OldFemaleXAoA = @NewFemaleXAoA;
  @NewFemaleXAoA = ();
  @OldFemaleAAoA = @NewFemaleAAoA;
  @NewFemaleAAoA = ();
  @OldMaleXAoA = @NewMaleXAoA;
  @NewMaleXAoA = ();
  @OldMaleAAoA = @NewMaleAAoA;
  @NewMaleAAoA = ();
}

#
#
#individual ancestry matrix for each position
@AncXAoA = ();
@AncAAoA = ();
for ($i = 0; $i < @OldFemaleXAoA; $i++){
  push @AncXAoA, [ @blank ];
  push @AncAAoA, [ @blank ];
}

for ($i = 0; $i < @OldFemaleXAoA; $i++){
  $anc = 1;
  $break = 0;
  for ($j = 0; $j < @PositionsX; $j++){
    while ( ($break < @{$OldFemaleXAoA[$i]} ) && ($OldFemaleXAoA[$i][$break] <= $PositionsX[$j])){
      $anc++;
      $break++;
    }
    if ( ($anc % 2) == 1){
      push @{$AncXAoA[$i]}, 1;
    }
    else{
      push @{$AncXAoA[$i]}, 0;
    }
  }
}
for ($i = 0; $i < @OldFemaleAAoA; $i++){
  $anc = 1;
  $break = 0;
  for ($j = 0; $j < @PositionsA; $j++){
    while ( ($break < @{$OldFemaleAAoA[$i]} ) && ($OldFemaleAAoA[$i][$break] <= $PositionsA[$j])){
      $anc++;
      $break++;
    }
    if ( ($anc % 2) == 1){
      push @{$AncAAoA[$i]}, 1;
    }
    else{
      push @{$AncAAoA[$i]}, 0;
    }
  }
}

#QTL simulation:  first get a random strength
#if ($r % $multiple == 0){  #use separate script for multi-cage analysis
@LocusStrengths = ();
$random = rand;
$LocusProp = $MinLocusProp + ($random * ($MaxLocusProp - $MinLocusProp));
$LocusStrength = $LocusProp / (1 - $LocusProp);
push @LocusStrengths, $LocusStrength;
push @PropVar, $LocusProp;

#simulate each x-linked locus... get a window position
for ($q = 0; $q < @AllQTLsX; $q++){
  @LocusPosX = ();
  push @LocusPosX, $AllQTLsX[$q];

#simulate phenotypes for each (diploid) female in last generation
  @phenotypes = ();
  for ($i = 0; $i < $NEachSex; $i++){
    push @phenotypes, 0;
  }
  for ($i = 0; $i < @LocusPosX; $i++){
#    for ($j = 0; $j < @PositionsX; $j++){
#      if ($PositionsX[$j] > $LocusPosX[$i]){
#	$pos = $j;
#	last;
#      }
    #    }
    $pos = $LocusPosX[$i];
    for ($j = 0; $j < @phenotypes; $j++){
      $allele = $j * 2;
      if ($AncXAoA[$allele][$pos] == 1){
	$phenotypes[$j] += $LocusStrengths[$i];
      }
      $allele++;
      if ($AncXAoA[$allele][$pos] == 1){
	$phenotypes[$j] += $LocusStrengths[$i];
      }
    }
  }

  for ($i = 0; $i < @phenotypes; $i++){
    $EnvVar = gaussian_rand();
    $EnvVar *= (0.5**0.5);  #makes it a similar effect on phenotypic variance as a codominant locus with an allelic effect of 1
    $phenotypes[$i] += $EnvVar;
  }
  
#phenotypic selection - make new ancestry matrices for individuals with highest and lowest phenotypic values
  @sorted = sort (@phenotypes);
  $LowerThresh = ($SelectionProp * @sorted) - 1;
  $LowerThresh = $sorted[$LowerThresh];
  $UpperThresh = @sorted - ($SelectionProp * @sorted);
  $UpperThresh = $sorted[$UpperThresh];
  
  @AncXHighAoA = ();
  @AncXLowAoA = ();
  for ($i = 0; $i < @phenotypes; $i++){
    if ($phenotypes[$i] <= $LowerThresh){
      $allele = $i * 2;
      @array = @{$AncXAoA[$allele]};
      push @AncXLowAoA, [ @array ];
      $allele++;
      @array = @{$AncXAoA[$allele]};
      push @AncXLowAoA, [ @array ];
    }
    elsif($phenotypes[$i] >= $UpperThresh){
      $allele = $i * 2;
      @array = @{$AncXAoA[$allele]};
      push @AncXHighAoA, [ @array ];
      $allele++;
      @array = @{$AncXAoA[$allele]};
      push @AncXHighAoA, [ @array ];
    }
  }

#simulate sequence read sampling at each site for high and low phenotype pools
  @AncPropsXLow = ();
  @AncPropsXHigh = ();
  for ($i = 0; $i < @PositionsX; $i++){
    $AncProp = 0;
    for ($j = 0; $j < $DepthsXLow[$i]; $j++){
      $allele = int(rand(@AncXLowAoA));
      if ($AncXLowAoA[$allele][$i] == 1){
	$AncProp += (1 / $DepthsXLow[$i]);
      }
    }
    push @AncPropsXLow, $AncProp;
    $AncProp = 0;
    for ($j = 0; $j < $DepthsXHigh[$i]; $j++){
      $allele = int(rand(@AncXHighAoA));
      if ($AncXHighAoA[$allele][$i] == 1){
	$AncProp += (1 / $DepthsXHigh[$i]);
      }
    }
    push @AncPropsXHigh, $AncProp;
    $AncDiff = $AncPropsXHigh[-1] - $AncPropsXLow[-1];
    $AncDiffsX[$i] += $AncDiff;
  }

#Smooth anc diffs
#  if ($r % $multiple == 0){
  @SmoothDiffs = ();
  for ($i = 0; $i < @AncDiffsX; $i++){
    $SmoothNum = $AncDiffsX[$i] * ($SmoothEachSide + 1);
    $SmoothDenom = $SmoothEachSide + 1;
    for ($j = 1; $j <= $SmoothEachSide; $j++){
      $weight = ( ($SmoothEachSide + 1) - $j);
      $pos = $i - $j;
      if ((($i - $j) >= 0) && ($WinChrsX[$i] eq $WinChrsX[$pos])){
	$SmoothDenom += $weight;
	$SmoothNum += ($AncDiffsX[$pos] * $weight);
      }
      $pos = $i + $j;
      if ((($i + $j) < @AncDiffsX) && ($WinChrsX[$i] eq $WinChrsX[$pos])){
	$SmoothDenom += $weight;
	$SmoothNum += ($AncDiffsX[$pos] * $weight);
      }
    }
    $SmoothNum = $SmoothNum / $SmoothDenom;
    push @SmoothDiffs, $SmoothNum
  }
  @AncDiffsX = @SmoothDiffs;
    
#Find "QTL peak" with highest ancestry difference between high and low phenotype pools
#Start at the real target and walk in both directions until hitting zero or the end of the chromosome
  $TargetWin = $LocusPosX[0];
  $TargetDiff = $AncDiffsX[$TargetWin];
  $MaxWin = $TargetWin;
  $MaxDiff = $TargetDiff;
  for ($i = $TargetWin - 1; $i >= 0; $i--){
    last if ($AncDiffsX[$i] < 0);
    if ($AncDiffsX[$i] > $MaxDiff){
      $MaxDiff = $AncDiffsX[$i];
      $MaxWin = $i;
    }
  }
  for ($i = $TargetWin + 1; $i < @AncDiffsX; $i++){
    last if ($AncDiffsX[$i] < 0);
    if ($AncDiffsX[$i] > $MaxDiff){
      $MaxDiff = $AncDiffsX[$i];
      $MaxWin = $i;
    }
  }
  push @PeakHeights, $MaxDiff;
#  $VertDev = $MaxDiff - $TargetDiff;
#  push @PrimDevs, $VertDev;
  push @PeakWindows, $MaxWin;
  
#Find the maximum secondary peak associated with this QTLs primary peak
  $MinSincePeak = $MaxDiff;
  $MaxDev = 0;
  for ($i = $MaxWin - 1; $i >= 0; $i--){
    last if ($AncDiffsX[$i] <= 0);
    last if ($CMStopsX[$i] > $CMStopsX[$i+1]);
    if ($AncDiffsX[$i] < $MinSincePeak){
      $MinSincePeak = $AncDiffsX[$i];
      next;
    }
    $VertDev = $AncDiffsX[$i] - $MinSincePeak;
    if ($VertDev > $MaxDev){
      $MaxDev = $VertDev;
    }
  }
  $MinSincePeak = $MaxDiff;
  for ($i = $MaxWin + 1; $i < @AncDiffsX; $i++){
    last if ($AncDiffsX[$i] <= 0);
    last if ($CMStopsX[$i] < $CMStopsX[$i-1]);
    if ($AncDiffsX[$i] < $MinSincePeak){
      $MinSincePeak = $AncDiffsX[$i];
      next;
    }
    $VertDev = $AncDiffsX[$i] - $MinSincePeak;
    if ($VertDev > $MaxDev){
      $MaxDev = $VertDev;
    }
  }
  push @SecDevs, $MaxDev;
  for ($i = 0; $i < @AncDiffsX; $i++){
    $AncDiffsX[$i] = 0;
  }
}


#then simulate an autosomal locus... get window position...
#if ($r % $multiple == 0){
for ($q = 0; $q < @AllQTLsA; $q++){
  @LocusPosA = ();
  push @LocusPosA, $AllQTLsA[$q];

#simulate phenotypes for each (diploid) female in last generation
  @phenotypes = ();
  for ($i = 0; $i < $NEachSex; $i++){
    push @phenotypes, 0;
  }
  for ($i = 0; $i < @LocusPosA; $i++){
#    for ($j = 0; $j < @PositionsA; $j++){
#      if ($PositionsA[$j] > $LocusPosA[$i]){
#	$pos = $j;
#	last;
#      }
#    }
    $pos = $LocusPosA[$i];
    for ($j = 0; $j < @phenotypes; $j++){
      $allele = $j * 2;
      if ($AncAAoA[$allele][$pos] == 1){
	$phenotypes[$j] += $LocusStrengths[$i];
      }
      $allele++;
      if ($AncAAoA[$allele][$pos] == 1){
	$phenotypes[$j] += $LocusStrengths[$i];
      }
    }
  }
  for ($i = 0; $i < @phenotypes; $i++){
    $EnvVar = gaussian_rand();
    $EnvVar *= (0.5**0.5);  #makes it a similar effect on phenotypic variance as a codominant locus with an allelic effect of 1
    $phenotypes[$i] += $EnvVar;
  }
  
#phenotypic selection - make new ancestry matrices for individuals with highest and lowest phenotypic values
  @sorted = sort (@phenotypes);
  $LowerThresh = ($SelectionProp * @sorted) - 1;
  $LowerThresh = $sorted[$LowerThresh];
  $UpperThresh = @sorted - ($SelectionProp * @sorted);
  $UpperThresh = $sorted[$UpperThresh];
  
  @AncAHighAoA = ();
  @AncALowAoA = ();
  for ($i = 0; $i < @phenotypes; $i++){
    if ($phenotypes[$i] <= $LowerThresh){
      $allele = $i * 2;
      @array = @{$AncAAoA[$allele]};
      push @AncALowAoA, [ @array ];
      $allele++;
      @array = @{$AncAAoA[$allele]};
      push @AncALowAoA, [ @array ];
    }
    elsif($phenotypes[$i] >= $UpperThresh){
      $allele = $i * 2;
      @array = @{$AncAAoA[$allele]};
      push @AncAHighAoA, [ @array ];
      $allele++;
      @array = @{$AncAAoA[$allele]};
      push @AncAHighAoA, [ @array ];
    }
  }
  
#simulate sequence read sampling at each site for high and low phenotype pools
  @AncPropsALow = ();
  @AncPropsAHigh = ();
  for ($i = 0; $i < @PositionsA; $i++){
    $AncProp = 0;
    for ($j = 0; $j < $DepthsALow[$i]; $j++){
      $allele = int(rand(@AncALowAoA));
      if ($AncALowAoA[$allele][$i] == 1){
	$AncProp += (1 / $DepthsALow[$i]);
      }
    }
    push @AncPropsALow, $AncProp;
    $AncProp = 0;
    for ($j = 0; $j < $DepthsAHigh[$i]; $j++){
      $allele = int(rand(@AncAHighAoA));
      if ($AncAHighAoA[$allele][$i] == 1){
	$AncProp += (1 / $DepthsAHigh[$i]);
      }
    }
    push @AncPropsAHigh, $AncProp;
    $AncDiff = $AncPropsAHigh[-1] - $AncPropsALow[-1];
    $AncDiffsA[$i] += $AncDiff;
  }
  
#Smooth anc diffs
#if ($r % $multiple == 0){
  @SmoothDiffs = ();
  for ($i = 0; $i < @AncDiffsA; $i++){
    $SmoothNum = $AncDiffsA[$i] * ($SmoothEachSide + 1);
    $SmoothDenom = $SmoothEachSide + 1;
    for ($j = 1; $j <= $SmoothEachSide; $j++){
      $weight = ( ($SmoothEachSide + 1) - $j);
      if (($i - $j) >= 0){
	$SmoothDenom += $weight;
	$pos = $i - $j;
	$SmoothNum += ($AncDiffsA[$pos] * $weight);
      }
      if (($i + $j) < @AncDiffsA){
	$SmoothDenom += $weight;
	$pos = $i + $j;
	$SmoothNum += ($AncDiffsA[$pos] * $weight);
      }
    }
    $SmoothNum = $SmoothNum / $SmoothDenom;
    push @SmoothDiffs, $SmoothNum
  }
  @AncDiffsA = @SmoothDiffs;
  
#Find "QTL peak" with highest ancestry difference between high and low phenotype pools
#Start at the real target and walk in both directions until hitting zero or the end of the chromosome
  $TargetWin = $LocusPosA[0];
  $TargetDiff = $AncDiffsA[$TargetWin];
  $MaxWin = $TargetWin;
  $MaxDiff = $TargetDiff;
  for ($i = $TargetWin - 1; $i >= 0; $i--){
    last if ($AncDiffsA[$i] < 0);
    last if ($CMStopsA[$i] > $CMStopsA[$i+1]);
    if ($AncDiffsA[$i] > $MaxDiff){
      $MaxDiff = $AncDiffsA[$i];
      $MaxWin = $i;
    }
  }
  for ($i = $TargetWin + 1; $i < @AncDiffsA; $i++){
    last if ($AncDiffsA[$i] < 0);
    last if ($CMStopsA[$i] < $CMStopsA[$i-1]);
    if ($AncDiffsA[$i] > $MaxDiff){
      $MaxDiff = $AncDiffsA[$i];
      $MaxWin = $i;
    }
  }
  push @PeakHeights, $MaxDiff;
#  $VertDev = $MaxDiff - $TargetDiff;
#  push @PrimDevs, $VertDev;

#Find the maximum secondary peak associated with this QTLs primary peak
  $MinSincePeak = $MaxDiff;
  $MaxDev = 0;
  for ($i = $MaxWin - 1; $i >= 0; $i--){
    last if ($AncDiffsA[$i] <= 0);
    last if ($CMStopsA[$i] > $CMStopsA[$i+1]);
    if ($AncDiffsA[$i] < $MinSincePeak){
      $MinSincePeak = $AncDiffsA[$i];
      next;
    }
    $VertDev = $AncDiffsA[$i] - $MinSincePeak;
    if ($VertDev > $MaxDev){
      $MaxDev = $VertDev;
    }
  }
  $MinSincePeak = $MaxDiff;
  for ($i = $MaxWin + 1; $i < @AncDiffsA; $i++){
    last if ($AncDiffsA[$i] <= 0);
    last if ($CMStopsA[$i] < $CMStopsA[$i-1]);
    if ($AncDiffsA[$i] < $MinSincePeak){
      $MinSincePeak = $AncDiffsA[$i];
      next;
    }
    $VertDev = $AncDiffsA[$i] - $MinSincePeak;
    if ($VertDev > $MaxDev){
      $MaxDev = $VertDev;
    }
  }
  push @SecDevs, $MaxDev;
  $MaxWin += @PositionsX;
  push @PeakWindows, $MaxWin;
  for ($i = 0; $i < @AncDiffsA; $i++){
    $AncDiffsA[$i] = 0;
  }
}
push @PeakHeightAoA, [ @PeakHeights ];
push @PeakWindowAoA, [ @PeakWindows ];
push @SecDevAoA, [ @SecDevs ];
@PeakHeights = ();
@PeakWindows = ();
@SecDevs = ();
print "Finished simulation replicate $r\n";
}

#output
open O, ">$OutputFile";
print O "PropVar\tHeight1\tPeakWin1\tSecDev1\tHeight2\tPeakWin2\tSecDev2\n";
for ($i = 0; $i < @PropVar; $i++){
  print O "$PropVar[$i]";
  for ($j = 0; $j < @{$PeakHeightAoA[$i]}; $j++){
    print O "\t$PeakHeightAoA[$i][$j]\t$PeakWindowAoA[$i][$j]\t$SecDevAoA[$i][$j]";
  }
  print O "\n";
}
close O;

