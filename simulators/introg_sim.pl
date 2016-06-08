#!/usr/bin/perl -w
use strict;
use cmdline_args;
#use Math::Random;

if (!$ARGV[0] || $ARGV[0] eq "help"){
    print <<EOF;

This is a forward simulation program designed to test performance of introgression mapping in the Drosophila genome.

It returns the median cM between QTL peak and causative locus, as well as the proportion of cases where the distance between them is within some threshold.

It requires the following arguments:
-l [number of causative loci, evenly distributed across 5 chromosome arms, must be a multiple of 5]
-g [number of generations of crosses before sequencing (g>=2)]
-n [number of individuals (females plus males) in the generations before sequencing]
-f [number of females in the last generation used for sequencing]
-s [proportion of females in each phenotypic tail selected for sequencing (0 < s <= 0.5)]
-d [sequencing depth at each window]
-c [number of independent crosses to sum results across]
-r [number of replicates to simulate]

EOF
exit 1;
}

#Get command line arguments
my $NLoci = "";
my $gens = "";
my $F2N = "";
my $LastNf = "";
my $SelectionProp = "";
my $SeqDepth = "";
my $NCrosses = "";
my $TrueReps = "";

my %args;
$args{'-l'} = [ \$NLoci, 1];
$args{'-g'} = [ \$gens, 1];
$args{'-n'} = [ \$F2N, 1];
$args{'-f'} = [ \$LastNf, 1];
$args{'-s'} = [ \$SelectionProp, 1];
$args{'-d'} = [ \$SeqDepth, 1];
$args{'-c'} = [ \$NCrosses, 1];
$args{'-r'} = [ \$TrueReps, 1];

cmdline_args::get_options(\%args, \@ARGV);

#Define parameters not set by the command line
my $increment = 0.0002;  #how densely to sample parental line ancestry: 0.0002 is about 5kb on a 25Mb arm
my $cMThresh = 0.5;  #Calculate "power" based on how many QTL peaks are within this genetic distance of the true locus
my @ChromCMs = (75, 107, 110); 
my $EnvVarSD = 1;  #Environmental variance:  value of 1 is analogous to 20% of the total genetic effort

#Output files
my $OutputFile = 'IntrogSimFull_L' . $NLoci . '_F' . $gens . '_N' . $F2N . '_LN' . $LastN . '_S' . $SelectionProp . '_D' . $SeqDepth . '_C' . $NCrosses . '_R'  . $TrueReps . '.txt';
my $PowerFile = 'IntrogSimPower_L' . $NLoci . '_F' . $gens . '_N' . $F2N . '_LN' . $LastN . '_S' . $SelectionProp . '_D' . $SeqDepth . '_C' . $NCrosses . '_R'  . $TrueReps . '.txt';
my $AncDiffFile = 'IntrogSimAnc_L' . $NLoci . '_F' . $gens . '_N' . $F2N . '_LN' . $LastN . '_S' . $SelectionProp . '_D' . $SeqDepth . '_C' . $NCrosses . '_R'  . $TrueReps . '.txt';

#USER MODIFICATION NOT INTENDED BELOW THIS POINT

#Other definitions
my $F1N = $F2N;  #Not currently varying the number of F1's independently
my $LastN = $LastNf * 2;
my $replicates = $TrueReps * $NCrosses; #how many replicates to actually simulate
my $multiple = $NCrosses;  #combine signals from this many independent crosses

#Set up locus positions and effects (current code assumes all loci are of equal magnitude and are distributed evenly across genome
my @LocusPosX = ();
my @LocusPos2 = ();
my @LocusPos3 = ();
my @LocusStrengthX = ();
my @LocusStrength2 = ();
my @LocusStrength3 = ();

if (($NLoci % 5) > 0){
  die "Number of loci must be a multiple of 5\n";
}
my $LociPerArm = int( ($NLoci+0.5) / 5);
my $i = 0;
my $pos = (1 / $LociPerArm) / 2;
my $effect = 1 / $LociPerArm;
for ($i = 1; $i <= $LociPerArm; $i++){
  push @LocusPosX, $pos;
  push @LocusStrengthX, $effect;
  $pos += (1 / $LociPerArm);
}
$pos = (1 / $LociPerArm) / 4;
for ($i = 1; $i <= ($LociPerArm*2); $i++){
  push @LocusPos2, $pos;
  push @LocusPos3, $pos;
  push @LocusStrength2, $effect;
  push @LocusStrength3, $effect;
  $pos += (1 / $LociPerArm) / 2;
}

#old example of the results of the above code
#my @LocusPosX = (0.1,0.3,0.5,0.7,0.9);
#my @LocusPos2 = (0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95);
#my @LocusPos3 = (0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95);
#my @LocusStrengthX = (0.2,0.2,0.2,0.2,0.2);
#my @LocusStrength2 = (0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2);
#my @LocusStrength3 = (0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2);

#build breakpoint matrix for F1 individuals
my @OldFemaleXAoA = ();
my @OldFemale2AoA = ();
my @OldFemale3AoA = ();
my @OldMaleXAoA = ();
my @OldMale2AoA = ();
my @OldMale3AoA = ();
my @NewFemaleXAoA = ();
my @NewFemale2AoA = ();
my @NewFemale3AoA = ();
my @NewMaleXAoA = ();
my @NewMale2AoA = ();
my @NewMale3AoA = ();

my @IndBreaks1 = (); #chromosomes from from 1st parental line
my @IndBreaks2 = (0); #chromosomes from from 2nd parental line

#Factorial subroutine
my $a = 0;
my $b = 0;
my $FactIn = 0;
sub factorial{
    $a = 1;
    $b = 1;
    while ($b <= $FactIn) {
	$a *= $b;
	$b++;
    }
    return $a;
}

#Probability of a given number of recombination events (Poisson) for each chromosome
#(allowing up to 10 recombination events per chromosome per generation)
#Cumulative probabilities are given for each number, for easy use with random numbers later.
my $TotalLength = $ChromCMs[0];
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
my @RecombProbs2 = ();
$TotalLength = $ChromCMs[1];
$ExpRecEvents = 0.01 * $TotalLength;
$PrSum = 0;
for ($i = 0; $i < 11; $i++){
  $FactIn = $i;
  $Pr = ((exp(-1 * $ExpRecEvents)) * ($ExpRecEvents ** $i)) / &factorial;
  $PrSum = $PrSum + $Pr;
  push @RecombProbs2, $PrSum;
}
my @RecombProbs3 = ();
$TotalLength = $ChromCMs[2];
$ExpRecEvents = 0.01 * $TotalLength;
$PrSum = 0;
for ($i = 0; $i < 11; $i++){
  $FactIn = $i;
  $Pr = ((exp(-1 * $ExpRecEvents)) * ($ExpRecEvents ** $i)) / &factorial;
  $PrSum = $PrSum + $Pr;
  push @RecombProbs3, $PrSum;
}

#set up arrays
my @positions = ();
my @AncDiffsX = ();
my @AncDiffs2 = ();
my @AncDiffs3 = ();
for ($i = 0; $i <= 1; $i += $increment){
  push @positions, $i;
  push @AncDiffsX, 0;
  push @AncDiffs2, 0;
  push @AncDiffs3, 0;
}

#Using pre-cooked random normal values with mean 0 and SD 1, to avoid requiring math::random on cluster
my @NormalDist = (8.790e-1, -1.280e+0,  4.100e-1,  7.030e-2, -6.960e-1, -2.600e-1,  7.610e-1,  2.580e-1, -1.040e+0,  4.190e-1,  3.140e-1,  1.450e+0, -7.090e-1, -8.910e-1,  3.920e-1, -6.530e-1,  4.280e-1, -4.330e-1,  5.540e-1,  8.730e-1,  1.880e-1, -3.220e-2,  5.030e-1, -7.110e-1,  6.900e-1,  6.480e-1, -1.740e+0, -5.530e-1, -6.820e-1, -2.170e-1,  2.970e-1,  1.510e-1,  5.390e-1,  1.150e+0, -2.860e-1,  5.120e-1, -2.640e-1,  1.430e-1,  1.440e-1,  1.630e+0, -1.300e-1,  1.380e+0,  1.080e+0, -6.560e-1,  9.700e-1, -6.860e-1, -1.360e-1, -1.480e-1,  1.340e+0,  2.120e-1, -7.840e-1, -9.730e-1, -3.020e-1, -1.900e+0,  2.060e+0, -1.170e-1, -1.630e-2,  6.470e-1,  1.140e+0,  8.710e-1, -7.770e-1,  2.420e-1,  1.030e+0, -1.630e+0,  9.240e-1, -6.800e-1,  9.570e-1,  3.580e-1,  8.450e-1,  9.470e-1,  3.600e-1,  9.810e-1,  1.370e+0, -8.080e-1,  1.790e+0,  7.760e-1,  1.870e+0,  2.360e-1, -5.460e-1, -2.780e-1, -8.380e-1,  1.500e+0, -1.280e-1,  9.830e-2, -1.420e+0,  8.910e-3, -5.750e-1,  3.980e-2,  1.630e+0, -5.770e-1,  1.650e+0, -5.930e-1, -1.060e+0, -1.480e+0,  9.680e-2,  1.760e-1, -1.870e+0, -2.120e-1, -4.640e-1,  4.390e-2,  3.350e-1, -8.030e-1,  5.780e-1,  1.290e-1,  1.160e+0, -7.920e-1, -6.150e-1,  6.800e-1, -4.910e-1, -1.730e-1,  3.170e-1, -1.550e+0,  1.300e+0, -3.870e-1, -6.240e-1,  5.430e-1,  4.540e-1,  9.980e-1, -2.550e+0, -1.340e+0, -9.210e-1,  4.050e-1,  7.190e-2,  1.040e+0,  1.840e+0, -1.710e+0,  9.870e-1,  1.490e+0,  6.620e-1, -2.820e-1,  9.310e-1, -4.530e-1,  3.860e-2,  8.080e-2,  6.150e-1,  3.480e-1, -2.140e+0, -1.010e+0,  1.420e+0, -5.200e-1,  5.210e-1, -1.710e+0, -1.710e+0,  8.710e-1, -1.280e+0,  6.210e-1,  1.690e+0,  1.180e+0, -7.890e-1, -3.410e-1, -1.210e+0,  5.340e-1, -2.910e-1,  9.150e-1, -6.810e-1, -3.260e-1,  7.860e-1, -1.520e+0,  1.040e+0, -3.080e-1, -2.780e-1, -3.400e-1, -1.230e+0,  1.550e-2, -1.190e+0, -4.890e-1, -7.260e-1,  7.450e-1,  3.670e-2, -5.310e-1,  2.530e+0,  8.220e-1, -3.780e-1, -9.200e-1, -8.170e-1, -4.830e-1, -3.640e-1, -1.540e+0,  2.230e-1,  8.060e-1,  2.000e+0,  9.580e-1,  1.720e-1, -4.870e-1,  1.180e-1, -2.390e-1,  8.110e-1,  1.620e+0,  2.350e-1,  6.270e-1, -3.700e-1,  1.440e+0, -4.340e-2, -1.090e-2, -3.940e-1,  1.330e+0, -7.850e-1,  2.940e-1, -1.000e+0, -8.200e-1,  4.040e-1,  1.100e+0,  2.250e-1, -2.230e+0, -4.670e-1,  1.750e-1,  2.500e+0,  8.110e-1,  1.600e-1,  1.810e+0, -3.850e-1,  7.110e-1, -6.400e-1,  4.840e-1, -1.460e+0, -7.000e-1, -6.990e-1, -4.390e-3, -2.620e-1,  1.820e-2,  4.020e-1,  7.200e-1, -7.400e-3, -5.880e-1, -4.100e-1,  1.160e+0,  2.710e+0, -4.640e-1,  1.760e-1,  1.540e-1,  1.960e+0, -2.850e-1, -7.040e-1, -4.160e-1,  4.410e-1, -8.650e-1,  8.820e-1,  9.050e-1,  7.650e-1, -2.130e+0, -9.470e-1, -6.790e-1,  2.640e-1,  3.630e-1,  1.010e-1, -2.290e+0, -1.590e+0,  1.900e+0, -7.510e-1,  5.310e-1, -1.090e+0, -1.940e-1,  4.860e-2, -1.550e+0, -1.030e+0,  1.540e+0, -1.610e+0,  1.140e+0, -1.570e+0,  1.240e+0, -6.850e-1,  1.480e+0, -1.270e+0, -5.570e-1,  1.380e+0,  1.660e-1, -1.290e-1,  2.050e+0, -1.410e+0, -5.310e-2, -9.560e-1, -1.820e-1,  2.150e-1, -1.620e+0, -1.680e+0, -1.170e+0,  2.410e+0,  9.110e-2, -5.360e-1,  1.520e+0,  1.200e+0, -8.940e-1,  3.680e-1, -1.210e+0,  7.240e-1, -4.920e-1, -8.390e-3,  6.680e-1, -2.560e-1,  1.930e+0, -3.030e-1, -1.490e+0, -1.240e-1,  1.370e-1,  1.160e+0, -2.450e-1, -3.180e-1, -3.570e-1, -2.130e-1, -9.540e-1, -3.170e-1,  6.970e-1, -1.450e+0, -2.120e-1, -7.700e-1,  1.270e+0,  8.900e-1,  2.710e-1, -1.320e+0, -6.920e-1, -4.970e-1,  5.230e-1,  4.680e-1, -3.730e-1,  1.330e-1, -6.740e-1,  1.350e+0, -1.340e+0, -1.700e-1, -1.680e+0, -1.730e+0, -3.260e-2,  1.220e+0, -3.830e-2, -4.110e-1, -6.120e-1,  6.980e-1, -2.870e-1, -2.970e-1,  4.550e-1,  8.730e-1,  1.560e+0,  1.410e+0, -1.200e+0,  5.610e-2, -1.780e+0,  6.960e-2, -1.530e-1, -1.360e+0,  5.090e-1, -2.260e-1, -4.240e-1,  5.270e-1, -7.460e-1,  2.220e-1,  1.020e+0,  1.550e+0, -9.180e-1, -3.410e-2,  2.570e-1, -8.480e-1, -7.010e-1, -2.350e-1, -2.450e-1, -1.650e+0,  3.900e-1, -1.730e+0, -1.310e-1, -1.140e+0,  1.660e-1, -2.720e-1, -1.420e+0, -5.490e-1, -8.190e-1, -8.750e-1,  1.020e+0, -4.680e-1, -4.740e-1,  9.610e-1, -7.360e-1, -9.950e-1, -6.580e-1, -1.560e-1, -8.800e-1, -1.820e+0, -2.420e-1, -1.990e-1, -3.320e-1,  1.630e+0,  7.070e-1, -4.790e-1, -7.540e-1, -2.080e+0,  2.230e-1,  6.830e-1, -1.660e-1, -7.260e-1,  1.200e-1,  2.570e-1, -7.740e-1,  1.740e+0, -8.320e-1,  5.540e-1, -7.050e-1, -1.150e+0,  8.560e-1,  1.190e+0,  2.010e+0, -1.470e+0, -9.340e-1,  1.700e+0, -2.360e-1, -1.610e+0, -9.100e-1, -6.250e-2, -1.240e+0, -2.090e+0,  2.690e+0,  2.300e-1, -1.290e+0,  2.890e-1,  2.100e-1, -1.880e-1, -1.650e+0,  1.090e+0, -1.060e+0,  2.130e+0,  2.960e-1,  5.490e-1, -9.700e-1,  1.950e-1,  1.060e+0, -1.120e+0, -3.410e-1, -1.920e-1, -7.300e-1, -1.070e+0, -1.200e+0,  8.440e-1,  4.570e-1,  1.660e-1, -8.140e-1, -1.230e+0,  2.260e-1, -2.370e-1, -9.250e-1,  8.680e-1, -8.200e-2,  6.180e-1, -7.770e-1, -1.660e-1, -2.400e+0, -4.330e-2,  1.380e+0,  1.600e+0,  1.490e+0, -1.520e+0, -4.520e-1, -6.750e-1,  1.210e+0,  1.890e-1, -5.230e-2, -1.040e+0,  7.580e-1, -6.450e-3, -1.030e+0,  5.170e-1,  4.850e-1, -3.280e-1,  9.490e-1,  5.720e-1,  1.090e-1, -6.890e-1, -8.770e-1, -1.220e+0, -8.860e-1, -3.890e-1, -7.400e-1,  1.260e-1, -1.540e+0,  1.130e+0,  2.060e+0, -1.180e+0, -5.470e-1,  8.600e-1, -4.660e-1,  7.140e-1, -1.450e+0,  7.570e-1,  2.250e-1,  1.150e+0,  9.990e-1, -6.420e-1, -6.420e-1,  1.030e+0, -1.140e+0,  1.950e-1, -1.090e+0,  9.270e-1, -1.050e-1,  1.300e+0, -4.780e-1, -6.690e-1, -7.310e-1,  5.160e-1, -1.130e+0,  1.130e+0, -1.140e+0, -9.430e-1,  1.300e+0,  7.250e-1,  2.230e+0, -1.550e+0, -2.250e-1,  9.700e-1,  3.630e-1, -6.410e-1,  3.020e-1, -1.850e+0,  2.340e-1,  5.680e-1, -7.170e-2,  1.310e+0,  4.180e-1, -1.870e+0, -1.550e+0, -6.410e-1,  6.060e-2, -1.170e+0,  1.130e+0,  9.150e-2, -1.320e+0,  1.180e+0,  7.430e-3,  6.810e-1, -7.530e-1, -1.210e+0,  1.090e+0,  1.340e-1,  1.330e+0,  1.780e-1,  6.140e-1, -3.240e-1, -2.930e-1,  6.220e-1, -1.720e-1, -1.460e+0,  3.350e-1, -1.500e+0,  2.970e-1,  1.840e+0,  1.160e-2, -3.240e-1, -1.820e+0, -8.690e-1, -3.130e-1,  2.020e+0, -5.090e-2,  2.400e+0,  9.490e-1, -7.140e-1, -4.220e-1,  3.980e-1,  2.560e-1,  1.230e-1, -1.950e+0, -4.150e-1, -6.640e-1,  3.210e-2,  1.020e+0,  8.790e-1,  7.370e-1, -5.630e-1, -2.870e-1,  1.170e+0, -9.890e-1, -1.450e+0,  1.740e-1, -7.970e-1,  9.510e-1,  4.250e-2,  1.750e-1, -6.500e-1, -8.710e-1,  6.190e-1, -1.050e+0,  7.590e-1,  5.740e-1, -1.120e+0,  5.280e-1,  1.550e+0, -1.610e+0, -1.750e+0,  8.940e-1, -1.540e+0,  1.470e+0, -6.420e-2, -3.610e-1,  1.800e-1, -6.060e-1, -1.460e+0, -7.440e-1, -1.700e+0, -6.240e-1,  1.120e+0, -4.910e-2, -7.890e-1,  5.220e-1,  5.160e-1, -2.060e+0,  7.350e-1, -1.790e+0,  5.960e-1,  7.290e-1, -1.000e-2, -7.970e-1,  2.630e-1, -1.940e-1, -1.510e+0, -5.240e-1, -7.900e-1, -4.280e-1,  4.880e-1, -2.070e-1, -3.790e-1, -3.140e-1, -2.100e-1, -3.040e-1,  1.950e-1,  6.010e-1, -5.470e-1,  2.370e-1,  1.030e+0,  1.280e+0,  1.920e-1, -1.320e+0,  7.510e-2,  6.620e-2,  4.460e-2,  1.180e+0, -5.150e-2,  8.190e-1, -7.760e-1,  1.430e+0, -1.850e+0,  1.290e+0,  9.950e-1, -6.230e-1,  8.080e-1,  2.570e-1,  3.630e-1,  5.560e-1, -1.970e+0,  1.610e-1, -2.210e+0, -1.190e+0,  2.300e+0,  1.810e+0, -2.660e-1,  4.380e-1, -1.780e+0,  4.450e-1,  1.740e-1, -1.840e+0,  2.210e+0, -3.920e-1,  3.640e-1,  6.470e-2, -7.110e-1, -4.200e-1, -6.590e-1, -2.640e+0, -7.180e-1, -2.120e+0, -1.700e+0, -4.940e-1, -1.310e+0, -1.950e+0,  6.710e-1, -1.340e+0, -2.570e-1, -4.900e-1,  4.970e-1, -1.350e+0,  1.960e+0,  8.060e-1, -1.300e+0, -4.960e-1,  1.770e-1,  5.220e-1, -1.860e+0, -5.520e-1, -6.570e-1, -3.350e-1,  8.870e-2, -1.690e-1,  9.420e-1,  1.690e+0, -4.330e-1,  9.500e-1,  1.240e+0,  1.510e+0,  3.820e-1, -7.740e-1,  2.370e-2, -1.260e+0, -8.370e-1, -1.160e-1, -1.300e+0,  8.030e-1, -1.360e+0,  7.680e-2,  9.120e-1,  1.920e-1,  7.460e-1, -1.100e+0,  5.020e-1, -3.140e-1,  2.480e+0, -5.380e-1, -2.190e-1, -3.180e-1,  8.180e-2,  2.140e-1, -1.380e+0, -5.950e-1, -5.390e-1, -6.190e-1,  8.100e-1,  1.970e+0, -7.470e-1,  4.680e-1,  2.810e-1, -1.950e-2,  4.940e-1, -8.360e-1,  6.010e-2, -2.020e-1, -1.210e+0, -8.310e-1,  3.850e-1,  1.070e+0, -9.960e-1, -1.830e-1,  1.130e-1,  1.250e-1, -8.710e-2,  3.680e-1,  6.710e-1,  1.690e+0, -4.900e-1,  6.080e-1, -1.780e+0, -8.100e-1, -5.120e-1, -1.040e+0,  7.050e-1, -6.290e-1,  4.580e-2,  4.540e-1, -1.050e+0, -1.410e+0, -2.530e-1,  8.720e-1, -8.020e-1, -1.900e-1,  1.420e+0, -3.450e-1, -2.290e-1,  1.030e-1, -1.270e+0, -8.520e-1, -6.730e-1,  1.280e-1,  1.460e+0, -1.010e+0, -7.510e-1,  1.680e-1,  1.090e+0, -4.780e-1, -6.140e-2, -3.110e-1, -6.910e-1, -7.360e-1,  3.050e-1, -5.710e-1,  5.200e-1, -2.030e+0, -1.750e+0,  5.810e-1, -1.880e+0, -4.460e-1,  9.970e-2, -9.300e-1,  1.620e-1,  4.690e-1, -1.140e+0,  1.040e+0, -6.490e-1,  4.910e-1, -1.960e-1, -3.670e-1, -5.520e-1,  5.450e-1,  1.330e+0, -5.560e-1,  4.740e-1,  1.370e+0,  9.860e-1, -2.170e+0, -5.070e-1, -7.660e-1,  1.030e+0, -1.220e+0, -8.120e-3,  1.290e+0,  6.240e-1,  1.390e+0,  8.980e-1, -5.160e-1, -1.080e+0, -8.200e-1,  1.890e-1,  1.660e+0, -2.520e-2,  1.340e+0, -1.470e-1, -1.710e-1,  3.990e-1,  6.840e-1, -2.450e-1, -1.000e+0,  1.490e+0,  5.140e-1, -2.710e-1, -8.910e-1,  2.120e-1, -3.580e-1, -2.930e-1,  2.520e-1, -2.320e+0, -4.370e-2, -1.560e-1,  1.300e+0,  1.980e-1, -1.570e+0, -8.920e-2, -7.300e-3, -1.850e+0,  6.000e-1,  3.200e-1,  9.260e-1,  1.160e+0,  5.210e-1,  1.450e+0,  2.200e-1,  1.090e-1,  9.050e-1, -6.720e-1, -5.630e-1, -1.010e+0,  4.440e-2, -4.300e-1, -2.030e-1, -6.870e-1,  1.460e+0, -2.480e+0, -1.300e+0,  6.230e-1,  1.630e+0, -5.700e-1,  4.420e-1, -7.000e-2, -2.410e-1,  3.120e-1,  1.890e+0,  5.420e-1, -2.190e-1, -1.770e-1,  7.480e-1,  3.330e-1, -1.580e+0, -6.000e-1,  1.890e+0,  1.560e+0, -1.180e-1,  2.320e-1, -1.840e+0, -2.070e+0,  4.090e-1,  7.770e-1,  7.580e-1, -6.390e-1, -1.160e+0,  7.280e-2, -1.130e-1, -5.830e-1,  9.720e-1, -1.950e-2,  2.810e-1,  7.040e-1,  4.910e-1, -1.020e+0, -2.280e-1,  2.750e-1, -1.820e-1, -7.990e-1, -1.620e+0, -2.100e-1, -6.060e-1,  6.570e-1, -9.670e-1, -8.120e-2,  1.720e-1,  8.210e-1,  7.710e-1,  1.270e+0, -6.490e-1,  1.650e+0,  6.310e-1, -1.210e+0, -1.240e+0,  1.260e+0, -1.670e-1, -7.930e-1, -3.490e-2, -4.280e-1,  8.400e-1,  4.500e-2,  1.020e+0,  3.670e-1, -7.310e-1,  1.190e+0,  2.980e-1, -3.900e-2,  3.430e-1,  2.260e+0, -6.860e-1,  1.980e-1, -4.270e-1,  9.000e-2,  8.920e-1,  2.830e-1, -8.910e-1,  7.430e-1,  4.590e-1, -3.530e-1, -1.490e+0, -2.560e+0, -5.720e-1, -2.900e-1,  8.730e-1,  8.890e-1, -2.220e-1,  4.430e-1,  4.480e-1, -1.180e+0, -8.200e-2,  1.760e+0,  5.730e-1,  2.880e-1,  4.240e-1, -4.370e-1, -4.770e-1, -7.150e-1,  2.250e+0,  1.150e-2,  4.590e-1,  2.460e-1, -4.470e-1,  2.960e-1, -1.900e-1, -2.160e-1,  2.890e-1, -1.520e+0, -1.250e-1,  1.090e+0,  4.170e-1,  1.010e+0,  1.590e+0,  2.030e-1,  3.290e-1, -3.400e-1, -1.650e-1,  1.860e-1, -2.270e+0,  3.990e-1, -6.120e-1,  2.310e+0, -2.470e-1,  3.540e-1, -1.400e-1,  9.840e-1, -6.330e-1,  8.360e-1,  1.760e-1, -5.780e-1, -3.670e-1,  1.750e+0,  1.760e-1, -3.860e-1,  2.150e-1,  1.460e+0, -1.100e+0,  5.110e-1, -1.310e+0,  6.940e-1, -1.630e+0,  4.980e-1, -1.060e+0);
for ($i = 0; $i < @NormalDist; $i++){
  $NormalDist[$i] *= $EnvVarSD;
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
my $MaxPos = 0;
my $RegionStart = 0;
my $RegionStop = 0;
my $PowerCount = 0;
my @RecombPos = ();
my @AncBreaks = ();
my @A1Breaks = ();
my @A2Breaks = ();
my @ancestries = ();
my @BreakpointNumbers = ();
my @AncPropsX = ();
my @AncProps2 = ();
my @AncProps3 = ();
my @AncestryTracker = ();
my @BreakpointTracker = ();
my @AncXAoA = ();
my @Anc2AoA = ();
my @Anc3AoA = ();
my @blank = ();
my @phenotypes = ();
my @sorted = ();
my @array = ();
my @AncXHighAoA = ();
my @Anc2HighAoA = ();
my @Anc3HighAoA = ();
my @AncXLowAoA = ();
my @Anc2LowAoA = ();
my @Anc3LowAoA = ();
my @AncPropsXHigh = ();
my @AncPropsXLow = ();
my @AncProps2High = ();
my @AncProps2Low = ();
my @AncProps3High = ();
my @AncProps3Low = ();
my @PeakDist = ();
my @PeakDistAoA = ();
my @AllPeakDists = ();
my @AncDiffAoA = ();
my @SelectedFemales = ();

#Replicate loop (not indented)
for ($r = 1; $r <= $replicates; $r++){

@OldFemaleXAoA = ();
@OldFemale2AoA = ();
@OldFemale3AoA = ();
@OldMaleXAoA = ();
@OldMale2AoA = ();
@OldMale3AoA = ();
@NewFemaleXAoA = ();
@NewFemale2AoA = ();
@NewFemale3AoA = ();
@NewMaleXAoA = ();
@NewMale2AoA = ();
@NewMale3AoA = ();

for ($i = 0; $i < ($F1N / 2); $i++){
  push @OldFemaleXAoA, [ @IndBreaks1 ];
  push @OldFemaleXAoA, [ @IndBreaks2 ];
  push @OldFemale2AoA, [ @IndBreaks1 ];
  push @OldFemale2AoA, [ @IndBreaks2 ];
  push @OldFemale3AoA, [ @IndBreaks1 ];
  push @OldFemale3AoA, [ @IndBreaks2 ];
  push @OldMale2AoA, [ @IndBreaks1 ];
  push @OldMale2AoA, [ @IndBreaks2 ];
  push @OldMale3AoA, [ @IndBreaks1 ];
  push @OldMale3AoA, [ @IndBreaks2 ];
}
for ($i = 0; $i < ($F1N / 4); $i++){
  push @OldMaleXAoA, [ @IndBreaks1 ];
  push @OldMaleXAoA, [ @IndBreaks2 ];
}
$NEachSex = $F2N / 2;

#For each individual in the next generation, choose parents, enact recombination from mother, transfer breakpoints
for ($g = 2; $g <= $gens; $g++){

#if it's an even-numbered generation, choose parents from previous generation, then perform phenotypic selection
 if (($g % 2) == 0){
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

#for chromosome 2, determine the number of recombination events (store that many random breakpoints) from maternal alleles
    $random = rand;
    @RecombPos = ();
    for ($j = 0; $j < @RecombProbs2; $j++){
      if ($random > $RecombProbs2[$j]){
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
      @AncBreaks = @{$OldFemale2AoA[$PAllele1]};
      push @NewFemale2AoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldFemale2AoA[$PAllele1]};
      @A2Breaks = @{$OldFemale2AoA[$PAllele2]};
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
      push @NewFemale2AoA, [ @AncBreaks ];
    }

#for chromosome 3, determine the number of recombination events (store that many random breakpoints) from maternal alleles
    $random = rand;
    @RecombPos = ();
    for ($j = 0; $j < @RecombProbs3; $j++){
      if ($random > $RecombProbs3[$j]){
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
      @AncBreaks = @{$OldFemale3AoA[$PAllele1]};
      push @NewFemale3AoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldFemale3AoA[$PAllele1]};
      @A2Breaks = @{$OldFemale3AoA[$PAllele2]};
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
      push @NewFemale3AoA, [ @AncBreaks ];
    }

#for females in the next generation, choose father
    $random = int(rand($parents));
    $random2 = rand;
    if ($random2 < 0.5){
      $PAllele1 = $random * 2;
    }
    else{
      $PAllele1 = $random * 2 + 1; 
    }
    $PAllele2 = $random;  #used for X (fewer male X's)

#transfer paternal chromosomes without recombination
    @AncBreaks = @{$OldMaleXAoA[$PAllele2]};
    push @NewFemaleXAoA, [ @AncBreaks ];
    @AncBreaks = @{$OldMale2AoA[$PAllele1]};
    push @NewFemale2AoA, [ @AncBreaks ];
    @AncBreaks = @{$OldMale3AoA[$PAllele1]};
    push @NewFemale3AoA, [ @AncBreaks ];
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

#for chromosome 2, determine the number of recombination events (store that many random breakpoints) from maternal alleles
    $random = rand;
    @RecombPos = ();
    for ($j = 0; $j < @RecombProbs2; $j++){
      if ($random > $RecombProbs2[$j]){
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
      @AncBreaks = @{$OldFemale2AoA[$PAllele1]};
      push @NewMale2AoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldFemale2AoA[$PAllele1]};
      @A2Breaks = @{$OldFemale2AoA[$PAllele2]};
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
      push @NewMale2AoA, [ @AncBreaks ];
    }

#for chromosome 3, determine the number of recombination events (store that many random breakpoints) from maternal alleles
    $random = rand;
    @RecombPos = ();
    for ($j = 0; $j < @RecombProbs3; $j++){
      if ($random > $RecombProbs3[$j]){
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
      @AncBreaks = @{$OldFemale3AoA[$PAllele1]};
      push @NewMale3AoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldFemale3AoA[$PAllele1]};
      @A2Breaks = @{$OldFemale3AoA[$PAllele2]};
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
      push @NewMale3AoA, [ @AncBreaks ];
    }

#for males in the next generation, choose father
    $random = int(rand($parents));
    $random2 = rand;
    if ($random2 < 0.5){
      $PAllele1 = $random * 2;
    }
    else{
      $PAllele1 = $random * 2 + 1; 
    }

#transfer paternal chromosomes without recombination
    @AncBreaks = @{$OldMale2AoA[$PAllele1]};
    push @NewMale2AoA, [ @AncBreaks ];
    @AncBreaks = @{$OldMale3AoA[$PAllele1]};
    push @NewMale3AoA, [ @AncBreaks ];
  }

#perform phenotypic selection (done outside this loop if it's the last generation)
  if ($g < $gens){
#individual ancestry matrix for each position
    @AncXAoA = ();
    @Anc2AoA = ();
    @Anc3AoA = ();
    for ($i = 0; $i < @NewFemaleXAoA; $i++){
      push @AncXAoA, [ @blank ];
      push @Anc2AoA, [ @blank ];
      push @Anc3AoA, [ @blank ];
    }
    
    for ($i = 0; $i < @NewFemaleXAoA; $i++){
      $anc = 1;
      $break = 0;
      for ($j = 0; $j < @positions; $j++){
	while ( ($break < @{$NewFemaleXAoA[$i]} ) && ($NewFemaleXAoA[$i][$break] <= $positions[$j])){
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
    for ($i = 0; $i < @NewFemale2AoA; $i++){
      $anc = 1;
      $break = 0;
      for ($j = 0; $j < @positions; $j++){
	while ( ($break < @{$NewFemale2AoA[$i]} ) && ($NewFemale2AoA[$i][$break] <= $positions[$j])){
	  $anc++;
	  $break++;
	}
	if ( ($anc % 2) == 1){
	  push @{$Anc2AoA[$i]}, 1;
	}
	else{
	  push @{$Anc2AoA[$i]}, 0;
	}
      }
    }
    for ($i = 0; $i < @NewFemale3AoA; $i++){
      $anc = 1;
      $break = 0;
      for ($j = 0; $j < @positions; $j++){
	while ( ($break < @{$NewFemale3AoA[$i]} ) && ($NewFemale3AoA[$i][$break] <= $positions[$j])){
	  $anc++;
	  $break++;
	}
	if ( ($anc % 2) == 1){
	  push @{$Anc3AoA[$i]}, 1;
	}
	else{
	  push @{$Anc3AoA[$i]}, 0;
	}
      }
    }
    
#simulate phenotypes for each (diploid) female in this generation
    @phenotypes = ();
    for ($i = 0; $i < $NEachSex; $i++){
      push @phenotypes, 0;
    }
    for ($i = 0; $i < @LocusPosX; $i++){
      $pos = ($LocusPosX[$i] / $increment);
      for ($j = 0; $j < @phenotypes; $j++){
	$allele = $j * 2;
	if ($AncXAoA[$allele][$pos] == 1){
	  $phenotypes[$j] += $LocusStrengthX[$i];
	}
	$allele++;
	if ($AncXAoA[$allele][$pos] == 1){
	  $phenotypes[$j] += $LocusStrengthX[$i];
	}
      }
    }
    for ($i = 0; $i < @LocusPos2; $i++){
      $pos = ($LocusPos2[$i] / $increment);
      for ($j = 0; $j < @phenotypes; $j++){
	$allele = $j * 2;
	if ($Anc2AoA[$allele][$pos] == 1){
	  $phenotypes[$j] += $LocusStrength2[$i];
	}
	$allele++;
	if ($Anc2AoA[$allele][$pos] == 1){
	  $phenotypes[$j] += $LocusStrength2[$i];
	}
      }
    }
    for ($i = 0; $i < @LocusPos3; $i++){
      $pos = ($LocusPos3[$i] / $increment);
      for ($j = 0; $j < @phenotypes; $j++){
	$allele = $j * 2;
	if ($Anc3AoA[$allele][$pos] == 1){
	  $phenotypes[$j] += $LocusStrength3[$i];
	}
	$allele++;
	if ($Anc3AoA[$allele][$pos] == 1){
	  $phenotypes[$j] += $LocusStrength3[$i];
	}
      }
    }
    for ($i = 0; $i < @phenotypes; $i++){
#  $EnvVar = Math::Random::random_normal(1,0,$EnvVarSD);  #Using pre-generated random numbers to avoid module requirement
      $random = rand(@NormalDist);
      $EnvVar = $NormalDist[$random];
      $phenotypes[$i] += $EnvVar;
    }
    
#phenotypic selection - select females with highest phenotypic values, then indicate their selected alleles
    @sorted = sort (@phenotypes);
    $LowerThresh = ($SelectionProp * @sorted) - 1;
    $LowerThresh = $sorted[$LowerThresh];
    $UpperThresh = @sorted - ($SelectionProp * @sorted);
    $UpperThresh = $sorted[$UpperThresh];
    
    @SelectedFemales = ();
    for ($i = 0; $i < @phenotypes; $i++){
      if($phenotypes[$i] >= $UpperThresh){
	push @SelectedFemales, 1;
      }
      else{ 
	push @SelectedFemales, 0;
      }
    }
#transfer new ancestry breakpoint matrices to old matrices to prepare for next generation
  @OldFemaleXAoA = @NewFemaleXAoA;
  @NewFemaleXAoA = ();
  @OldFemale2AoA = @NewFemale2AoA;
  @NewFemale2AoA = ();
  @OldFemale3AoA = @NewFemale3AoA;
  @NewFemale3AoA = ();
  @OldMaleXAoA = @NewMaleXAoA;
  @NewMaleXAoA = ();
  @OldMale2AoA = @NewMale2AoA;
  @NewMale2AoA = ();
  @OldMale3AoA = @NewMale3AoA;
  @NewMale3AoA = ();  
  }
 }

###if we're in an odd number generation, choose female parents from those phenotypically selected, and males from parental line 2
 else{  
#for females in the next generation, choose mother (must be phenotypically selected)
  for ($i = 0; $i < $NEachSex; $i++){
    $random = int(rand($parents));
    while ($SelectedFemales[$random] == 0){
      $random = int(rand($parents));
    }
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

#for chromosome 2, determine the number of recombination events (store that many random breakpoints) from maternal alleles
    $random = rand;
    @RecombPos = ();
    for ($j = 0; $j < @RecombProbs2; $j++){
      if ($random > $RecombProbs2[$j]){
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
      @AncBreaks = @{$OldFemale2AoA[$PAllele1]};
      push @NewFemale2AoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldFemale2AoA[$PAllele1]};
      @A2Breaks = @{$OldFemale2AoA[$PAllele2]};
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
      push @NewFemale2AoA, [ @AncBreaks ];
    }

#for chromosome 3, determine the number of recombination events (store that many random breakpoints) from maternal alleles
    $random = rand;
    @RecombPos = ();
    for ($j = 0; $j < @RecombProbs3; $j++){
      if ($random > $RecombProbs3[$j]){
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
      @AncBreaks = @{$OldFemale3AoA[$PAllele1]};
      push @NewFemale3AoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldFemale3AoA[$PAllele1]};
      @A2Breaks = @{$OldFemale3AoA[$PAllele2]};
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
      push @NewFemale3AoA, [ @AncBreaks ];
    }

#for females in the next generation, transfer pure paternal chromosomes from parental line 2
    push @NewFemaleXAoA, [ @IndBreaks2 ];
    push @NewFemale2AoA, [ @IndBreaks2 ];
    push @NewFemale3AoA, [ @IndBreaks2 ];
  }

#for males in the next generation, choose mother
  for ($i = 0; $i < $NEachSex; $i++){
    $random = int(rand($parents));
    while ($SelectedFemales[$random] == 0){
      $random = int(rand($parents));
    }
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

#for chromosome 2, determine the number of recombination events (store that many random breakpoints) from maternal alleles
    $random = rand;
    @RecombPos = ();
    for ($j = 0; $j < @RecombProbs2; $j++){
      if ($random > $RecombProbs2[$j]){
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
      @AncBreaks = @{$OldFemale2AoA[$PAllele1]};
      push @NewMale2AoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldFemale2AoA[$PAllele1]};
      @A2Breaks = @{$OldFemale2AoA[$PAllele2]};
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
      push @NewMale2AoA, [ @AncBreaks ];
    }

#for chromosome 3, determine the number of recombination events (store that many random breakpoints) from maternal alleles
    $random = rand;
    @RecombPos = ();
    for ($j = 0; $j < @RecombProbs3; $j++){
      if ($random > $RecombProbs3[$j]){
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
      @AncBreaks = @{$OldFemale3AoA[$PAllele1]};
      push @NewMale3AoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldFemale3AoA[$PAllele1]};
      @A2Breaks = @{$OldFemale3AoA[$PAllele2]};
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
      push @NewMale3AoA, [ @AncBreaks ];
    }

#for males in the next generation, transfer pure paternal chromosomes from parental line 2
    push @NewMale2AoA, [ @IndBreaks2 ];
    push @NewMale3AoA, [ @IndBreaks2 ];
  }

#transfer new ancestry breakpoint matrices to old matrices to prepare for next generation
  @OldFemaleXAoA = @NewFemaleXAoA;
  @NewFemaleXAoA = ();
  @OldFemale2AoA = @NewFemale2AoA;
  @NewFemale2AoA = ();
  @OldFemale3AoA = @NewFemale3AoA;
  @NewFemale3AoA = ();
  @OldMaleXAoA = @NewMaleXAoA;
  @NewMaleXAoA = ();
  @OldMale2AoA = @NewMale2AoA;
  @NewMale2AoA = ();
  @OldMale3AoA = @NewMale3AoA;
  @NewMale3AoA = ();
 }
}

###
#for ($i = 0; $i < @OldFemale2AoA; $i++){
#  for ($j = 0; $j < @{$OldFemale2AoA[$i]}; $j++){
#    print "$OldFemale2AoA[$i][$j]\t";
#  }
#  print "\n";
#}
###

#individual ancestry matrix for each position
@AncXAoA = ();
@Anc2AoA = ();
@Anc3AoA = ();
for ($i = 0; $i < @NewFemaleXAoA; $i++){
  push @AncXAoA, [ @blank ];
  push @Anc2AoA, [ @blank ];
  push @Anc3AoA, [ @blank ];
}

for ($i = 0; $i < @NewFemaleXAoA; $i++){
  $anc = 1;
  $break = 0;
  for ($j = 0; $j < @positions; $j++){
    while ( ($break < @{$NewFemaleXAoA[$i]} ) && ($NewFemaleXAoA[$i][$break] <= $positions[$j])){
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
for ($i = 0; $i < @NewFemale2AoA; $i++){
  $anc = 1;
  $break = 0;
  for ($j = 0; $j < @positions; $j++){
    while ( ($break < @{$NewFemale2AoA[$i]} ) && ($NewFemale2AoA[$i][$break] <= $positions[$j])){
      $anc++;
      $break++;
    }
    if ( ($anc % 2) == 1){
      push @{$Anc2AoA[$i]}, 1;
    }
    else{
      push @{$Anc2AoA[$i]}, 0;
    }
  }
}
for ($i = 0; $i < @NewFemale3AoA; $i++){
  $anc = 1;
  $break = 0;
  for ($j = 0; $j < @positions; $j++){
    while ( ($break < @{$NewFemale3AoA[$i]} ) && ($NewFemale3AoA[$i][$break] <= $positions[$j])){
      $anc++;
      $break++;
    }
    if ( ($anc % 2) == 1){
      push @{$Anc3AoA[$i]}, 1;
    }
    else{
      push @{$Anc3AoA[$i]}, 0;
    }
  }
}

#simulate phenotypes for each (diploid) female in last generation
@phenotypes = ();
for ($i = 0; $i < $NEachSex; $i++){
  push @phenotypes, 0;
}
for ($i = 0; $i < @LocusPosX; $i++){
  $pos = ($LocusPosX[$i] / $increment);
  for ($j = 0; $j < @phenotypes; $j++){
    $allele = $j * 2;
    if ($AncXAoA[$allele][$pos] == 1){
      $phenotypes[$j] += $LocusStrengthX[$i];
    }
    $allele++;
    if ($AncXAoA[$allele][$pos] == 1){
      $phenotypes[$j] += $LocusStrengthX[$i];
    }
  }
}
for ($i = 0; $i < @LocusPos2; $i++){
  $pos = ($LocusPos2[$i] / $increment);
  for ($j = 0; $j < @phenotypes; $j++){
    $allele = $j * 2;
    if ($Anc2AoA[$allele][$pos] == 1){
      $phenotypes[$j] += $LocusStrength2[$i];
    }
    $allele++;
    if ($Anc2AoA[$allele][$pos] == 1){
      $phenotypes[$j] += $LocusStrength2[$i];
    }
  }
}
for ($i = 0; $i < @LocusPos3; $i++){
  $pos = ($LocusPos3[$i] / $increment);
  for ($j = 0; $j < @phenotypes; $j++){
    $allele = $j * 2;
    if ($Anc3AoA[$allele][$pos] == 1){
      $phenotypes[$j] += $LocusStrength3[$i];
    }
    $allele++;
    if ($Anc3AoA[$allele][$pos] == 1){
      $phenotypes[$j] += $LocusStrength3[$i];
    }
  }
}
for ($i = 0; $i < @phenotypes; $i++){
#  $EnvVar = Math::Random::random_normal(1,0,$EnvVarSD);
#Using pre-generated random numbers to avoid module requirement
  $random = rand(@NormalDist);
  $EnvVar = $NormalDist[$random];
  $phenotypes[$i] += $EnvVar;
}

#phenotypic selection - make new ancestry matrices for individuals with highest and lowest phenotypic values
@sorted = sort (@phenotypes);
$LowerThresh = ($SelectionProp * @sorted) - 1;
$LowerThresh = $sorted[$LowerThresh];
$UpperThresh = @sorted - ($SelectionProp * @sorted);
$UpperThresh = $sorted[$UpperThresh];

@AncXHighAoA = ();
@Anc2HighAoA = ();
@Anc3HighAoA = ();
@AncXLowAoA = ();
@Anc2LowAoA = ();
@Anc3LowAoA = ();
for ($i = 0; $i < @phenotypes; $i++){
  if ($phenotypes[$i] <= $LowerThresh){
    $allele = $i * 2;
    @array = @{$AncXAoA[$allele]};
    push @AncXLowAoA, [ @array ];
    @array = @{$Anc2AoA[$allele]};
    push @Anc2LowAoA, [ @array ];
    @array = @{$Anc3AoA[$allele]};
    push @Anc3LowAoA, [ @array ];
    $allele++;
    @array = @{$AncXAoA[$allele]};
    push @AncXLowAoA, [ @array ];
    @array = @{$Anc2AoA[$allele]};
    push @Anc2LowAoA, [ @array ];
    @array = @{$Anc3AoA[$allele]};
    push @Anc3LowAoA, [ @array ];
  }
  elsif($phenotypes[$i] >= $UpperThresh){
    $allele = $i * 2;
    @array = @{$AncXAoA[$allele]};
    push @AncXHighAoA, [ @array ];
    @array = @{$Anc2AoA[$allele]};
    push @Anc2HighAoA, [ @array ];
    @array = @{$Anc3AoA[$allele]};
    push @Anc3HighAoA, [ @array ];
    $allele++;
    @array = @{$AncXAoA[$allele]};
    push @AncXHighAoA, [ @array ];
    @array = @{$Anc2AoA[$allele]};
    push @Anc2HighAoA, [ @array ];
    @array = @{$Anc3AoA[$allele]};
    push @Anc3HighAoA, [ @array ];
  }
}

#simulate sequence read sampling at each site for high and low phenotype pools
@AncPropsXLow = ();
@AncPropsXHigh = ();
@AncProps2Low = ();
@AncProps2High = ();
@AncProps3Low = ();
@AncProps3High = ();
#@AncDiffsX = ();
#@AncDiffs2 = ();
#@AncDiffs3 = ();
for ($i = 0; $i < @positions; $i++){
  $AncProp = 0;
  for ($j = 0; $j < $SeqDepth; $j++){
    $allele = int(rand(@AncXLowAoA));
    if ($AncXLowAoA[$allele][$i] == 1){
      $AncProp += (1 / $SeqDepth);
    }
  }
  push @AncPropsXLow, $AncProp;
  $AncProp = 0;
  for ($j = 0; $j < $SeqDepth; $j++){
    $allele = int(rand(@Anc2LowAoA));
    if ($Anc2LowAoA[$allele][$i] == 1){
      $AncProp += (1 / $SeqDepth);
    }
  }
  push @AncProps2Low, $AncProp;
  $AncProp = 0;
  for ($j = 0; $j < $SeqDepth; $j++){
    $allele = int(rand(@Anc3LowAoA));
    if ($Anc3LowAoA[$allele][$i] == 1){
      $AncProp += (1 / $SeqDepth);
    }
  }
  push @AncProps3Low, $AncProp;
  $AncProp = 0;
  for ($j = 0; $j < $SeqDepth; $j++){
    $allele = int(rand(@AncXHighAoA));
    if ($AncXHighAoA[$allele][$i] == 1){
      $AncProp += (1 / $SeqDepth);
    }
  }
  push @AncPropsXHigh, $AncProp;
  $AncProp = 0;
  for ($j = 0; $j < $SeqDepth; $j++){
    $allele = int(rand(@Anc2HighAoA));
    if ($Anc2HighAoA[$allele][$i] == 1){
      $AncProp += (1 / $SeqDepth);
    }
  }
  push @AncProps2High, $AncProp;
  $AncProp = 0;
  for ($j = 0; $j < $SeqDepth; $j++){
    $allele = int(rand(@Anc3HighAoA));
    if ($Anc3HighAoA[$allele][$i] == 1){
      $AncProp += (1 / $SeqDepth);
    }
  }
  push @AncProps3High, $AncProp;

  $AncDiff = $AncPropsXHigh[-1] - $AncPropsXLow[-1];
  $AncDiffsX[$i] += $AncDiff;
  $AncDiff = $AncProps2High[-1] - $AncProps2Low[-1];
  $AncDiffs2[$i] += $AncDiff;
  $AncDiff = $AncProps3High[-1] - $AncProps3Low[-1];
  $AncDiffs3[$i] += $AncDiff;
}

#For chromosome 2, add AncDiffs2 to a new output AoA (to visualize individual cross outcomes)
if ($r <= 100){
  push @AncDiffAoA, [ @AncDiffs2 ];
}

#For each chromosome arm, find "QTL peak" with highest ancestry difference between high and low phenotype pools, calculate distance to real locus
if ($r % $multiple == 0){

@PeakDist = ();
for ($k = 0; $k < @LocusPosX; $k++){
  if ($k == 0){
    $RegionStart = 0;
  }
  else{
    $RegionStart = ($LocusPosX[$k] + $LocusPosX[$k-1]) / 2;
    $RegionStart = int(($RegionStart + ($increment / 2)) * @positions);
  }
  if ($k == (@LocusPosX - 1)){
    $RegionStop = @positions - 1;
  }
  else{
    $RegionStop = ($LocusPosX[$k] + $LocusPosX[$k+1]) / 2;
    $RegionStop = int(($RegionStop + ($increment / 2)) * @positions);
  }

  $MaxDiff = 0;
  for ($i = $RegionStart; $i <= $RegionStop; $i++){
    if ($AncDiffsX[$i] > $MaxDiff){
      $MaxDiff = $AncDiffsX[$i];
      $MaxPos = $positions[$i];
    }
  }
  $j = abs($MaxPos - $LocusPosX[$k]) * $ChromCMs[0];
  push @PeakDist, $j;
  push @AllPeakDists, $j;
  if ($j < $cMThresh){
    $PowerCount++;
  }
}

for ($k = 0; $k < @LocusPos2; $k++){
  if ($k == 0){
    $RegionStart = 0;
  }
  else{
    $RegionStart = ($LocusPos2[$k] + $LocusPos2[$k-1]) / 2;
    $RegionStart = int(($RegionStart + ($increment / 2)) * @positions);
  }
  if ($k == (@LocusPos2 - 1)){
    $RegionStop = @positions - 1;
  }
  else{
    $RegionStop = ($LocusPos2[$k] + $LocusPos2[$k+1]) / 2;
    $RegionStop = int(($RegionStop + ($increment / 2)) * @positions);
  }

  $MaxDiff = 0;
  for ($i = $RegionStart; $i <= $RegionStop; $i++){
    if ($AncDiffs2[$i] > $MaxDiff){
      $MaxDiff = $AncDiffs2[$i];
      $MaxPos = $positions[$i];
    }
  }
  $j = abs($MaxPos - $LocusPos2[$k]) * $ChromCMs[1];
  push @PeakDist, $j;
  push @AllPeakDists, $j;
  if ($j < $cMThresh){
    $PowerCount++;
  }
}

for ($k = 0; $k < @LocusPos3; $k++){
  if ($k == 0){
    $RegionStart = 0;
  }
  else{
    $RegionStart = ($LocusPos3[$k] + $LocusPos3[$k-1]) / 2;
    $RegionStart = int(($RegionStart + ($increment / 2)) * @positions);
  }
  if ($k == (@LocusPos3 - 1)){
    $RegionStop = @positions - 1;
  }
  else{
    $RegionStop = ($LocusPos3[$k] + $LocusPos3[$k+1]) / 2;
    $RegionStop = int(($RegionStop + ($increment / 2)) * @positions);
  }

  $MaxDiff = 0;
  for ($i = $RegionStart; $i <= $RegionStop; $i++){
    if ($AncDiffs3[$i] > $MaxDiff){
      $MaxDiff = $AncDiffs3[$i];
      $MaxPos = $positions[$i];
    }
  }
  $j = abs($MaxPos - $LocusPos3[$k]) * $ChromCMs[2];
  push @PeakDist, $j;
  push @AllPeakDists, $j;
  if ($j < $cMThresh){
    $PowerCount++;
  }
}

push @PeakDistAoA, [ @PeakDist ];

print "Replicate $r distances:  ";
for ($i = 0; $i < @PeakDist; $i++){
  printf "%7.5f",$PeakDist[$i];
  print "\t";
}
print "\n";

for ($i = 0; $i < @AncDiffsX; $i++){
  $AncDiffsX[$i] = 0;
  $AncDiffs2[$i] = 0;
  $AncDiffs3[$i] = 0;
}
}
}

#output
open O, ">$OutputFile";
for ($i = 0; $i < @PeakDistAoA; $i++){
  for ($j = 0; $j < @{$PeakDistAoA[$i]}; $j++){
    print O $PeakDistAoA[$i][$j];
    if ($j < (@{$PeakDistAoA[$i]} - 1)){
      print O "\t";
    }
    else{
      print O "\n";
    }
  }
}
close O;

@AllPeakDists = sort(@AllPeakDists);
my $MedianDist = 0;
if ((@AllPeakDists % 2) == 0){
  $k = @AllPeakDists / 2;
  $MedianDist = ($AllPeakDists[$k] + $AllPeakDists[$k-1]) / 2;
}
else{
  $k = (@AllPeakDists-1) / 2;
  $MedianDist = $AllPeakDists[$k];
}
my $power = $PowerCount / @AllPeakDists;

open P, ">$PowerFile";
print P "Loci\tGens\tN\tLastN\tSelProp\tSeqDepth\tNCrosses\tMedianCMdist\tPowerWithin0.5cM\n";
print P "$NLoci\t$gens\t$F2N\t$LastN\t$SelectionProp\t$SeqDepth\t$NCrosses\t$MedianDist\t$power";
close P;

open A, ">$AncDiffFile";
for ($i = 0; $i < @AncDiffAoA; $i++){
  for ($j = 0; $j < @{$AncDiffAoA[$i]}; $j++){
    print A "$AncDiffAoA[$i][$j]\t";
  }
  print A "\n";
}
close A;

