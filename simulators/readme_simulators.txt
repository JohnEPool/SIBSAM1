Simulators ReadMe

bsa_sim.pl
introg_sim.pl

Note - you do not need to run these scripts in order to use SIBSAM.

These perl scripts simulate Bulk Segregant Analysis and Introgression
Mapping experiments, respectively.  They are intended to facilitate
the optimization of experimental design. They simulate a cross between 
two parental strains, and in each following generation, they monitor
parental strain ancestry along the genomes of each individual in the
mapping population.  XX females and XY males are simulated, with two
pairs of autosomes (motivated by Drosophila melanogaster chromosomes 2
and 3).  Random mating and female recombination occur each generation. 
Phenotypes for selection (enabled in the last generation for BSA, and
every 2nd generation for IM) are generated based on genetic and
environmental/random contributions, as detailed below.  Genotyping
effort (e.g. the number of ancestry-informative sequencing reads per
marker or window) is also simulated.

The scripts require the file cmdline_args.pm in the same directory in
order to operate. The experimental variables that can be defined through
the command line include:

-l The number of causative loci of equal effect to simulate.  These will
be evenly distributed across chromosome arms X, 2L, 2R, 3L, and 3R. 
These program versions thus require a number of causative loci divisible
by 5.  Contact John Pool if you are interested in a different
configuration.

-g The number of generations in the mapping experiment before
sequencing.  For example, if g=12, 12th generation offspring are
sequenced.  For introgression mapping, this number should be a multiple
of 2.

-n The number of individuals in the mapping population each generation
prior to the last generation.

-f The number of individuals used for phenotyping in the last
generation.  The user can allow more individuals to be generated in the
last generation, and might assume that only one sex will be phenotyped.

-s The proportion of individuals in each phenotypic tail selected for
sequencing. (0 < s <= 0.5)

-d The allelic sampling depth at each marker/window.  This could be the
number of sequence reads per window that are informative with regard to
parental ancestry.

-c The number of independent mapping crosses to sum ancestry difference
results across. Here, it is assumed that all crosses have the same QTLs.

-r The number of simulation replicates to perform.

Additional user-definable parameters exist within the perl script (just
below the above command line options):

$increment tells the simulator how densely to place markers/windows
along the chromosome. The chromosome is represented as a 0 to 1
interval, with markers every $increment

$cMThresh is your target for how close an ancestry difference peak
should be to the real QTL location, in centiMorgans.  The simulator will
report how often this condition is met.

@ChromCMs specifies the centiMorgan length of the X and the 2 autosomes.

$EnvVarSD controls the strength of random variance due to environmental
factors or measurement error.  It is scaled such that a value of 1
implies 20% of the total genetic effect.

The simulators return 3 output files:

A "power" file giving the median distance between statistic peak and QTL
across all simulated QTLs from all replicates, as well as the proportion
of such distances that fell within cMTresh.

A "full" results file giving the above distance for each QTL in each 
simulated replicate.

An "ancestry difference" file giving sample output of QTL surfaces.
This is given only for the first 100 replicates for chromosome 2 (the
first autosome).  It is meant to visualize QTL mapping outcomes.

John Pool
jpool@wisc.edu

8 June 2016





