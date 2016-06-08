SIBSAM (Simulation-Based Inference for Bulk Segregant Analysis Mapping) - ReadMe

This document is a companion (not a substitute) for the paper describing this method.
As further described in that publication, SIBSAM is a QTL mapping method to facilitate
Bulk Segregant Analysis in Drosophila and species with comparable laboratory fecundity.

Applying SIBSAM involves a multi-stage process:
1.  Preparing SIBSAM input.
2.  Running null simulations with no QTLs.
3.  Running one locus simulations (genomes with one QTL).
4.  If needed, running cluster simulations with multiple linked QTLs.
5.  Generating and interpreting the final SIBSAM output.

PREPARING SIBSAM INPUT

The SIBSAM input file includes seven columns.  Following a header row, each subsequent
row gives information about one genomic window.  A template is given in the example 
folder called "EF8N_pigm_windows.txt".  SIBSAM expects this file to end in "_windows.txt".
The first column gives the chromosome arm that each window occurs on.  The second and 
third columns give the bp start and stop positions of that window.  Windows should be 
given in numeric order along each chromosome arm.

The fourth column gives the centiMorgan stop position of each window.  Thus, local
recombination rate estimates (or crude assumptions) are needed in order to run SIBSAM.
Note that SIBSAM is set up with Drosophila in mind, and does not simulate recombination
in males.  If this situation does not describe your organism, you can modify the sex-
averaged recombination rates you provide for each chromosome accordingly.  If you think
it would be important to have SIBSAM configured a different way for your application, you
are welcome to contact John Pool to ask about it (jpool@wisc.edu)

The fifth and sixth columns give information about allelic sampling depth in each of the
two bulk pools representing phenotypic extreme individuals from the mapping population.  
Here, the quantity of interest is the number of reads that are informative with regard to 
parental strain ancestry.  Dr. Justin Lack developed a script to collect this information 
from pool SAM files and parental strain FastA diploid consensus sequence files 
(informative_depth.pl).  We count a read pair once if it calls one or more sites that
have a fixed difference between parental strains, and we count it as one half if it calls
a site that is fixed in one parent but polymorphic in the other strain.  Simulations will
sample the same number of reads for each window as you report for the empirical data.
Each sampled read is assumed to be informative for ancestry difference calculations.

The seventh and final column gives the ancestry difference values calculated from the two
bulk pools - that is the average difference in the frequency of the "high" parental strain
allele between the high and low bulk pools.  One estimation procedure for ancestry 
difference is described in the application of SIBSAM to Drosophila pigmentation described 
by Bastide et al.  A pair of scripts applying this approach, written by Dr. Amir Yassin, 
are provided with SIBSAM.  ancestry_test_snp.pl evaluates an mpileup file at each site 
that differs between parental strains, while ancestry_windows.pl averages SNP ancestry
difference values for each of a set of defined windows (see "input-related files").  

RUNNING NULL SIMULATIONS

Perform null simulations by running sibsam_null.pl.  Like all SIBSAM scripts, it requires
the window input file and cmdline_args.pm to run.  Provide it with the following 
experimental and simulation parameters in the command line:
-g number of generations in the BSA experiment
-n number of individuals in each intermediate generation
-l number of individuals used for phenotypic selection in the final generation
-s proportion of individuals in each phenotypic extreme pool for sequencing
-c number of independent crosses to sum results across
-r number of simulated replicates
-f file name stem (the part before "_windows.txt" in your input file),
   which will also be the beginning of output file names
-b batch ID, optional flag, in case you are running these simulations in parallel

Analyze null simulations by running sibsam_null_analysis.pl.  This file must be located
in the same directory as one or more null simulation files named (filestem)_nullsims*
This script requires only the -f file stem command line argument.  It produces two output 
files, one labeled (filestem)_primpeaks_pvalues.txt.  For each primary peak that met 
minimal preliminary cutoffs, this file lists the window numbers, ancestry difference 
heights, and P values of these QTL peaks.  The other output file, (filestem)_secpeaks.txt, 
contains a preliminary list of secondary peaks that have not been tested for significance.

RUNNING SINGLE QTL SIMULATIONS

Single QTL simulations are run by sibsam_1locus.pl using exactly the same input flags as
for the null simulations.  Note that combining data from multiple crosses via -c has not
been implemented from this stage forward in SIBSAM, so this flag should be set to 1.  This
script requires primary peak P value file generated above.

Analysis of single QTL simulations is via sibsam_1locus_analysis.pl.  This script 
requires both output files generated by sibsam_null_analysis.pl. Its usage is the same as
for the null analysis script, with only the -f file stem command line argument. It 
produces a new primary peak output file and a new secondary peak output file.  The
primary peak output, (filestem)_primpeaks_fullresults.txt, adds output regarding QTL
strength and genomic confidence intervals for each previously-identified significant 
primary peak.  The secondary peak output file, (filestem)_secpeaks_pvalues.txt, provides
a P value for each secondary peak, as well as other information needed for cluster
simulations below. 

RUNNING MULTI-QTL SIMULATIONS

Simulations of multiple linked QTLs only need to be conducted if there is a significant
secondary peak associated with a given primary peak.  The first task is to generate a 
simulation input file for each peak cluster.  sibsam_cluster_prepare.pl does this in an 
automated way for all qualifying peak clusters at once.  It usage requires just the
-f file stem command line argument, and it will create input files name as 
(filestem)_cluster(X).txt (where X is the primary peak number associated with this cluster).
It also requires the two output files generated by the single QTL analysis script above.

Each cluster then needs to be simulated (separately) using sibsam_cluster.pl.  This script
requires the same command line arguments as the above two simulators, with the addition of
an -x flag denoting the cluster number to be simulated.  It returns updated estimates of
strength and location for each peak in the cluster (this time accounting for the presence 
of the other linked QTLs) in a file labeled (filestem)_cluster(X)results.txt.

INTEGRATING AND INTERPRETING SIBSAM RESULTS

The script sibsam_summarize.pl is designed to collect all relevant output files described
above, requiring only the -f file stem command line argument.  It produces the output file
(filestem)_finalresults.txt.  This is the only file users interested in the final results
of SIBSAM need to consult (all other files represent intermediate stages in the analysis 
pipeline).  It includes the location and confidence interval of each peak (in bp, rather
than window numbers), the P value and whether it's a primary or secondary peak, and the
estimate of effect size and its confidence interval.  The number of single QTL or cluster
simulations that were accepted during the rejection step is also indicated.

John Pool
jpool@wisc.edu
8 June 2016