#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my $InSAM = '';
my $HighParentAlleleIn = '';
my $LowParentAlleleIn = '';
my @Acodes = ('A','R','W','M');
my @Tcodes = ('T','Y','W','K');
my @Gcodes = ('G','R','S','K');
my @Ccodes = ('C','Y','S','M');
my @score = ('1','0.5','0.5','0.5');

my @chroms = ('2L','X','3L','2R','3R');
my @contigs = ('3','4','5','7','8');
my @pools = ('I','II','III','IV');
my $k = 0;
my $cmd = '';

for ($k = 0; $k < @pools; $k++) {

my $outfile = 'EF73NxZI418N_wing-' . $pools[$k] . '_ReadAncestry.txt'; #to fix...
open C, ">$outfile";
print C "Chromosome\tRead_1_Start\tRead_2_Start\tHighCounts\tLowCounts\tHomozygous_informative\tHeterozygous_informative\tComps\tHighProp\tLowProp\tWinner\n";

my $z = 0;

for ($z = 0; $z < @chroms; $z++){

	$InSAM = 'EF73NxZI418N-' . $pools[$k] . '-wingheader_Chr' . $chroms[$z] . '.sam'; #to fix...
	$HighParentAlleleIn = 'EF73N_Chr' . $chroms[$z] . '_diploid.fas'; #to fix...
	$LowParentAlleleIn = 'ZI418N_Chr' . $chroms[$z] . '_diploid.fas'; #to fix...
	
	$cmd = 'gunzip ' . $InSAM . '.gz';
	system($cmd);
	  
    my @line = ();
	my @highallele = ();
	my @lowallele = ();
	my $pos = '';
	my $a = 0;
	my $b = 0;
	
	open U, "<../$HighParentAlleleIn";
	while (<U>){
  		chomp;
  		last if m/^$/;
  		@line = split;
  		if ($line[0] ne '#'){
  			push @highallele, (split //, $line[0]);
		}
	}
	close U;
	print "Got parent 1\n";

	undef @line;

	open G, "<../$LowParentAlleleIn";
	while (<G>){
  		chomp;
  		last if m/^$/;
  		@line = split;
  		if ($line[0] ne '#'){
  			push @lowallele, (split //, $line[0]);
		}
	}
	close G;
	print "Got parent 2\n";

	undef @line;
	
	my $readname = '';
	my $readstart = '';
	my $pairstart = '';
	my @readseq = ();
	my @readtwoname = ();
	my @readtwostart = ();
	my @readtwo = ();
	my @line2 = ();
	open H, "<$InSAM";
	while (<H>){
  		chomp;
  		last if m/^$/;
  		@line2 = split;
  		next if ($line2[0] =~ m/@/);
#  		next if ($line2[7] > $line2[3]);
  		if ($line2[6] eq '=') {
  			push @readtwoname, $line2[0];
  			push @readtwostart, $line2[3];
  			push @readtwo, $line2[9];
  		}
  	}
  	close H;
  	$cmd = 'purge';
  	system($cmd);

	my $spot = 0;
	my $count=0;
	
	open D, "<$InSAM";
	while (<D>){
		my @readseq = ();
  		chomp;
  		last if m/^$/;
  		@line = split;
  		next if ($line[0] =~ m/@/);
  		next if (($line[7] < $line[3]) && ($line[7] != 0));
  		next if (($line[7] - $line[3]) > 10000);
  		$count++;
  		if ($count > 10000) {
		  	$cmd = 'purge';
  			system($cmd);
  			$count = 0;
  		}
  		if (($line[6] eq '=') || ($line[6] eq '*')) {
  			$readname = $line[0];
  			$readstart = $line[3];
  			$pairstart = $line[7];
  			push @readseq, (split //, $line[9]);
		
		my $comps = 0;
		my $high = 0;
		my $low = 0;
		my $hetinf = 0;
		my $hominf = 0;
		my $x = 0;
		$pos = ($readstart - 2);
		for ($x = 0; $x < @readseq; $x++){
			$pos++;
			if (($readseq[$x] ne 'N') && ($lowallele[$pos] ne 'N') && ($highallele[$pos] ne 'N')) {
				$comps++;
				if ($lowallele[$pos] ne $highallele[$pos]) {
					if (($lowallele[$pos] =~ m/A|T|C|G/) && ($highallele[$pos] =~ m/A|T|C|G/)){
						$hominf++;
						if ($readseq[$x] eq $lowallele[$pos]) {
							$low++;
						}
						elsif ($readseq[$x] eq $highallele[$pos]) {
							$high++;
						}
					}
					else {
					$hetinf++;
						if ($readseq[$x] eq 'A') {
							$a = 0;
							for ($a = 0; $a < @Acodes; $a++) {
								if ($lowallele[$pos] eq $Acodes[$a]) {
									$low += $score[$a];
									#last;
								}
							}
							$b=0;
							for ($b = 0; $b < @Acodes; $b++) {
								if ($highallele[$pos] eq $Acodes[$b]) {
									$high += $score[$b];
									#last;
								}
							}
						}
						if ($readseq[$x] eq 'T') {
							$a = 0;
							for ($a = 0; $a < @Tcodes; $a++) {
								if ($lowallele[$pos] eq $Tcodes[$a]) {
									$low += $score[$a];
									#last;
								}
							}
							$b=0;
							for ($b = 0; $b < @Tcodes; $b++) {
								if ($highallele[$pos] eq $Tcodes[$b]) {
									$high += $score[$b];
									#last;
								}
							}
						}
						if ($readseq[$x] eq 'C') {
							$a = 0;
							for ($a = 0; $a < @Ccodes; $a++) {
								if ($lowallele[$pos] eq $Ccodes[$a]) {
									$low += $score[$a];
									#last;
								}
							}
							$b=0;
							for ($b = 0; $b < @Ccodes; $b++) {
								if ($highallele[$pos] eq $Ccodes[$b]) {
									$high += $score[$b];
									#last;
								}
							}
						}
						if ($readseq[$x] eq 'G') {
							$a = 0;
							for ($a = 0; $a < @Gcodes; $a++) {
								if ($lowallele[$pos] eq $Gcodes[$a]) {
									$low += $score[$a];
									#last;
								}
							}
							$b=0;
							for ($b = 0; $b < @Gcodes; $b++) {
								if ($highallele[$pos] eq $Gcodes[$b]) {
									$high += $score[$b];
									#last;
								}
							}
						}
					}
				}
			}
		}

		if ($line[6] eq '=') {
			my $f = $spot;
			my $exception = 0;
			for ($f = $spot; $f < @readtwoname; $f++) {
				if ($readtwostart[$f] > $pairstart) {
					if ($exception == 0) {
						$f -= 100;
						$exception++;
					}
					else {
						$exception = 0;
						last;
					}
				}
  				if (($readtwoname[$f] eq $readname) && ($readtwostart[$f] == $pairstart)) {
  					my @readtwoseq = ();
					$pos = ($readtwostart[$f] - 2);
					push @readtwoseq, (split //, $readtwo[$f]);
					$x = 0;
					for ($x = 0; $x < @readtwoseq; $x++){
					$pos++;
					if (($readtwoseq[$x] ne 'N') && ($lowallele[$pos] ne 'N') && ($highallele[$pos] ne 'N')) {
						$comps++;
						if ($lowallele[$pos] ne $highallele[$pos]) {
							if (($lowallele[$pos] =~ m/A|T|C|G/) && ($highallele[$pos] =~ m/A|T|C|G/)){
							$hominf++;
								if ($readtwoseq[$x] eq $lowallele[$pos]) {
									$low++;
								}
							elsif ($readtwoseq[$x] eq $highallele[$pos]) {
								$high++;
							}
						}
						else {
						$hetinf++;
						if ($readseq[$x] eq 'A') {
							$a = 0;
							for ($a = 0; $a < @Acodes; $a++) {
								if ($lowallele[$pos] eq $Acodes[$a]) {
									$low += $score[$a];
									#last;
								}
							}
							$b=0;
							for ($b = 0; $b < @Acodes; $b++) {
								if ($highallele[$pos] eq $Acodes[$b]) {
									$high += $score[$b];
									#last;
								}
							}
						}
						if ($readseq[$x] eq 'T') {
							$a = 0;
							for ($a = 0; $a < @Tcodes; $a++) {
								if ($lowallele[$pos] eq $Tcodes[$a]) {
									$low += $score[$a];
									#last;
								}
							}
							$b=0;
							for ($b = 0; $b < @Tcodes; $b++) {
								if ($highallele[$pos] eq $Tcodes[$b]) {
									$high += $score[$b];
									#last;
								}
							}
						}
						if ($readseq[$x] eq 'C') {
							$a = 0;
							for ($a = 0; $a < @Ccodes; $a++) {
								if ($lowallele[$pos] eq $Ccodes[$a]) {
									$low += $score[$a];
									#last;
								}
							}
							$b=0;
							for ($b = 0; $b < @Ccodes; $b++) {
								if ($highallele[$pos] eq $Ccodes[$b]) {
									$high += $score[$b];
									#last;
								}
							}
						}
						if ($readseq[$x] eq 'G') {
							$a = 0;
							for ($a = 0; $a < @Gcodes; $a++) {
								if ($lowallele[$pos] eq $Gcodes[$a]) {
									$low += $score[$a];
									#last;
								}
							}
							$b=0;
							for ($b = 0; $b < @Gcodes; $b++) {
								if ($highallele[$pos] eq $Gcodes[$b]) {
									$high += $score[$b];
									#last;
								}
								}
							}
                    	}
					}
				}
			}
		$spot = $f;	
		last;
		}
    }
    	}

        my $highprop = '';
        my $lowprop = '';
        my $winner = '';
        if ($comps > 0) {
        $highprop = ($high/$comps);
        $lowprop = ($low/$comps);
        if ($highprop > $lowprop) {
        	$winner = 'H';
        }
        elsif ($highprop < $lowprop) {
        	$winner = 'L';	
        }
        else {
        	$winner = '0';
        }
        }
        else {
            $highprop = '0';
        	$lowprop = '0';
        	$winner = 'nocomps';
        }
        print C "$chroms[$z]\t$readstart\t$pairstart\t$high\t$low\t$hominf\t$hetinf\t$comps\t$highprop\t$lowprop\t$winner\n";
    }
}
    close D;
    close H;
    undef @line;
    undef @line2;
}
close C;
}