#!/usr/bin/perl -w
#Estimating parental ancestry proportion in mapping data (20% tails + 25% parental differences)
#use strict;

my @c= ('X','2L','2R','3L','3R');
my $i=0;
my $sites = 100;

for ($i=0; $i < @c; $i++){

  my $outfile = ">AI_Ancestry_$c[$i]_"."$sites"."_SNP.txt";
  open(O,$outfile);

  my $k=0;
  my @nuc= ('A','T','C','G');
  use List::Util qw( min max sum );

  my $mpileup = "ZI403N_UG11_DX-L1_DX-L2_DX-D1_DX-D2_$c[$i].sync";
  open(M, $mpileup);
  my $c=$c[$i];
  while (<M>){
    chomp;
    ($chr,$l, $ref, $PL, $PH, $FL1, $FL2, $FH1, $FH2)=split("\t",$_); #PL = low altitude parent, PH = high altitude parent, FL = low altitude phenotype in offspring, FH = high altitude phenotype in offspring
    my @PL = split(/:/,$PL);
    my @PH = split(/:/,$PH);
    my @FL1 = split(/:/,$FL1);
    my @FL2 = split(/:/,$FL2);
    my @FH1 = split(/:/,$FH1);
    my @FH2 = split(/:/,$FH2);
#For each site, skip if there are missing data 'N'
    if (major_allele(@PL) ==5 || major_allele(@PH) == 5 || major_allele(@FL1) == 5 || major_allele(@FH1) == 5 || major_allele(@FL2) == 5 || major_allele(@FH2) == 5){
      next;
    }
#Skip if both parentals only have reads for same allele
    if ((major_allele(@PL) eq major_allele(@PH)) && (major_allele_freq(@PL) == 1) && (major_allele_freq(@PH) == 1)) {
      next;
    }
    else {
#Estimate the frequency in PL of the major allele of PH "allele 1"
      $SumL = $PL[0] + $PL[1] + $PL[2] + $PL[3];
      $L = $PL[major_allele(@PH)]/$SumL;
      $Parent_diff = major_allele_freq(@PH) - $L;

#Skip site if the proportion of reads for allele 1 is <0.1
      if ($Parent_diff >= 0.25) {
		my @FL = ($FL1[0]+$FL2[0],$FL1[1]+$FL2[1],$FL1[2]+$FL2[2],$FL1[3]+$FL2[3]);
		my @FH = ($FH1[0]+$FH2[0],$FH1[1]+$FH2[1],$FH1[2]+$FH2[2],$FH1[3]+$FH2[3]);
#Freq of allele 1 in FL
		$SumFL = $FL[0] + $FL[1] + $FL[2] + $FL[3];
		$AL = $FL[major_allele(@PH)]/$SumFL;
#Freq of allele 1 in FH1
		$SumFH = $FH[0] + $FH[1] + $FH[2] + $FH[3];
		$AH = $FH[major_allele(@PH)]/$SumFH;
#Ancestry proportion of PH in FL1 
		$anc_PL_FL = ($AL - $L)/$Parent_diff;
#Ancestry proportion of PH in FH1
		$anc_PL_FH = ($AH - $L)/$Parent_diff;
#Ancestry proportion difference between FH and FL
		$ANC = $anc_PL_FH - $anc_PL_FL;
		print "$c\t$l\t$ANC\n";
		print O "$c\t$l\t$ANC\n";
		}
      else {
		next;
      	}
      }  
    }
  close O;
  close M;
}


exit;

sub major_allele {
$Major_allele = '';
my @nucp = @_[0 .. 3];
my $Sum = 0;
$Sum = $_[0] + $_[1] + $_[2] + $_[3];
if ($Sum == 0){
$Major_allele = 5;
}
else{
for ($a=0; $a<@nucp; $a++){
if ($nucp[$a] == max (@nucp)){
$Major_allele = $a;
}
}
}
return $Major_allele;
}

sub major_allele_freq {
$Major_allele_freq = 0;
my @nucp = @_[0 .. 3];
my $Sum = 0;
foreach (@nucp){
$Sum += $_;
}
if ($Sum == 0){
$Major_allele_freq = -1;
}
else{
for ($b=0; $b<@nucp; $b++){
if ($nucp[$b] == max (@nucp)){
$Major_allele_freq = $_[$b]/$Sum;
}
}
}
return $Major_allele_freq;
}
