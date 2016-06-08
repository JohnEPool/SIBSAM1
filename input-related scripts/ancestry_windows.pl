#!/usr/bin/perl -w
#Ancestry estimation for introgression analysis
#use strict;

my @pop= ('AI_EF86N','AII_EF8','AIV_EF15','BVI_EF73','BVIII_CO13','CXII_CO9','CXIII_CO10','DIX_UG22','DX_UG11');
my @c= ('X','2L','2R','3L','3R');
my $i=0;
my $sites = 100;
for ($p=0; $p < @pop; $p++){

for ($i=0; $i < @c; $i++){
  my @Start=();
  my @Stop=();
  
  my $file_w = "windows_ZI_Chr$c[$i].txt";
  open(W, $file_w);
  while (<W>){
    chomp;
    my @W =split("\t",$_);
	push @Start, $W[0];
	push @Stop, $W[1];
	}
  close W;
  my $outfile = ">$pop[$p]/ZI_AI_Ancestry_$c[$i]_SNP.txt";
  open(O,$outfile);
	
  my $file = "$pop[$p]/AI_Ancestry_$c[$i]_"."$sites"."_SNP.txt";
  open(F, $file);

  my $l = 0;
  my $c = $c[$i];
  
  while (<F>){
    chomp;
    my $SumA = 0;
    my @F =split("\t",$_);
    for ($j=0;$j < @Start; $j++){
	  if ($F[1] >= $Start[$j] && $F[1] <= $Stop[$j]){
		$w= $j+1;
	  	print "$w\t$_\n";
	  	print O "$w\t$_\n";
		last;
		}
	  }
  }
  close F;
}
  close O;
}
exit;
