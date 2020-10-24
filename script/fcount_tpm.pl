#!/usr/bin/env perl 

## add the TPM column to the featureCounts output.
## follow the simple historical approach (feature length = eff. length) 

use strict; 
use warnings; 

if (scalar @ARGV != 1) {
  print STDERR "Usage: ./fcount_tpm.pl <featureCounts_out>\n";
  exit 1
}

my $fcout = shift @ARGV; 

open FCOUT,"<",$fcout or die "$!"; 

## FC output has two header lines 
my $h1 = <FCOUT>; 
my $h2 = <FCOUT>; 
printf "$h1$h2";

my $scaling = 0; 
my @lines; 
my @rpk; 

while (<FCOUT>) { 
  chomp; 
  push @lines,$_; 
  my @t = split /\t+/; 
  my $rpk = $t[6]/$t[5]; 
  $scaling += $rpk; 
  push @rpk,$rpk; 
} 

$scaling /= 1000000;

for (my $i=0; $i<scalar @lines; $i++) {
  printf "%s\t%.3f\n",$lines[$i],$rpk[$i]/$scaling; 
} 

close FCOUT; 

  
