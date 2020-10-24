#!/usr/bin/env perl 

use strict; 
use warnings; 
## use Data::Dumper; 

## use translated sequences of putative reference CDS mapped onto genomes 
## and identify improper CDS; change feature type to `misc` - which applies e.g. to pseudogenes

if (scalar @ARGV != 2) {
  print STDERR "Usage: ./modify_reference_fasta.pl <protein_fa> <ref_fa>\n";
  exit 1
}

my $prot_fa = shift @ARGV; 
my $ref_fa = shift @ARGV; 

open PROT_FA,"<",$prot_fa; 
open REF_FA,"<",$ref_fa; 

my $seq_name;
my $protein = {};
my $names = {}; 
 
while (<PROT_FA>) {
  if (m/^>(.*?)\n/) {  
    $seq_name = $1;
  } else {
    chomp; 
    $protein->{$seq_name} = (defined $protein->{$seq_name}) ? join "",$protein->{$seq_name},$_ : $_;
  } 
}

## print Dumper $protein; 

foreach $seq_name (keys %{$protein}) { 
  $seq_name =~ m/^(.*?)_/; 
  my $cds_name = $1;
  ## we screen for 3 possible types of problems: 
  ## 1) start AA is not IVLM; 2) there's no stop codon; 3) there's a stop (X) in the middle of the seq 
  if ($protein->{$seq_name} !~ m/X$/) { 
    printf STDERR "Sequence %s (%s) does not end with a stop codon; reference sequence %s type is changed to misc.\n",$seq_name,$protein->{$seq_name},$cds_name; 
    $names->{$cds_name} = "misc"; 
  } elsif ($protein->{$seq_name} !~ m/^[ILVM]/) { 
    printf STDERR "Sequence %s (%s) does not start with I/V/L/M; reference sequence %s type is changed to misc.\n",$seq_name,$protein->{$seq_name},$cds_name; 
    $names->{$cds_name} = "misc"; 
  } elsif ($protein->{$seq_name} =~ m/X[A-Z]/) { 
    printf STDERR "Sequence %s (%s) has a stop codon in the middle; reference sequence %s type is changed to misc.\n",$seq_name,$protein->{$seq_name},$cds_name; 
    $names->{$cds_name} = "misc";
  } 
} 

while (<REF_FA>) { 
  ## change ref type to misc for all identified problematic CDS 
  if (m/^>(.*?).CDS/) { 
    my $cds_name = $1; 
    s/^>$cds_name\.CDS/>$cds_name.misc/ if (defined $names->{$cds_name}); 
    print; 
  } else { 
    print; 
  }
}  

close PROT_FA; 
close REF_FA;
