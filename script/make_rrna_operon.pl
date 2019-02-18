#!/usr/bin/env perl 

## take two GFF files, Prokka & ncRNA 
## return rRNA operons and singleton tRNAs as BED intervals
## operon is extended until 50 bp before nearest non-rRNA/non-tRNA feature  

use strict; 
use warnings; 
#use Data::Dumper; 

if ($#ARGV != 1) {
  die "USAGE: make_rrna_operon.pl <prokka_gff> <gene_gff>\n";
}

my $prokka_gff = shift @ARGV; 
my $gene_gff   = shift @ARGV; 

open PRKGFF,"<",$prokka_gff or die "$!"; 
open GENGFF,"<",$gene_gff or die "$!"; 

my $gene = {}; ## 
my $crispr_count = 1; 
my @rrna; 

while (<PRKGFF>) {
  if (! m/\t/) {
    next; 
  } elsif (m/CRISPR/) {
    chomp; 
    my @tt = split /\t+/; 
    my $id = "CRISPR_" . $crispr_count; 
    $gene->{$id}->{chr} = $tt[0]; 
    $gene->{$id}->{beg} = $tt[3]-1; 
    $gene->{$id}->{end} = $tt[4];
    $gene->{$id}->{type} = "other"; 
    $crispr_count++; 
  } else { 
    chomp; 
    my @tt = split /\t+/; 
    $tt[8] =~ m/ID=(.*?);/; 
    my $id = $1; 
    $gene->{$id}->{chr} = $tt[0]; 
    $gene->{$id}->{beg} = $tt[3]-1; 
    $gene->{$id}->{end} = $tt[4];
    if ($tt[2] eq "CDS") {
      $gene->{$id}->{type} = "CDS"; 
    } elsif ($tt[2] eq "tRNA") {  
      $gene->{$id}->{type} = "tRNA"; 
    } elsif ($tt[2] eq "rRNA") {  
      $gene->{$id}->{type} = "rRNA"; 
      push @rrna,$id; 
    } else { 
      $gene->{$id}->{type} = "other"; 
    } 
  } 
}

while (<GENGFF>) { 
  chomp; 
  my @tt = split /\t+/; 
  $tt[8] =~ m/ID=(.*?);/;
  my $id = $1;
  if ($tt[8] !~ m/gene_biotype=rRNA;/ && $tt[8] !~ m/gene_biotype=tRNA;/) { 
    $tt[8] =~ m/gene_biotype=(.*?);/;
    $gene->{$id}->{chr} = $tt[0];
    $gene->{$id}->{beg} = $tt[3]-1;
    $gene->{$id}->{end} = $tt[4];
    $gene->{$id}->{type} = $1;
  }  
}

## print Dumper $gene; 

foreach my $rrna (@rrna) { 
  ## distance from nearest non-rRNA, non-tRNA gene to this rRNA gene from the left/right 
  my $d_left  = 1000000; 
  my $d_right = 1000000; 
  my $rbeg = $gene->{$rrna}->{beg};
  my $rend = $gene->{$rrna}->{end};
  
  foreach my $key (keys %{$gene}) { 
    if ($gene->{$key}->{chr} eq $gene->{$rrna}->{chr} && $gene->{$key}->{type} ne "rRNA" && $gene->{$key}->{type} ne "tRNA") { 
       my $cbeg = $gene->{$key}->{beg}; 
       my $cend = $gene->{$key}->{end}; 

       $d_left  = $rbeg-$cend if ($cend < $rbeg && $rbeg-$cend < $d_left);  
       $d_right = $cbeg-$rend if ($cbeg > $rend && $cbeg-$rend < $d_right); 
    } 
  }
  
  $gene->{$rrna}->{lbound} = ($d_left > 50) ? $gene->{$rrna}->{beg} - $d_left + 50 : $gene->{$rrna}->{beg}; 
  $gene->{$rrna}->{rbound} = ($d_right > 50) ? $gene->{$rrna}->{end} + $d_right - 50 : $gene->{$rrna}->{end}; 
}

foreach my $key (keys %{$gene}) { 
  if ($gene->{$key}->{type} eq "rRNA") {
    printf "%s\t%d\t%d\trRNA.%s\n",$gene->{$key}->{chr},$gene->{$key}->{lbound},$gene->{$key}->{rbound},$key;
  } elsif ($gene->{$key}->{type} eq "tRNA") { 
    printf "%s\t%d\t%d\ttRNA.%s\n",$gene->{$key}->{chr},$gene->{$key}->{beg},$gene->{$key}->{end},$key;
  } 
} 

close PRKGFF; 
close GENGFF; 
