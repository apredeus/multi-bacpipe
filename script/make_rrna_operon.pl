#!/usr/bin/env perl 

## take united GFF file - united Blastn/Prokka output 
## return rRNA operons and singleton tRNAs as BED intervals
## operon is extended until 50 bp before nearest non-rRNA/non-tRNA feature  

use strict; 
use warnings; 
#use Data::Dumper; 

if ($#ARGV != 0) {
  die "USAGE: make_rrna_operon.pl <united_gff>\n";
}

my $gff   = shift @ARGV; 
open GFF,"<",$gff or die "$!"; 

my $genes = {}; ## 
my @rrna; 

while (<GFF>) {
  if (! m/\t/) {
    next; 
  } else { 
    chomp; 
    my @t = split /\t+/; 
    $t[8] =~ m/ID=(.*?);/; 
    my $id = $1; 
    $genes->{$id}->{chr} = $t[0]; 
    $genes->{$id}->{beg} = $t[3]-1; 
    $genes->{$id}->{end} = $t[4];

    if ($t[2] eq "tRNA") {  
      $genes->{$id}->{type} = "tRNA"; 
    } elsif ($t[2] eq "rRNA") {  
      $genes->{$id}->{type} = "rRNA"; 
      push @rrna,$id; 
    } else { 
      $genes->{$id}->{type} = "other"; 
    } 
  } 
}

## print Dumper $genes; 

foreach my $rrna (@rrna) { 
  ## distance from nearest non-rRNA, non-tRNA gene to this rRNA gene from the left/right 
  my $d_left  = 1000000; 
  my $d_right = 1000000; 
  my $rbeg = $genes->{$rrna}->{beg};
  my $rend = $genes->{$rrna}->{end};
  
  foreach my $key (keys %{$genes}) { 
    if ($genes->{$key}->{chr} eq $genes->{$rrna}->{chr} && $genes->{$key}->{type} ne "rRNA" && $genes->{$key}->{type} ne "tRNA") { 
       my $cbeg = $genes->{$key}->{beg}; 
       my $cend = $genes->{$key}->{end}; 

       $d_left  = $rbeg-$cend if ($cend < $rbeg && $rbeg-$cend < $d_left);  
       $d_right = $cbeg-$rend if ($cbeg > $rend && $cbeg-$rend < $d_right); 
    } 
  }
  
  $genes->{$rrna}->{lbound} = ($d_left > 50) ? $genes->{$rrna}->{beg} - $d_left + 50 : $genes->{$rrna}->{beg}; 
  $genes->{$rrna}->{rbound} = ($d_right > 50) ? $genes->{$rrna}->{end} + $d_right - 50 : $genes->{$rrna}->{end}; 
}

foreach my $key (keys %{$genes}) { 
  if ($genes->{$key}->{type} eq "rRNA") {
    printf "%s\t%d\t%d\trRNA.%s\n",$genes->{$key}->{chr},$genes->{$key}->{lbound},$genes->{$key}->{rbound},$key;
  } elsif ($genes->{$key}->{type} eq "tRNA") { 
    printf "%s\t%d\t%d\ttRNA.%s\n",$genes->{$key}->{chr},$genes->{$key}->{beg},$genes->{$key}->{end},$key;
  } 
} 

close GFF; 
