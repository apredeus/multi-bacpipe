#!/usr/bin/env perl 

## filter blast hits obtained by mapping external reference to all study genomes 
## write it out as a GFF for further testung 

use strict; 
use warnings; 
#use Data::Dumper; 

if (scalar @ARGV != 2) {
  print STDERR "Usage: ./blast_to_gff.pl <ref_blast_out> <modified_ref_fa>\n";
  exit 1
}

my $blast_out = shift @ARGV; 
my $ref_fa = shift @ARGV; 

open BLAST_OUT,"<",$blast_out or die "$!"; 
open REF_FA,"<",$ref_fa or die "$!"; 

my %length;
my $blast = {};

## read reference fasta (possibly folded), get length of each sequence 

my $seq_name;
 
while (<REF_FA>) {
  if (m/^>(.*?)\n/) {  
    $seq_name = $1;
  } else {
    $length{$seq_name} += length($_)-1;
  } 
}

## parse all blast hits, get rid of the poor quality ones
## filtering is exactly the same as will be done later 
## hits need to be unique **per replicon**!  
 
while (<BLAST_OUT>) { 
  my @t = split /\t+/;
  $t[0] =~ m/^(.*)\.(.*?)$/;
  my $name = $1; 
  my $type = $2;
  my $chr  = $t[1]; 
  die "ERROR: external reference can be of CDS/ncRNA/misc type only!" if ($type ne "CDS" && $type ne "ncRNA" && $type ne "misc");  
  my $len_ratio = $t[3]/$length{$t[0]}; 

  ## retain all high identity matches in blast hash 
  if ($t[2] > 90 && $len_ratio > 0.9) {   
    my $hit = (defined $blast->{$chr}->{$name}) ? scalar(keys(%{$blast->{$chr}->{$name}})) + 1 : 1;
 
    if ($t[9] > $t[8]) { 
      $blast->{$chr}->{$name}->{$hit}->{type} = $type; 
      $blast->{$chr}->{$name}->{$hit}->{chr} = $t[1]; 
      $blast->{$chr}->{$name}->{$hit}->{beg} = $t[8]; 
      $blast->{$chr}->{$name}->{$hit}->{end} = $t[9]; 
      $blast->{$chr}->{$name}->{$hit}->{strand} = "+";
 
      $blast->{$chr}->{$name}->{$hit}->{ident} = $t[2]; 
      $blast->{$chr}->{$name}->{$hit}->{len} = $len_ratio; 
    } else { 
      $blast->{$chr}->{$name}->{$hit}->{type} = $type; 
      $blast->{$chr}->{$name}->{$hit}->{chr} = $t[1]; 
      $blast->{$chr}->{$name}->{$hit}->{beg} = $t[9]; 
      $blast->{$chr}->{$name}->{$hit}->{end} = $t[8]; 
      $blast->{$chr}->{$name}->{$hit}->{strand} = "-";
 
      $blast->{$chr}->{$name}->{$hit}->{ident} = $t[2]; 
      $blast->{$chr}->{$name}->{$hit}->{len} = $len_ratio; 
    } 
  } 
}

## now parse the resulting hash and get rid of hits you don't want 
foreach my $chr (keys %{$blast})  { 
  foreach my $name (keys %{$blast->{$chr}}) { 
    if (scalar keys %{$blast->{$chr}->{$name}} > 1) { 
      ## 2+ hits per name
      my @hits = keys %{$blast->{$chr}->{$name}}; 
      my $best_ident = 0; 
      my $best_len = 0;
      ## establish best length and best id
      foreach my $hit (@hits) { 
        $best_ident = ($blast->{$chr}->{$name}->{$hit}->{ident} > $best_ident) ? $blast->{$chr}->{$name}->{$hit}->{ident} : $best_ident;
        $best_len = ($blast->{$chr}->{$name}->{$hit}->{len} > $best_len) ? $blast->{$chr}->{$name}->{$hit}->{len} : $best_len;
      } 
      ## remove non-qualifying hits from the blast hash
      foreach my $hit (@hits) { 
        my $ident = $blast->{$chr}->{$name}->{$hit}->{ident};
        my $len = $blast->{$chr}->{$name}->{$hit}->{len};
        delete $blast->{$chr}->{$name}->{$hit} if ($len < $best_len); ## if matches are not the same, keep longest
      }
    } 
  }
} 

## print final GFF here 
foreach my $chr (keys %{$blast})  {
  foreach my $name (keys %{$blast->{$chr}}) {
    foreach my $hit (keys %{$blast->{$chr}->{$name}}) {
      my $chr = $blast->{$chr}->{$name}->{$hit}->{chr};
      my $beg = $blast->{$chr}->{$name}->{$hit}->{beg};
      my $end = $blast->{$chr}->{$name}->{$hit}->{end};
      my $type = $blast->{$chr}->{$name}->{$hit}->{type};
      my $strand = $blast->{$chr}->{$name}->{$hit}->{strand};
      my $lt = join "_",$name,$chr,$beg;  
      printf "%s\tBacpipe\t%s\t%d\t%d\t.\t%s\t.\tID=%s;\n",$chr,$type,$beg,$end,$strand,$lt;
    }  
  } 
}

close BLAST_OUT; 
close REF_FA;
