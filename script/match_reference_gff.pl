#!/usr/bin/env perl 

## similar to unify_strain_gff.pl, used on ref_strains; different in few ways: 
## - no new GFF file is produced, only blast results are filtered and matched to existing ref features 
## - match has to be pretty close to qualify and also of the same kind  
## in the end, $TAG.match.tsv is generated and later used by make_ortholog_table.pl  

use strict; 
use warnings; 
#use Data::Dumper; 

if (scalar @ARGV != 3) {
  print STDERR "Usage: ./match_reference_gff.pl <clean_gff> <ref_blast_out> <modified_ref_fa>\n";
  exit 1
}

my $clean_gff = shift @ARGV; 
my $blast_out = shift @ARGV; 
my $ref_fa = shift @ARGV; 
my $tag = $blast_out;
$tag =~ s/\.ref_blast.out//g;  

my $match_table = $tag.".match.tsv";

open GFF,"<",$clean_gff or die "$!"; 
open BLAST,"<",$blast_out or die "$!"; 
open FA,"<",$ref_fa or die "$!"; 
open MATCH,">",$match_table or die "$!"; 

my %length;
my $ref_locus = {};
my $blast = {};

## read reference fasta (make sure it's not folded!)
## get length of each sequence 
for (;;) {
   my $seq_name = <FA>;
   last if not defined $seq_name;
   my $seq = <FA>;
   last if not defined $seq;
   $seq_name =~ m/>(.*?)\n/; 
   my $seq_id = $1; 
   my $seq_len = length($seq)-1; 
   $length{$seq_id} = $seq_len; 
}

## parse all blast hits, get rid of the poor quality ones
## hash genes uses both names and prokka locus tags as keys 
## if there are multiple hits passing identity/length threshold, do the following: 
## - if hits are identical, keep all, add _2, _3 etc to gene name; 
## - if hits are not idential, keep only the longest (usually full length).
  
while (<BLAST>) { 
  my @t = split /\t+/;
  $t[0] =~ m/^(.*)\.(.*?)$/;
  my $name = $1; 
  my $type = $2;
  die "ERROR: external reference can be of CDS/ncRNA/misc type only!" if ($type ne "CDS" && $type ne "ncRNA" && $type ne "misc");  
  my $len_ratio = $t[3]/$length{$t[0]}; 

  ## retain all high identity matches in blast hash 
  if ($t[2] > 90 && $len_ratio > 0.9) {   
    my $hit = (defined $blast->{$name}) ? scalar(keys(%{$blast->{$name}})) + 1 : 1;
 
    if ($t[9] > $t[8]) { 
      $blast->{$name}->{$hit}->{type} = $type; 
      $blast->{$name}->{$hit}->{chr} = $t[1]; 
      $blast->{$name}->{$hit}->{beg} = $t[8]; 
      $blast->{$name}->{$hit}->{end} = $t[9]; 
      $blast->{$name}->{$hit}->{strand} = "+";
 
      $blast->{$name}->{$hit}->{ident} = $t[2]; 
      $blast->{$name}->{$hit}->{len} = $len_ratio; 
    } else { 
      $blast->{$name}->{$hit}->{type} = $type; 
      $blast->{$name}->{$hit}->{chr} = $t[1]; 
      $blast->{$name}->{$hit}->{beg} = $t[9]; 
      $blast->{$name}->{$hit}->{end} = $t[8]; 
      $blast->{$name}->{$hit}->{strand} = "-";
 
      $blast->{$name}->{$hit}->{ident} = $t[2]; 
      $blast->{$name}->{$hit}->{len} = $len_ratio; 
    } 
  } 
}

#print Dumper $blast;

foreach my $name (keys %{$blast}) { 
  ## 2+ hits per name require additional parsing, see above 
  if (scalar keys %{$blast->{$name}} != 1) { 
    my @hits = keys %{$blast->{$name}}; 
    my $best_ident = 0; 
    my $best_len = 0;
    ## establish best length and best identity %
    foreach my $hit (@hits) { 
      $best_ident = $blast->{$name}->{$hit}->{ident} if ($blast->{$name}->{$hit}->{ident} > $best_ident);
      $best_len = $blast->{$name}->{$hit}->{len} if ($blast->{$name}->{$hit}->{len} > $best_len);
    } 
    foreach my $hit (@hits) { 
      delete $blast->{$name}->{$hit} if ($blast->{$name}->{$hit}->{len} != $best_len);
    }
  } 
}

## clean GFF can only have 5 types of values in col 3: CDS, ncRNA, tRNA, rRNA, other.
## they do not have pseudogenes because we take the biggest CDS with that lt 
## that way Roary can make most use of it (otherwise the whole feature is dropped) 
 
while (<GFF>) {
  chomp; 
  my @t = split /\t+/;
  $t[8] =~ m/ID=(.*?);/;
  my $lt = $1; ## all entries in clean GFF would have an lt, so no checks 
  my $type = $t[2];     
  my $chr = $t[0];     
  my $beg = $t[3];     
  my $end = $t[4];
  my $strand = $t[6];
     
  foreach my $name (keys %{$blast}) { 
    foreach my $hit (keys %{$blast->{$name}}) {
      my $type_match = 0; 
      $type_match = 1 if ($type eq $blast->{$name}->{$hit}->{type} || $blast->{$name}->{$hit}->{type} eq "misc"); 
      if ($chr eq $blast->{$name}->{$hit}->{chr} && $strand eq $blast->{$name}->{$hit}->{strand} && $type_match) { 
        my $hit_beg = $blast->{$name}->{$hit}->{beg}; 
        my $hit_end = $blast->{$name}->{$hit}->{end};
        my $max_beg = ($beg > $hit_beg) ? $beg : $hit_beg;  
        my $min_end = ($end < $hit_end) ? $end : $hit_end; 
        my $ovl1 = ($min_end - $max_beg)/($end-$beg); 
        my $ovl2 = ($min_end - $max_beg)/($hit_end-$hit_beg); 
        if ($ovl1 > 0.5 && $ovl2 > 0.5) {                      ## you can adjust these later; I think they are ok for now  
          if (!defined $ref_locus->{$lt}->{name}) { 
            $ref_locus->{$lt}->{name} = $name;
            $ref_locus->{$lt}->{ident} =  $blast->{$name}->{$hit}->{ident};
            $ref_locus->{$lt}->{len} =  $blast->{$name}->{$hit}->{len};
          } else {  
            if ($blast->{$name}->{$hit}->{len} > $ref_locus->{$lt}->{len} || $blast->{$name}->{$hit}->{ident} > $ref_locus->{$lt}->{ident}) { 
              $ref_locus->{$lt}->{name} = $name;
              $ref_locus->{$lt}->{ident} =  $blast->{$name}->{$hit}->{ident};
              $ref_locus->{$lt}->{len} =  $blast->{$name}->{$hit}->{len};
            } 
          } 
        }
      }
    }
  }
}

foreach my $lt (keys %{$ref_locus}) { 
  printf MATCH "%s\t%s\n",$lt,$ref_locus->{$lt}->{name};
}
 
close FA;
close GFF;  
close BLAST; 
close MATCH;  
