#!/usr/bin/env perl 

## LT2 example shows us that some annotations have 
## gene coordinates that are different from CDS coordinates 
## this should take care of such cases 

use strict;
use warnings; 

my $gff = shift @ARGV; 
my $fna = shift @ARGV; 
my $tag = shift @ARGV; 

open GFF,"<",$gff or die "$!"; 
open FNA,"<",$fna or die "$!"; 
open RRY,">","$tag.gff" or die "$!"; 

my $cds_coord={};
my @gff; 

while (<GFF>) { 
  if (m/\tCDS\t/) {
    my @tt = split /\t+/; 
    $tt[8] =~ m/Parent=(.*?);/; 
    my $gene = $1;
    ## this is done like this because NCBI is dumb
    ## and allows multiple identical gene IDs within the same GFF
    ## e.g. in SL1344 there's a gene0 on chr and each of the plasmids
    $cds_coord->{$tt[0]}->{$gene} = $tt[3]."\t".$tt[4]; 
  } elsif (m/\tgene\t/ && m/gene_biotype=protein_coding;/) { 
    push @gff,$_; 
  } 
} 

foreach my $j (@gff) {
  my @tt = split /\t+/,$j; 
 
  $j =~ s/\tgene\t/\tCDS\t/; 
  $j =~ m/ID=(.*?);.*Name=(.*?);.*locus_tag=(.*)/; 
  my $gene = $1; 
  my $name = $2; 
  my $lt = $3;
  ## if Name is the same as a locus tag, we don't want it 
  $j =~ s/Name=.*?;//g if ($name eq $lt); 
  ## replace gene ID with locus tag; 
  $j =~ s/ID=.*?;/ID=$lt;/; 
  ## replace gene coordinates with CDS coordinates
  $j =~ s/\t$tt[3]\t$tt[4]\t/\t$cds_coord->{$tt[0]}->{$gene}\t/; 
  print RRY "$j"; 
} 

print RRY "##FASTA\n";

while (<FNA>) { 
  if (m/^>/) { 
    my @tt = split /\s+/; 
    print RRY "$tt[0]\n"; 
  } else { 
    print RRY; 
  } 
} 

close GFF; 
close FNA;  
close RRY; 
