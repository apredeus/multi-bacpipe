#!/usr/bin/env perl 

## LT2 example shows us that some annotations have 
## gene coordinates that are different from CDS coordinates 
## this should take care of such cases 

use strict;
use warnings; 

if ($#ARGV != 1) {
  die "USAGE: ncbi_to_roary.pl <ncbi_gff> <ncbi_fna>\n";
}

my $gff = shift @ARGV; 
my $fna = shift @ARGV; 

open GFF,"<",$gff or die "$!"; 
open FNA,"<",$fna or die "$!"; 

my $cds_coord = {};
my $cds_to_gene = {}; 
my @gff; 
my $cds_count = 0; 
my $gene_count = 0; 
my $pseudo_count = 0; 
my $rep_count = 0; 


while (<GFF>) { 
  if (m/\tCDS\t/) {
    my @tt = split /\t+/; 
    $tt[8] =~ m/Parent=(.*?);/; 
    my $gene = $1;
    ## the code below is done like this because NCBI is dumb
    ## and allows multiple identical gene IDs within the same GFF
    ## e.g. in SL1344 there's a gene0 on chr and each of the plasmids
    if (defined $cds_coord->{$tt[0]}->{$gene}) { 
      my @kk = split /\t+/,$cds_coord->{$tt[0]}->{$gene}; 
      my $len1 = $tt[4] - $tt[3]; 
      my $len2 = $kk[1] - $kk[0];
      ## if there are two fragments with the same ID, use the longest one  
      $cds_coord->{$tt[0]}->{$gene} = $tt[3]."\t".$tt[4] if ($len1 > $len2);
      $rep_count++; 
    } else { 
      $cds_coord->{$tt[0]}->{$gene} = $tt[3]."\t".$tt[4]; 
    } 
    $cds_count++; 
  } elsif (m/\tgene\t/ && m/gene_biotype=protein_coding;/) {
    push @gff,$_; 
    $gene_count++; 
  } elsif (m/\tpseudogene\t/) {
    ## we don't report pseudogenes in the final GFF  
    $pseudo_count++; 
  } 
} 

print STDERR "Processed GFF file: $gene_count protein-coding genes, $pseudo_count pseudogenes, $cds_count CDS entries, of them $rep_count repetitive (>1 CDS per 1 gene).\n"; 
print STDERR "Pseudogenes (entries with \"pseudogene\" in the 3rd column of GFF file) are not printed in the final Roary-friendly GFF.\n"; 
print STDERR "In case pseudogenes are listed as a \"gene\", there is often more than 1 CDS entry; This script would choose the longest fragment.\n"; 

foreach my $j (@gff) {
  my @tt = split /\t+/,$j; 
 
  $j =~ s/\tgene\t/\tCDS\t/; 
  $j =~ m/ID=(.*?);/;
  my $gene = $1; 
  $j =~ m/Name=(.*?);/;
  my $name = $1;
  ## some insanely moronic NCBI samples had old_locus_tag field etc  
  $j =~ m/;locus_tag=(\w+)/;
  my $lt = $1;
  ## if Name is the same as a locus tag, we don't want it 
  $j =~ s/Name=.*?;//g if ($name eq $lt); 
  ## replace gene ID with locus tag; 
  $j =~ s/ID=.*?;/ID=$lt;/; 
  ## replace gene coordinates with CDS coordinates - crucial for rare cases like LT2, when UTR are included in gene 
  $j =~ s/\t$tt[3]\t$tt[4]\t/\t$cds_coord->{$tt[0]}->{$gene}\t/;
  $j =~ s/part=1\/2;//;  
  print "$j" if ($j !~ m/part=2\/2/); 
} 

print "##FASTA\n";

while (<FNA>) { 
  if (m/^>/) { 
    my @tt = split /\s+/; 
    print "$tt[0]\n"; 
  } else { 
    print; 
  } 
} 

close GFF; 
close FNA;  
