#!/usr/bin/env perl 

## similar to reference_gff_cleanup.pl, but with some differences 
## basically all "gene" features but rRNA/tRNA are included, and "gene" borders are those defined by gene feature in GFF (eg with UTR for LT2)  
## features that have the same locus_tag (usually halves of a pseudogene) are UNIFIED into one feature that would be quantified as 1 entity
## "gene" entries with no locus_tag are GIVEN a new locus tag based on their ID (usually gene*) - now LT looks like NOTAG_gene123 
## "ncRNA" are all things that match *rna (non-case-spec), but not rRNA/tRNA. 

use strict;
use warnings; 

if ($#ARGV != 0) {
  die "USAGE: simple_gff_cleanup.pl <ncbi_gff>\n";
}

my $gff = shift @ARGV; 

open GFF,"<",$gff or die "$!"; 

my $genes = {};
my $prod = {};
my ($coding,$pseudo,$ncrna,$trna,$rrna,$other,$notag) = (0)x7; 


while (<GFF>) { 
  if (m/\tgene\t/) {
    my @t = split /\t+/; 
    my $id = ($t[8] =~ m/ID=(.*?);/) ? $1 : "NONE"; 
    my $name = ($t[8] =~ m/Name=(.*?);/) ? $1 : "NONE"; 
    my $biotype = ($t[8] =~ m/;gene_biotype=(\w+)/) ? $1 : "NONE";
    ## quite few annotations have ncRNAs and rRNA/tRNA without a locus tag
    my $lt   = ($t[8] =~ m/;locus_tag=(\w+)/) ? $1 : "NOTAG_".$id; 
    $notag++ if ($lt =~ m/_$id$/); 
    ## a bit of extra safety 
    $biotype = "pseudogene" if ($t[8] =~ m/pseudo=true/); 
    ## this accounts for duplicate locus tags in 
    if (! defined $genes->{$lt}) { 
      $genes->{$lt}->{lt} = $lt;  
      $genes->{$lt}->{id} = $id;  
      $genes->{$lt}->{name} = $name;  
      $genes->{$lt}->{biotype} = $biotype;
      $genes->{$lt}->{chr} = $t[0];   
      $genes->{$lt}->{beg} = $t[3];   
      $genes->{$lt}->{end} = $t[4];   
      $genes->{$lt}->{strand} = $t[6];
    } else { 
      $genes->{$lt}->{beg} = ($genes->{$lt}->{beg} <= $t[3]) ? $genes->{$lt}->{beg} : $t[3];  
      $genes->{$lt}->{end} = ($genes->{$lt}->{end} >= $t[4]) ? $genes->{$lt}->{end} : $t[4];  
      ## everything else is defined and stays the same 
    } 
  } elsif (m/Parent=gene/) {
    ## products are either in product or Note field; use former, if not there, use latter
    my @t = split /\t+/; 
    my $id = ($t[8] =~ m/Parent=(.*?);/) ? $1 : "NONE"; 
    my $product = ($t[8] =~ m/;product=(.*?);/) ? $1 : "NONE"; 
    my $note = ($t[8] =~ m/;Note=(.*?);/) ? $1 : "NONE"; 

    ## gene0 can be found on chr and plasmid sometimes, so
    $prod->{$t[0]}->{$id}->{product} = $product; 
    $prod->{$t[0]}->{$id}->{note} = $note; 
  } elsif (m/\tpseudogene\t/) {
    my @t = split /\t+/; 
    my $id = ($t[8] =~ m/ID=(.*?);/) ? $1 : "NONE"; 
    my $name = ($t[8] =~ m/Name=(.*?);/) ? $1 : "NONE"; 
    my $biotype = "pseudogene";
    my $lt   = ($t[8] =~ m/;locus_tag=(\w+)/) ? $1 : "NOTAG_".$id; 
    $notag++ if ($lt =~ m/_$id$/); 

    ## this accounts for duplicate locus tags in 
    if (!defined $genes->{$lt}) { 
      $genes->{$lt}->{lt} = $lt;  
      $genes->{$lt}->{id} = $id;  
      $genes->{$lt}->{name} = $name;  
      $genes->{$lt}->{biotype} = $biotype;
      $genes->{$lt}->{chr} = $t[0];   
      $genes->{$lt}->{beg} = $t[3];   
      $genes->{$lt}->{end} = $t[4];   
      $genes->{$lt}->{strand} = $t[6];
    } else {
      $genes->{$lt}->{beg} = ($genes->{$lt}->{beg} <= $t[3]) ? $genes->{$lt}->{beg} : $t[3];  
      $genes->{$lt}->{end} = ($genes->{$lt}->{end} >= $t[4]) ? $genes->{$lt}->{end} : $t[4]; 
      ## everything else is defined and stays the same 
    } 
  } 
} 

print STDERR "$gff annotation file processed; found $notag gene entries without a locus tag, for which new locus tags were generated.\n";

my @keys = sort { $genes->{$a}->{chr} cmp $genes->{$b}->{chr} || $genes->{$a}->{beg} <=> $genes->{$b}->{beg} } keys %{$genes};

foreach my $lt (@keys) {
  if ($lt ne "NONE") { 
    my $out = sprintf "%s\tBacpipe\tgene\t%s\t%s\t.\t%s\t.\t",$genes->{$lt}->{chr},$genes->{$lt}->{beg},$genes->{$lt}->{end},$genes->{$lt}->{strand};
    my $gene_length = $genes->{$lt}->{end} - $genes->{$lt}->{beg}; 
    $out = join ('',$out,"ID=",$lt,";");
    my $name = ($genes->{$lt}->{name} eq "NONE") ? $lt : $genes->{$lt}->{name}; 
    $out = join ('',$out,"Name=",$name,";");
    my $product = (defined $prod->{$genes->{$lt}->{chr}}->{$genes->{$lt}->{id}}->{product}) ? $prod->{$genes->{$lt}->{chr}}->{$genes->{$lt}->{id}}->{product} : "NONE"; 
    my $note    = (defined $prod->{$genes->{$lt}->{chr}}->{$genes->{$lt}->{id}}->{note}) ? $prod->{$genes->{$lt}->{chr}}->{$genes->{$lt}->{id}}->{note} : "NONE"; 
    
    ## if product is defined, then use product; if note is defined instead of product, use note;
    ## if neither is defined, use nothing. 
    my $biotype = $genes->{$lt}->{biotype}; 
    if ($biotype =~ m/rna/i && $biotype ne "rRNA" && $biotype ne "tRNA") { 
      $biotype = "noncoding_rna"; 
    }
    if ($biotype eq "protein_coding") { 
      $coding++; 
    } elsif ($biotype eq "pseudogene") { 
      $pseudo++;
    } elsif ($biotype eq "noncoding_rna") { 
      $ncrna++; 
    } elsif ($biotype eq "rRNA") { 
      $rrna++; 
    } elsif ($biotype eq "tRNA") { 
      $trna++; 
    } else { 
      $other++; 
    } 
    
    $out = join ('',$out,"gene_biotype=",$biotype,";");
   
    if ($product ne "NONE") { 
      $out = join ('',$out,"product=",$product,";");
    } elsif ($note ne "NONE") { 
      $out = join ('',$out,"product=",$note,";");
    } 
    $out = join ('',$out,"locus_tag=",$lt,"\n");
    ## print out every gene entry, even rRNA and tRNA (they won't have expression since reads are removed from BAM)
    ## merge pseudogene halves into one "gene" entry with the same locus tag if LT of halves are the same 
    ## above treatment of pseudogenes ONLY APPLIES WHEN --SIMPLE IS IN USE! 
    if ($gene_length <= 50000 && $biotype ne "rRNA" && $biotype ne "tRNA") { 
      print $out; 
    } elsif ($gene_length > 50000) { 
      print STDERR "WARNING: Skipping gene $lt; length $gene_length is over 50000 - possible annotation bug?\n"; 
    } elsif ($biotype eq "rRNA" || $biotype eq "tRNA") { 
      print STDERR "WARNING: Skipping gene $lt (rRNA/tRNA).\n"; 
    } 
  }
} 

print STDERR "Parsed GFF stats: $coding protein coding, $pseudo pseudogenes, $ncrna noncoding RNAs, $trna tRNAs, $rrna rRNAs, $other others.\n"; 

close GFF; 
