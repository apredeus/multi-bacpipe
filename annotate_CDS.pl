#!/usr/bin/env perl 

## script for parsing Roary CDS orthology output
## and for combining it with prophage information 
## we aim to make the final table that would include 
## 1) name (not necessarily unique), 2) loc (chr/prophage/plasmid); 
## 3-n) all IDs (locus tags) for N used + n reference strains. 

use strict; 
use warnings; 
use Data::Dumper; 

my $pres_abs = shift @ARGV; 
my $refdir = shift @ARGV; 
my $config = shift @ARGV; 

open ROARY,"<",$pres_abs or die "$!"; 
open CONFIG,"<",$config or die "$!"; 

my $cds_loc = {}; 
my $cds_coord = {}; 
my @ref_strains; 
my @study_strains; 

while (<CONFIG>) { 
  ## define prophage ranges for every strain
  ## also, define what is the chromosome 
  chomp; 
  my $strain = (split /\t+/)[1];
  if (! defined $cds_loc->{$strain} && $_ !~ m/^Reference\t/) { 
    my $prophage = join "",$refdir,"/",$strain,"/",$strain,".prophage.bed";
    my $fasta_index = join "",$refdir,"/",$strain,"/",$strain,".genome.fa.fai";
    my $gff = join "",$refdir,"/",$strain,"/",$strain,".CDS.gff";
    my $chr_name = `cat $fasta_index | sort -k2,2nr | head -n 1 | cut -f 1`;
    chomp $chr_name;
    $cds_loc->{$strain}->{chr_name}=$chr_name; 
    push @study_strains,$strain; 

    open PROP,"<",$prophage; 
    while (<PROP>) { 
      chomp;
      my @tt = split /\t+/;
      $cds_loc->{$strain}->{$tt[3]}->{chr}=$tt[0];
      $cds_loc->{$strain}->{$tt[3]}->{beg}=$tt[1];
      $cds_loc->{$strain}->{$tt[3]}->{end}=$tt[2];
    }  
    close PROP; 

    open GFF,"<",$gff; 
    while (<GFF>) { 
      chomp; 
      my @tt = split /\t+/; 
      $tt[8] =~ m/ID=(.*?);/; 
      my $locus_tag = $1;
      ## all locus_tags are unique (confirmed in Roary) 
      $cds_coord->{$locus_tag}->{chr}=$tt[0]; 
      $cds_coord->{$locus_tag}->{beg}=$tt[3]; 
      $cds_coord->{$locus_tag}->{end}=$tt[4]; 
    } 
    close GFF;  
  } elsif (! defined $cds_loc->{$strain} && $_ =~ m/^Reference\t/) {
      $cds_loc->{$strain}->{type}="reference"; 
      push @ref_strains,$strain; 
  } 
} 

## both have to be sorted alphabetically to make strain order reproducible 
@ref_strains = sort { $a cmp $b } @ref_strains; 
@study_strains = sort { $a cmp $b } @study_strains;

my $header = <ROARY>;
chomp $header;  
my @names = split /,/,$header; 

for (my $i = 0; $i < scalar @names; $i++) {
  $names[$i] =~ s/"//g; 
  $cds_loc->{$names[$i]}->{index}=$i if (defined $cds_loc->{$names[$i]});
}   
   
#print Dumper($cds_loc);
while (<ROARY>) {
  chomp; 
  my @tt = split /,/;
  foreach my $tt (@tt) { 
    $tt =~ s/"//g; 
  } 
  my $gene_name = ($tt[1] eq "") ? $tt[0] : $tt[1];
  my $gene_loc="chromosome";                ## chromosome, plasmid, prophage 
  my $output = ""; 
  foreach my $strain (@study_strains) { 
    my $index = $cds_loc->{$strain}->{index};
    my $locus_tag = ($tt[$index] eq "") ? "NONE" : $tt[$index]; 
    if ($locus_tag ne "NONE") {
      my $chr = $cds_coord->{$locus_tag}->{chr};
      my $beg = $cds_coord->{$locus_tag}->{beg};
      my $end = $cds_coord->{$locus_tag}->{end};
      print STDERR "LT is $locus_tag\n" if (! defined $chr); 
      $gene_loc = "plasmid" if ($chr ne $cds_loc->{$strain}->{chr_name});
      foreach my $key (keys %{$cds_loc->{$strain}}) { 
        if ($key ne "chr_name" && $key ne "index" &&  $cds_loc->{$strain}->{$key}->{chr} eq $chr) {
          my $prop_beg = $cds_loc->{$strain}->{$key}->{beg};
          my $prop_end = $cds_loc->{$strain}->{$key}->{end};
          $gene_loc = "prophage" if (($prop_beg <= $beg && $prop_end >= $beg) || ($prop_beg <= $end && $prop_end >= $end));
        }
      }
    }
    $output = join "\t",$output,$locus_tag;
  }  
  foreach my $strain (@ref_strains) {
    my $index = $cds_loc->{$strain}->{index}; 
    #print "DEBUG1: $strain $index\n"; 
    my $locus_tag = ($tt[$index] eq "") ? "NONE" : $tt[$index]; 
    $output = join "\t",$output,$locus_tag;
  }  
  printf "%s\t%s%s\n",$gene_name,$gene_loc,$output; 
} 

close ROARY; 
close CONFIG;  
