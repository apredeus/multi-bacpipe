#!/usr/bin/env perl 

## this script is useful for generation of pangenome-based GFF file for each study strain

use strict; 
use warnings; 

if ($#ARGV != 2) { 
  die "USAGE: make_GFF_from_anno.pl <working_dir> <strain_tag> <products_file>\n";
} 

my $wdir = shift @ARGV; 
my $tag = shift @ARGV; 
my $products = shift @ARGV; 

my %prods; 
my %names; 

my $exp_table = join '/',$wdir,'exp_tables/Master_table.counts.tsv'; 

open EXP,"<",$exp_table or die "$!"; 
open PRODS,"<",$products or die "$!"; 

my $exp_header = <EXP>; 
my $col_num; 
my @ids = split /\t+/,$exp_header; 
for (my $i=0; $i < scalar @ids; $i++) { 
  $col_num = $i if ($ids[$i] eq $tag); 
} 

while (<EXP>) { 
  chomp; 
  my @tt = split /\t+/; 
  $names{$tt[$col_num]} = $tt[0]; 
} 

while (<PRODS>) { 
  chomp; 
  my @tt = split /\t+/; 
  $prods{$tt[0]} = $tt[1]; 
} 

my $united_gff = join '/',$wdir,'study_strains',$tag,"$tag.united.gff";

open UNITED_GFF,"<",$united_gff or die "$!"; 
 
while (<UNITED_GFF>) { 
  chomp; 
  my @t = split /\t/; 
  $t[8] =~ m/ID=(.*?);/;
  my $lt = $1;  
  my $name = "";
  my $product = ""; 
  $name = "Name=$names{$lt};" if (defined $names{$lt}); 
  $product = $prods{$lt} if (defined $prods{$lt}); 
  printf "%s%s%s\n",$_,$name,$product; 
}

close UNITED_GFF; 

close EXP; 
close PRODS; 
