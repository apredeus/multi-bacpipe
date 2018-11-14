#!/usr/bin/env perl 

use strict; 
use warnings; 

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

my $gff_path = join '/',$wdir,'study_strains',$tag,$tag;
my $prokka_gff = join '.',$gff_path,"prokka/$tag.prokka.gff";  
my $ncrna_gff = $gff_path . ".ncRNA.gff";  

#print STDERR "DEBUG: $prokka_gff\n"; 
#print STDERR "DEBUG: $ncrna_gff\n"; 

open PROKKA_GFF,"<",$prokka_gff or die "$!"; 
open NCRNA_GFF,"<",$ncrna_gff or die "$!"; 
 
while (<PROKKA_GFF>) { 
  chomp; 
  next if (m/^#/ || m/^$/ || ! m/\t/); 
  
  if (m/\tCDS\t/) { 
    s/\tID=/%/; 
    my @tt = split /%/; 
    $tt[1] =~ m/^(.*?);/; 
    my $id = $1;
    my $name = (defined $names{$id}) ? $names{$id} : $id; 
    my $product = (defined $prods{$id}) ? $prods{$id} : "product=hypothetical protein;";  
    printf "%s\tID=%s;Name=%s;%sgene_biotype=protein_coding;\n",$tt[0],$id,$name,$product; 
  } else {
    ## if rRNA/tRNA/CRISPR, just print as is  
    print "$_\n"; 
  }  
}

while (<NCRNA_GFF>) {
  chomp;
  s/\tID=/%/;
  my @tt = split /%/;
  $tt[1] =~ m/^(.*?);/;
  my $id = $1;
  my $name = (defined $names{$id}) ? $names{$id} : $id; 
  my ($product,$biotype);
 
  if ($name =~ m/mia-/) {
    $product = "Small protein-coding ORF $name"; 
    $biotype = "protein_coding"; 
    $tt[0] =~ s/\tncRNA\t/\tCDS\t/g; 
  } else { 
    $product = "Non-coding RNA $name";  
    $biotype = "noncoding_rna"; 
  } 
  printf "%s\tID=%s;Name=%s;product=%s;gene_biotype=%s;\n",$tt[0],$id,$name,$product,$biotype; 
}

close PROKKA_GFF; 
close NCRNA_GFF; 

close EXP; 
close PRODS; 