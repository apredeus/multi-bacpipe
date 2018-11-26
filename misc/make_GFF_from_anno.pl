#!/usr/bin/env perl 

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
    my $product = (defined $prods{$id}) ? $prods{$id} : "note=hypothetical protein;";  
    printf "%s\tID=%s;Name=%s;%sgene_biotype=protein_coding;\n",$tt[0],$id,$name,$product; 
  } elsif (m/\ttmRNA\t/) {
    s/product=/note=/; 
    print "$_\n"; 
  } elsif (m/note=CRISPR/) { 
    s/\tnote=/\tName=CRISPR;note=/g; 
    print "$_\n"; 
  } else {
    my $name;
    m/\tID=(.*?);/;
    my $lt = $1; 
    if (m/product=16S/) { 
      $name="16S"; 
    } elsif (m/product=23S/) { 
      $name="23S"; 
    } elsif (m/product=5S/) {
      $name="5S"; 
    } elsif (m/product=(tRNA-.*?)\(/) {
      $name="$1"; 
    } 
    s/product=/note=/; 
    s/ID=$lt/ID=$lt;Name=$name/;
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
  printf "%s\tID=%s;Name=%s;note=%s;gene_biotype=%s;\n",$tt[0],$id,$name,$product,$biotype; 
}

close PROKKA_GFF; 
close NCRNA_GFF; 

close EXP; 
close PRODS; 
