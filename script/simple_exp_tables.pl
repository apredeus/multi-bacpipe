#!/usr/bin/env perl 

use strict; 
use warnings; 

if (scalar @ARGV != 2) {
  print "Usage: ./simple_exp_tables.pl <full_wdir> <tag>\n";
  exit 1
}

my $wdir = shift @ARGV;
my $tag  = shift @ARGV;

my $prophage = join "",$wdir,"/study_strains/",$tag,"/",$tag,".prophage.bed";
my $gene_gff = join "",$wdir,"/study_strains/",$tag,"/",$tag,".gene.gff";
my $count_table = join "",$wdir,"/exp_tables/",$tag,".counts.tsv";
my $tpm_table = join "",$wdir,"/exp_tables/",$tag,".TPM.tsv";

my $count_out = join "",$wdir,"/exp_tables/",$tag,".annotated.counts.tsv";
my $tpm_out = join "",$wdir,"/exp_tables/",$tag,".annotated.TPM.tsv";

if (!-f $gene_gff || !-f $count_table || !-f $tpm_table) { 
  die "ERROR: One of the required files - $tag.gene.gff, $tag.counts.tsv, or $tag.TPM.tsv - is not found!"
} 

if (-f $prophage) { 
  print "Found interval file $tag.prophage.bed, will use it in annotation.\n";
} else { 
  print "Interval file $tag.prophage.bed not found, will NOT annotate chr/prophage/plasmid!\n";
} 

open GENE_GFF,"<",$gene_gff or die "$!"; 
open COUNT,"<",$count_table or die "$!"; 
open TPM,"<",$tpm_table or die "$!"; 
open COUNT_OUT,">",$count_out or die "$!"; 
open TPM_OUT,">",$tpm_out or die "$!"; 

my $genes = {}; 
my $cds_loc = {}; 

## parse our gene.gff file, get gene type and name
while (<GENE_GFF>) { 
  chomp; 
  my @t = split /\t+/; 
  $t[8] =~ m/ID=(.*?);Name=(.*?);gene_biotype=(.*?);/; 
  my $id = $1; 
  my $name = $2; 
  my $type = $3; ## don't change if == other, rRNA, or tRNA
  if ($type eq "protein_coding") {
    $type = "CDS"; 
  } elsif ($type eq "noncoding_rna") {
    $type = "ncRNA"; 
  } elsif ($type eq "pseudogene") { 
    $type = "pseudo";
  } 
  $genes->{$id}->{chr} = $t[0];  
  $genes->{$id}->{beg} = $t[3]-1;  ## oh come on   
  $genes->{$id}->{end} = $t[4];  
  $genes->{$id}->{name} = $name;  
  $genes->{$id}->{type} = $type;  
} 

## parse interval file
if (-f $prophage) { 
  my $fasta_index = join "",$wdir,"/study_strains/",$tag,"/",$tag,".genome.fa.fai";
  my $chr_name = `cat $fasta_index | sort -k2,2nr | head -n 1 | cut -f 1`;
  chomp $chr_name;
  $cds_loc->{$tag}->{chr_name} = $chr_name;
  open PROP,"<",$prophage;
  while (<PROP>) {
    chomp;
    my @t = split /\t+/;
    $cds_loc->{$tag}->{$t[3]}->{chr}=$t[0];
    $cds_loc->{$tag}->{$t[3]}->{beg}=$t[1];
    $cds_loc->{$tag}->{$t[3]}->{end}=$t[2];
  }
  close PROP;
} 

my $count_header = <COUNT>; 
my $tpm_header = <TPM>; 
$count_header =~ s/$tag\t//g; 

if (-f $prophage) {  
  print COUNT_OUT "Locus_tag\tGene_name\tGene_type\tLocation\t$count_header";
  print TPM_OUT   "Locus_tag\tGene_name\tGene_type\tLocation\t$count_header";
} else { 
  print COUNT_OUT "Locus_tag\tGene_name\tGene_type\t$count_header";
  print TPM_OUT   "Locus_tag\tGene_name\tGene_type\t$count_header";
} 

while (<COUNT>) {
  my $counts = $_; 
  $counts =~ s/^(.*?)\t//; 
  my $id = $1; 

  if (-f $prophage) {
    my $chr = $genes->{$id}->{chr};
    my $beg = $genes->{$id}->{beg};
    my $end = $genes->{$id}->{end};
    my $name = $genes->{$id}->{name};
    my $type = $genes->{$id}->{type};
    ## if replicon name is not that of chromosome, it's one of the plasmids
    my $gene_loc = "chromosome";
    $gene_loc = "plasmid" if ($chr ne $cds_loc->{$tag}->{chr_name});
    foreach my $key (keys %{$cds_loc->{$tag}}) {
      if ($key ne "chr_name" &&  $cds_loc->{$tag}->{$key}->{chr} eq $chr) {
        my $prop_beg = $cds_loc->{$tag}->{$key}->{beg};
        my $prop_end = $cds_loc->{$tag}->{$key}->{end};
        $gene_loc = "prophage" if (($prop_beg <= $beg && $prop_end >= $beg) || ($prop_beg <= $end && $prop_end >= $end));
      }
    }
    printf COUNT_OUT "%s\t%s\t%s\t%s\t%s",$id,$name,$type,$gene_loc,$counts; 
  } else { 
    my $chr = $genes->{$id}->{chr};
    my $name = $genes->{$id}->{name};
    my $type = $genes->{$id}->{type};
    printf COUNT_OUT "%s\t%s\t%s\t%s",$id,$name,$type,$counts; 
  } 
}

while (<TPM>) {
  my $tpms = $_; 
  $tpms =~ s/^(.*?)\t//; 
  my $id = $1;
 
  if (-f $prophage) {
    my $chr = $genes->{$id}->{chr};
    my $beg = $genes->{$id}->{beg};
    my $end = $genes->{$id}->{end};
    my $name = $genes->{$id}->{name};
    my $type = $genes->{$id}->{type};
    ## if replicon name is not that of chromosome, it's one of the plasmids
    my $gene_loc = "chromosome";
    $gene_loc = "plasmid" if ($chr ne $cds_loc->{$tag}->{chr_name});
    foreach my $key (keys %{$cds_loc->{$tag}}) {
      if ($key ne "chr_name" &&  $cds_loc->{$tag}->{$key}->{chr} eq $chr) {
        my $prop_beg = $cds_loc->{$tag}->{$key}->{beg};
        my $prop_end = $cds_loc->{$tag}->{$key}->{end};
        $gene_loc = "prophage" if (($prop_beg <= $beg && $prop_end >= $beg) || ($prop_beg <= $end && $prop_end >= $end));
      }
    }
    printf TPM_OUT "%s\t%s\t%s\t%s\t%s",$id,$name,$type,$gene_loc,$tpms; 
  } else { 
    my $chr = $genes->{$id}->{chr};
    my $name = $genes->{$id}->{name};
    my $type = $genes->{$id}->{type};
    printf TPM_OUT "%s\t%s\t%s\t%s",$id,$name,$type,$tpms; 
  } 
} 

close GENE_GFF; 
close COUNT; 
close TPM; 
close COUNT_OUT; 
close TPM_OUT; 
