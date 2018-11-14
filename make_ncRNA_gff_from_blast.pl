#!/usr/bin/env perl 

use strict; 
use warnings; 

my $nc_fa = shift @ARGV; 
my $cds_gff = shift @ARGV; 
my $locus_tag = shift @ARGV; 
my $blast_out = shift @ARGV; 

open FA,"<",$nc_fa or die "$!"; 
open GFF,"<",$cds_gff or die "$!"; 
open BLAST,"<",$blast_out or die "$!"; 

my %length;
my $cds = {};
my $count = 1;  

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

while (<GFF>) { 
  chomp; 
  my @tt = split /\t+/; 
  $tt[8] =~ m/ID=(.*?);/;
  my $lt = $1; 
  $cds->{$lt}->{chr} = $tt[0]; 
  $cds->{$lt}->{beg} = $tt[3]; 
  $cds->{$lt}->{end} = $tt[4]; 
  $cds->{$lt}->{strand} = $tt[6]; 
}

while (<BLAST>) { 
  my @t = split /\t+/; 
  #print "DEBUG: $t[3] $t[0] $length{$t[0]}\n"; 
  my $len_ratio = $t[3]/$length{$t[0]}; 
  if ($t[2] > 90 && $len_ratio > 0.9) {  ## we only want high identity matches
    my ($chr,$beg,$end,$strand); 
    ## ncRNA will have different numbering of IDs/locus tags  
    my $id = sprintf "%s_nc_%03d",$locus_tag,$count;
    if ($t[9] > $t[8]) { 
      $chr = $t[1]; 
      $beg = $t[8]; 
      $end = $t[9]; 
      $strand = "+"; 
    } else { 
      $chr = $t[1]; 
      $beg=$t[9]; 
      $end=$t[8]; 
      $strand = "-"; 
    } 

    ## I know this is a dumb way to do it, but oh well
    my $overlapped = 0;  
    foreach my $key (keys %{$cds}) {
      if ($chr eq $cds->{$key}->{chr} && $strand eq $cds->{$key}->{strand} && 
      $beg >= $cds->{$key}->{beg} && $end <= $cds->{$key}->{end} && $t[0] =~ m/mia-/) { 
        printf STDERR "Removing overlapping sORF: (%s) %s %d %d %s which overlaps CDS at %s %d %d %s\n",
        $t[0],$chr,$beg,$end,$strand,$cds->{$key}->{chr},$cds->{$key}->{beg},$cds->{$key}->{end},$cds->{$key}->{strand};
        $overlapped = 1;
      } 
    }

    if ($overlapped == 0) {     
      printf "%s\tBlastn\tncRNA\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s;gene_biotype=noncoding_rna;\n",$t[1],$beg,$end,$strand,$id,$t[0];  
      $count++;
    } 
  } 
}

close FA;
close GFF;  
close BLAST;  
