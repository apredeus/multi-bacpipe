#!/usr/bin/env perl 

use strict; 
use warnings; 

my $nc_fa = shift @ARGV; 
my $locus_tag = shift @ARGV; 
my $blast_out = shift @ARGV; 


open FA,"<",$nc_fa or die "$!"; 
open BLAST,"<",$blast_out or die "$!"; 

my %length;
my $count=1;  

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

while (<BLAST>) { 
  my @t = split /\t+/; 
  #print "DEBUG: $t[3] $t[0] $length{$t[0]}\n"; 
  my $len_ratio = $t[3]/$length{$t[0]}; 
  if ($t[2] > 90 && $len_ratio > 0.9) {  ## we only want high identity matches
    my ($beg,$end,$strand); 
    ## ncRNA will have different numbering of IDs/locus tags  
    my $id = sprintf "%s_nc_%03d",$locus_tag,$count;
    if ($t[9] > $t[8]) { 
      $beg = $t[8]; 
      $end = $t[9]; 
      $strand = "+"; 
    } else { 
      $beg=$t[9]; 
      $end=$t[8]; 
      $strand = "-"; 
    }  
    printf "%s\tBlastn\tncRNA\t%d\t%d\t.\t%s\t.\tID=%s;Name=%s;gene_biotype=noncoding_rna;\n",$t[1],$beg,$end,$strand,$id,$t[0];  
    $count++;
  } 
}

close FA; 
close BLAST;  
