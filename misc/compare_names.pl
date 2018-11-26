#!/usr/bin/env perl 

use strict; 
use warnings; 

## compare names in a given reference with names assigned by Roary 
## report discrepancies, fix if necessary 

if ($#ARGV != 1) { 
  die "USAGE: compare_names.pl <ref_gff> <master_exp_table>\n";
} 

my $ref_gff = shift @ARGV; 
my $exp_table = shift @ARGV; 

open GFF,"<",$ref_gff or die "$!"; 
open EXP,"<",$exp_table or die "$!"; 

my $tag = $ref_gff; 
$tag =~ s/.gff//; 
my %ref_name; 

while (<GFF>) {
  next if (m/^#/ || ! m/\t/); 
  if (m/ID=(.*?);Name=(.*?);/) { 
    $ref_name{$1} = $2; 
    ## print "DEBUG: assigned name $2 to locus tag $1\n"; 
  } 
} 

my $exp_header = <EXP>; 
my $col_num; 
my @ids = split /\t+/,$exp_header; 
for (my $i=0; $i < scalar @ids; $i++) { 
  $col_num = $i if ($ids[$i] eq $tag); 
}

while (<EXP>) { 
  chomp; 
  my @tt = split /\t+/; 
  my $lt = $tt[$col_num];
  if (defined $ref_name{$lt} && $ref_name{$lt} ne $tt[0] && $tt[0] !~ m/group_/) {  
    print "LT: $lt Reference: $ref_name{$lt}, roary: $tt[0]\n"; 
  } 
}

close EXP; 
close GFF;   
