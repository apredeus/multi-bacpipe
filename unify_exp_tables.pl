#!/usr/bin/env perl 

use strict; 
use warnings; 
use Data::Dumper; 

## make master expression table with Roary ortholog table 
## this also uses phages annotation from PHAST and GFF files 
## to annotate genes as belonging to chromosome, prophage, or a plasmid

my $map = shift @ARGV; 
my $bed1 = shift @ARGV; 
my $bed2 = shift @ARGV; 
my $gff1 = shift @ARGV; 
my $gff2 = shift @ARGV; 
my $expmat1 = shift @ARGV; ## assume bacpipe format - 3-col annotation for genes
my $expmat2 = shift @ARGV; 

my $phage1={}; 
my $phage2={};
my $loc1={}; 
my $loc2={};
my $exp1={};
my $exp2={};

my $chr_name1="D23580_liv_chro";  # I'll cheat here  
my $chr_name2="CP002487.fasta"; 

open  MAP,"<",$map  or die "$!";   
open BED1,"<",$bed1 or die "$!";   
open BED2,"<",$bed2 or die "$!";   
open GFF1,"<",$gff1 or die "$!";   
open GFF2,"<",$gff2 or die "$!";   
open EXP1,"<",$expmat1 or die "$!";   
open EXP2,"<",$expmat2 or die "$!";   

while (<BED1>) { 
  chomp; 
  my @tt = split /\t+/; 
  $phage1->{$tt[3]}->{chr}=$tt[0];
  $phage1->{$tt[3]}->{begin}=$tt[1];
  $phage1->{$tt[3]}->{end}=$tt[2];
} 

while (<BED2>) { 
  chomp; 
  my @tt = split /\t+/; 
  $phage2->{$tt[3]}->{chr}=$tt[0];
  $phage2->{$tt[3]}->{begin}=$tt[1];
  $phage2->{$tt[3]}->{end}=$tt[2];
}

while (<GFF1>) { 
  chomp; 
  my @tt = split /\t+/;
  $tt[8] =~ m/ID=(.*?);.*/;
  my $id = $1; 
  $tt[8] =~ m/.*;Name=(.*?);.*/;
  my $name = ($1 ne "") ? $1 : "NONE"; 
  $loc1->{$id}->{chr}=$tt[0]; 
  $loc1->{$id}->{begin}=$tt[3]-1;  ## we're doing all 0-based  
  $loc1->{$id}->{end}=$tt[4]; 
  $loc1->{$id}->{name}=$name;
} 

while (<GFF2>) { 
  chomp; 
  my @tt = split /\t+/; 
  $tt[8] =~ m/ID=(.*?);.*/;
  my $id = $1; 
  $tt[8] =~ m/.*;Name=(.*?);.*/;
  my $name = ($1 ne "") ? $1 : "NONE"; 
  $loc2->{$id}->{chr}=$tt[0]; 
  $loc2->{$id}->{begin}=$tt[3]-1;  ## we're doing all 0-based  
  $loc2->{$id}->{end}=$tt[4]; 
  $loc2->{$id}->{name}=$name;
} 

my $h1 = <EXP1>; 
chomp $h1; 
my $h2 = <EXP2>; 
chomp $h2; 
my @h1 = split /\t+/, $h1;
my @h2 = split /\t+/, $h2;
my $nsmp1 = scalar @h1 - 1;        ## how many samples in matrix 1
my $nsmp2 = scalar @h2 - 1;        ## same for exp. matrix 2
$h1 = join "\t",@h1[3..$nsmp1]; 
$h2 = join "\t",@h2[3..$nsmp2]; 

while (<EXP1>) { 
  chomp; 
  my @tt = split /\t+/;
  my $id = shift @tt; 
  my $name = shift @tt;  
  shift @tt; 
  my $exp = join "\t",@tt; 
  $exp1->{$id}->{name}=$name; 
  $exp1->{$id}->{exp}=$exp; 
}

while (<EXP2>) { 
  chomp; 
  my @tt = split /\t+/;
  my $id = shift @tt; 
  my $name = shift @tt;  
  shift @tt; 
  my $exp = join "\t",@tt; 
  $exp2->{$id}->{name}=$name; 
  $exp2->{$id}->{exp}=$exp; 
}
 
## let's do it finally! 
## assume that map has format name-(tab)-id1-(tab)-id2

print "D23_id\tD23_anno\tL474_id\tL474_anno\tGene_name\t$h1\t$h2\n"; 

while (<MAP>) { 
  chomp;
  my $ann1=""; 
  my $ann2=""; 
  my $expvec1=""; 
  my $expvec2=""; 
  my $id1 = (split /\t+/)[0]; 
  my $id2 = (split /\t+/)[1];
  my $names = (split /\t+/)[2]; 
 
  if ($id1 ne "NONE") { 
    my $chr1 = $loc1->{$id1}->{chr}; 
    my $beg1 = $loc1->{$id1}->{begin}; 
    my $end1 = $loc1->{$id1}->{end}; 
    if ($chr1 eq $chr_name1) { 
      $ann1 = "chrom"; 
      foreach my $phid (keys %{$phage1}) { 
	my $pb1 = $phage1->{$phid}->{begin}; 
        my $pe1 = $phage1->{$phid}->{end}; 
        $ann1 = "prophage" if (($beg1 >= $pb1 && $beg1 <= $pe1) || ($end1 >= $pb1 && $end1 <= $pb1));
      } 
    } else { 
      $ann1 = "plasmid"; 
    }
  } else { 
    $ann1 = "NONE"; 
  }
  
  if ($id2 ne "NONE") { 
    my $chr2 = $loc2->{$id2}->{chr}; 
    my $beg2 = $loc2->{$id2}->{begin}; 
    my $end2 = $loc2->{$id2}->{end}; 
    if ($chr2 eq $chr_name2) { 
      $ann2 = "chrom"; 
      foreach my $phid (keys %{$phage2}) { 
	my $pb2 = $phage2->{$phid}->{begin}; 
        my $pe2 = $phage2->{$phid}->{end}; 
        $ann2 = "prophage" if (($beg2 >= $pb2 && $beg2 <= $pe2) || ($end2 >= $pb2 && $end2 <= $pb2));
      } 
    } else { 
      $ann2 = "plasmid"; 
    }
  } else { 
    $ann2 = "NONE"; 
  }
  
  if ($id1 ne "NONE") {
    $expvec1 = $exp1->{$id1}->{exp};
  } else { 
    my $veclen=$nsmp1-2;
    my @zero = ("0") x $veclen; 
    $expvec1 = join "\t",@zero;
  }
  
  if ($id2 ne "NONE") {
    $expvec2 = $exp2->{$id2}->{exp};
  } else {
    my $veclen=$nsmp2-2;
    my @zero = ("0") x $veclen; 
    $expvec2 = join "\t",@zero;
  } 
  
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$id1,$ann1,$id2,$ann2,$names,$expvec1,$expvec2;  
} 
  

close BED1;
close BED2;
close GFF1;
close GFF2;
close EXP1;
close EXP2; 
close MAP; 
