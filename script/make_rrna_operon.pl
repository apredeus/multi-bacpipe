#!/usr/bin/env perl 

## take united GFF file - united Blastn/Prokka output 
## return rRNA operons and singleton tRNAs as BED intervals
## operon is extended until 50 bp before nearest non-rRNA/non-tRNA feature  

use strict; 
use warnings; 
use Data::Dumper; 

if ($#ARGV != 1) {
  die "USAGE: make_rrna_operon.pl <prokka_gff> <gene_gff>\n";
}

my $prokka_gff = shift @ARGV; 
my $gene_gff = shift @ARGV; 
open PRK_GFF,"<",$prokka_gff or die "$!"; 
open GEN_GFF,"<",$gene_gff or die "$!";  ## "gene" GFF - could be $TAG.gene.gff for --simple, or $TAG.united.gff for --multi. 

my $genes = {}; ## all features in both GFF files 
my @rtrna;      ## rRNA/tRNA as defined by Prokka
my $clen = {};  ## estimate of contigs length (as judged from GFF annotations, since we don't use fasta)  

while (<PRK_GFF>) {
  if (! m/\t/) {
    next; 
  } else { 
    chomp; 
    my @t = split /\t+/; 
    if ($t[8] =~ m/ID=(.*?);/) { 
      my $id = $1; 
      $genes->{$id}->{chr} = $t[0]; 
      $genes->{$id}->{beg} = $t[3]-1; 
      $genes->{$id}->{end} = $t[4];
      
      if (defined $clen->{$t[0]}) { 
        $clen->{$t[0]} = $t[4] if ($t[4] > $clen->{$t[0]}); 
      } else { 
        $clen->{$t[0]} = $t[4];
      }
  
      if ($t[2] eq "tRNA") {  
        $genes->{$id}->{type} = "tRNA"; 
        push @rtrna,$id; 
      } elsif ($t[2] eq "rRNA") {  
        $genes->{$id}->{type} = "rRNA"; 
        push @rtrna,$id; 
      } else { 
        $genes->{$id}->{type} = "other"; 
      }
    } 
  } 
}

while (<GEN_GFF>) {
  if (! m/\t/) {
    next;
  } else {
    chomp;
    my @t = split /\t+/;
    if ($t[8] =~ m/ID=(.*?);/) {
      ## in case of Prokka/unified GFF collusions happen that made the script produce bugs  
      my $id = "GENE_".$1;
      $genes->{$id}->{chr} = $t[0];
      $genes->{$id}->{beg} = $t[3] - 1;
      $genes->{$id}->{end} = $t[4];

      if (defined $clen->{$t[0]}) { 
        $clen->{$t[0]} = $t[4] if ($t[4] > $clen->{$t[0]}); 
      } else { 
        $clen->{$t[0]} = $t[4];
      }
  
      if ($t[2] eq "tRNA" || $t[8] =~ m/gene_biotype=tRNA/) {
        $genes->{$id}->{type} = "tRNA";
      } elsif ($t[2] eq "rRNA" || $t[8] =~ m/gene_biotype=rRNA/) {
        $genes->{$id}->{type} = "rRNA";
      } else {
        $genes->{$id}->{type} = "other";
      }
    }
  }
}

foreach my $rtrna (@rtrna) { 
  ## distance from nearest non-rRNA, non-tRNA gene to this rRNA/tRNA gene from the left/right
  ## upd: also, limit the left border with 0 and right with the contig length 
  my $d_left  = 1000000; 
  my $d_right = 1000000;
  my $rbeg = $genes->{$rtrna}->{beg};
  my $rend = $genes->{$rtrna}->{end};
  my $rchr = $genes->{$rtrna}->{chr};
 
  foreach my $key (keys %{$genes}) { 
    my $cbeg = $genes->{$key}->{beg}; 
    my $cend = $genes->{$key}->{end}; 
    if ($genes->{$key}->{chr} eq $rchr && $genes->{$key}->{type} ne "rRNA" && $genes->{$key}->{type} ne "tRNA" && $cend-$cbeg <= 50000) { 
      my $max_beg = ($rbeg > $cbeg) ? $rbeg : $cbeg; 
      my $min_end = ($rend < $cend) ? $rend : $cend; 
      if ( $min_end >= $max_beg ) {
      ## this means non-rtRNA feature overlaps rRNA/tRNA
      ## this should never happen, but with blasted references it sometimes does (esp with ncRNAs)
        $d_left = 0; 
        $d_right = 0; 
        last; 
      } else { 
        $d_left  = $rbeg-$cend if ($cend < $rbeg && $rbeg-$cend < $d_left);   ## find closest non-rtRNA feature on the left
        $d_right = $cbeg-$rend if ($cbeg > $rend && $cbeg-$rend < $d_right);  ## find closest non-rtRNA feature on the right 
      }
    } 
  }
  ## if less than 50 bp away, keep the exact borders of rtRNA; if more, extend to -50 bp to next non-rtRNA feature
  my $lbound = ($d_left > 50) ? $genes->{$rtrna}->{beg} - $d_left + 50 : $genes->{$rtrna}->{beg}; 
  my $rbound = ($d_right > 50) ? $genes->{$rtrna}->{end} + $d_right - 50 : $genes->{$rtrna}->{end};
  $lbound = 0 if ($lbound < 0); 
  $rbound = $clen->{$rchr} if ($rbound > $clen->{$rchr} && defined $clen->{$rchr});  
  printf "%s\t%d\t%d\trtRNA.%s\n",$genes->{$rtrna}->{chr},$lbound,$rbound,$rtrna;
}

close PRK_GFF; 
close GEN_GFF; 
