#!/usr/bin/env perl 

## annotate ncRNA using name for orthology

use strict; 
use warnings;
use Data::Dumper;  

my $refdir = shift @ARGV; 
my $config = shift @ARGV; 

open CONFIG,"<",$config or die "$!"; 

my $nc_loc = {}; 
my $nc_coord = {};
my $name_count = {};  
my $max_count = {}; 
my @ref_strains; 
my @study_strains; 

while (<CONFIG>) { 
  ## define prophage ranges for every strain
  ## also, define chromosome name 
  chomp; 
  my $strain = (split /\t+/)[1];
  if (! defined $nc_loc->{$strain} && $_ !~ m/^Reference\t/) { 
    my $prophage = join "",$refdir,"/",$strain,"/",$strain,".prophage.bed";
    my $fasta_index = join "",$refdir,"/",$strain,"/",$strain,".genome.fa.fai";
    my $gff = join "",$refdir,"/",$strain,"/",$strain,".ncRNA.gff";
    my $chr_name = `cat $fasta_index | sort -k2,2nr | head -n 1 | cut -f 1`;
    chomp $chr_name;
    $nc_loc->{$strain}->{chr_name}=$chr_name; 
    push @study_strains,$strain; 

    open PROP,"<",$prophage; 
    while (<PROP>) { 
      chomp;
      my @tt = split /\t+/;
      $nc_loc->{$strain}->{$tt[3]}->{chr}=$tt[0];
      $nc_loc->{$strain}->{$tt[3]}->{beg}=$tt[1];
      $nc_loc->{$strain}->{$tt[3]}->{end}=$tt[2];
    }  
    close PROP; 

    open GFF,"<",$gff; 
    while (<GFF>) { 
      chomp; 
      my @tt = split /\t+/; 
      $tt[8] =~ m/ID=(.*?);.*Name=(.*?);/; 
      my $locus_tag = $1;
      my $name = $2; 
      ## populate data structure based on unique locus tags 
      $nc_coord->{$locus_tag}->{name}=$name; 
      $nc_coord->{$locus_tag}->{chr}=$tt[0]; 
      $nc_coord->{$locus_tag}->{beg}=$tt[3]; 
      $nc_coord->{$locus_tag}->{end}=$tt[4];

      ## name_count keeps track of ncRNA counts and locus tags in each genome 
      push @{$name_count->{$strain}->{$name}->{lt}},$locus_tag; 
      if (! defined $name_count->{$strain}->{$name}->{count}) { 
        $name_count->{$strain}->{$name}->{count}=1;
      } else { 
        $name_count->{$strain}->{$name}->{count}++;
      }

      ## get max count for each name across all strains
      my $cnt = $name_count->{$strain}->{$name}->{count}; 
      if (! defined $max_count->{$name}) {
        $max_count->{$name} = $cnt; 
      } else { 
        $max_count->{$name} = ($max_count->{$name} < $cnt) ? $cnt : $max_count->{$name}; 
      }  
    } 
    close GFF;  
  } elsif (! defined $nc_loc->{$strain} && $_ =~ m/^Reference\t/) {
      $nc_loc->{$strain}->{type}="reference"; 
      push @ref_strains,$strain;
      ## TODO: add parsing of reference strains (?)   
  } 
}

## both have to be sorted alphabetically to make strain order reproducible 
my @strains = sort { $a cmp $b } @study_strains;
@ref_strains = sort { $a cmp $b } @ref_strains; 
push @strains,@ref_strains;

my $new_header = "ncRNA_name\tncRNA_loc";  

foreach my $strain (@strains) { 
  $new_header = join "\t",$new_header,$strain; 
} 
print "$new_header\n"; 

foreach my $name (keys %{$max_count}) { 
  for (my $i=0; $i < $max_count->{$name}; $i++) {
    my $output = ""; 
    my $gene_loc = "chromosome";                    ## default 
    foreach my $strain (@strains) {
      my $locus_tag; 
      ## check yoself before you wreck yoself
      if (defined $name_count->{$strain}->{$name}->{lt}) {
        my @lt = @{$name_count->{$strain}->{$name}->{lt}}; 
        if (defined $lt[$i]) { 
          $locus_tag = $lt[$i]; 
        } else { 
          $locus_tag = "NONE"; 
        }  
      } else { 
        $locus_tag = "NONE";          ## this should also account for ref. strains 
      } 
      
      if ($locus_tag ne "NONE") {
        my $chr = $nc_coord->{$locus_tag}->{chr};
        my $beg = $nc_coord->{$locus_tag}->{beg};
        my $end = $nc_coord->{$locus_tag}->{end};
        
        ## print STDERR "LT is $locus_tag\n" if (! defined $chr);
        $gene_loc = "plasmid" if ($chr ne $nc_loc->{$strain}->{chr_name});
        foreach my $key (keys %{$nc_loc->{$strain}}) {
          if ($key ne "chr_name" && $nc_loc->{$strain}->{$key}->{chr} eq $chr) {
            my $prop_beg = $nc_loc->{$strain}->{$key}->{beg};
            my $prop_end = $nc_loc->{$strain}->{$key}->{end};
            $gene_loc = "prophage" if (($prop_beg <= $beg && $prop_end >= $beg) || ($prop_beg <= $end && $prop_end >= $end));
          }
        }
      }
      $output = join "\t",$output,$locus_tag;
    }
    print "$name\t$gene_loc$output\n"; 
  }
}

close CONFIG;  
