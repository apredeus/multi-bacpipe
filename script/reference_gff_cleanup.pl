#!/usr/bin/env perl 

## similar to simple_gff_cleanup.pl, but with some differences
## feature types are now stated in 3rd column  
## basically all "gene" features are included, and "gene" borders are REDEFINED to CDS borders if type is CDS. 
## features that have the same locus_tag (usually halves of a pseudogene) are SPLIT: shorter one is dropped, longer one is kept
## "gene" entries with no locus_tag are simply skipped (no locus tag, no cartoons)  
## "ncRNA" are all things that match *rna (non-case-spec), but not rRNA/tRNA. 

use strict;
use warnings; 

if ($#ARGV != 0) {
  die "USAGE: reference_gff_cleanup.pl <ncbi_gff>\n";
}

my $gff = shift @ARGV; 

open GFF,"<",$gff or die "$!"; 

my $genes = {};
my $prod = {};
my ($coding,$pseudo,$ncrna,$trna,$rrna,$other) = (0)x6; 


while (<GFF>) { 
  if (m/\tgene\t/ || m/\tpseudogene\t/) {
    my @t = split /\t+/; 
    my $id = ($t[8] =~ m/ID=(.*?);/) ? $1 : "NONE"; 
    my $name = ($t[8] =~ m/Name=(.*?);/) ? $1 : "NONE"; 
    my $biotype = ($t[8] =~ m/;gene_biotype=(\w+)/) ? $1 : "NONE";
    my $lt   = ($t[8] =~ m/;locus_tag=(\w+)/) ? $1 : "NONE"; 
    ## a bit of extra safety 
    $biotype = "pseudogene" if ($t[8] =~ m/pseudo=true/); 
    ## this accounts for duplicate locus tags in 
    if (! defined $genes->{$lt}) { 
      $genes->{$lt}->{lt} = $lt;  
      $genes->{$lt}->{id} = $id;  
      $genes->{$lt}->{name} = $name;  
      $genes->{$lt}->{biotype} = $biotype;
      $genes->{$lt}->{chr} = $t[0];   
      $genes->{$lt}->{beg} = $t[3];   
      $genes->{$lt}->{end} = $t[4];   
      $genes->{$lt}->{strand} = $t[6];
    } else {
      ## get the coords for the longer feature 
      my $length_old = $genes->{$lt}->{end} - $genes->{$lt}->{beg}; 
      my $length_new = $t[4] - $t[3];
      if ($length_new > $length_old) { 
        $genes->{$lt}->{beg} = $t[3];  
        $genes->{$lt}->{end} = $t[4]; 
      }  
      ## everything else is defined and stays the same 
    } 
  } elsif (m/Parent=gene/) {
    ## products are either in product or Note field; use former, if not there, use latter
    my @t = split /\t+/; 
    my $id = ($t[8] =~ m/Parent=(.*?);/) ? $1 : "NONE"; 
    my $product = ($t[8] =~ m/;product=(.*?);/) ? $1 : "NONE"; 
    my $note = ($t[8] =~ m/;Note=(.*?);/) ? $1 : "NONE";
     
    ## gene0 can be found on chr and plasmid sometimes, so ref via chr
    ## coordinates for CDS to replace gene, in case gene has UTR defined with it (see LT2)  
    if (! defined $prod->{$t[0]}->{$id}) { 
      $prod->{$t[0]}->{$id}->{product} = $product; 
      $prod->{$t[0]}->{$id}->{note} = $note; 
      $prod->{$t[0]}->{$id}->{type} = $t[2]; 
      if ($t[2] =~ m/rna/i && $t[2] ne "rRNA" && $t[2] ne "tRNA") { 
        $prod->{$t[0]}->{$id}->{type} = "ncRNA"; 
      }
      $prod->{$t[0]}->{$id}->{beg} = $t[3]; 
      $prod->{$t[0]}->{$id}->{end} = $t[4];
    } else { 
      my $length_old = $prod->{$t[0]}->{$id}->{end} - $prod->{$t[0]}->{$id}->{beg};
      my $length_new = $t[4] - $t[3];
      if ($length_new > $length_old) { 
        $prod->{$t[0]}->{$id}->{beg} = $t[3];  
        $prod->{$t[0]}->{$id}->{end} = $t[4]; 
      }  
    } 
  } 
} 

my @keys = sort { $genes->{$a}->{chr} cmp $genes->{$b}->{chr} || $genes->{$a}->{beg} <=> $genes->{$b}->{beg} } keys %{$genes};

foreach my $lt (@keys) {
  if ($lt ne "NONE") {
    my $chr = $genes->{$lt}->{chr}; 
    my $id = $genes->{$lt}->{id}; 
    my $type = $prod->{$chr}->{$id}->{type};
    my $beg = ($type eq "CDS") ? $prod->{$chr}->{$id}->{beg} : $genes->{$lt}->{beg}; 
    my $end = ($type eq "CDS") ? $prod->{$chr}->{$id}->{end} : $genes->{$lt}->{end};
    my $strand = $genes->{$lt}->{strand}; 
     
    my $out = sprintf "%s\tBacpipe\t%s\t%d\t%d\t.\t%s\t.\t",$chr,$type,$beg,$end,$strand;
    $out = join ('',$out,"ID=",$lt,";");
    ## we skip Name= if it's not defined - no Name=locus tag here. That's for Roary.  
    my $name = $genes->{$lt}->{name}; 
    $out = join ('',$out,"Name=",$name,";") if ($name ne "NONE");
    my $product = (defined $prod->{$chr}->{$id}->{product}) ? $prod->{$chr}->{$id}->{product} : "NONE"; 
    my $note = (defined $prod->{$chr}->{$id}->{note}) ? $prod->{$chr}->{$id}->{note} : "NONE";
 
    ## if product is defined, then use product; if note is defined instead of product, use note;
    ## if neither is defined, use nothing. 
    my $biotype = $genes->{$lt}->{biotype}; 
    if ($biotype =~ m/rna/i && $biotype ne "rRNA" && $biotype ne "tRNA") { 
      $biotype = "noncoding_rna"; 
    }
    if ($biotype eq "protein_coding") { 
      $coding++; 
    } elsif ($biotype eq "pseudogene") { 
      $pseudo++;
    } elsif ($biotype eq "noncoding_rna") { 
      $ncrna++; 
    } elsif ($biotype eq "rRNA") { 
      $rrna++; 
    } elsif ($biotype eq "tRNA") { 
      $trna++; 
    } else { 
      $other++; 
    } 
    
    $out = join ('',$out,"gene_biotype=",$biotype,";");
   
    if ($product ne "NONE") { 
      $out = join ('',$out,"product=",$product,";");
    } elsif ($note ne "NONE") { 
      $out = join ('',$out,"product=",$note,";");
    } 
    $out = join ('',$out,"locus_tag=",$lt,"\n");
    ## print out every gene entry, even rRNA and tRNA (they won't have expression since reads are removed from BAM)
    ## merge pseudogene halves into one "gene" entry with the same locus tag if LT of halves are the same 
    ## above treatment of pseudogenes ONLY APPLIES WHEN --SIMPLE IS IN USE! 
    print $out; 
  }
} 

print STDERR "GFF output stats: $coding protein coding, $pseudo pseudogenes, $ncrna noncoding RNAs, $trna tRNAs, $rrna rRNAs, $other others.\n"; 

close GFF; 
