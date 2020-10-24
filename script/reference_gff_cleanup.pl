#!/usr/bin/env perl 

## similar to simple_gff_cleanup.pl, but with some differences
## feature types are now stated in 3rd column  
## basically all "gene" features are included, and "gene" borders are REDEFINED to CDS borders if type is CDS. 
## features that have the same locus_tag (usually halves of a pseudogene) are SPLIT: shorter one is dropped, longer one is kept
## "gene" entries with no locus_tag are simply skipped (no locus tag, no cartoons)  
## "ncRNA" are all things that match *rna (non-case-spec), but not rRNA/tRNA. 

use strict;
use warnings; 

if (scalar @ARGV != 1) {
  die "Usage: ./reference_gff_cleanup.pl <ncbi_gff>\n";
}

my $gff = shift @ARGV; 

open GFF,"<",$gff or die "$!"; 

my $genes = {};
my $prod = {};
my ($in_gene,$in_pseudo) = (0)x2; 
my ($in_cds,$in_ncrna,$in_trna,$in_rrna,$in_other) = (0)x5; 
my ($out_cds,$out_pseudo,$out_ncrna,$out_trna,$out_rrna,$out_other) = (0)x6; 

while (<GFF>) { 
  if (m/\tgene\t/ || m/\tpseudogene\t/) {
    ## parse Level 1 features 
    my @t = split /\t+/; 
    my $id = ($t[8] =~ m/ID=(.*?);/i) ? $1 : "NONE"; 
    my $name = ($t[8] =~ m/Name=(.*?);/i) ? $1 : "NONE"; 
    my $biotype = ($t[8] =~ m/gene_biotype=(\w+)/i) ? $1 : "NONE";
    my $lt   = ($t[8] =~ m/;locus_tag=([-\w]+)/i) ? $1 : "NONE"; 
    
    if ($t[8] =~ m/pseudo=true/i || $t[2] eq "pseudogene") { 
      $biotype = "pseudogene";
      $in_pseudo++; 
    } else { 
      $in_gene++; 
    }  
    ## this accounts for duplicate locus tags in case of pseudogenes 
    if (! defined $genes->{$lt} && $lt ne "NONE") { 
      $genes->{$lt}->{lt} = $lt;  
      $genes->{$lt}->{id} = $id;  
      $genes->{$lt}->{name} = $name; 
      $biotype = "noncoding_rna" if ($biotype =~ m/rna/i && $biotype ne "rRNA" && $biotype ne "tRNA"); 
      if ($biotype ne "protein_coding" && $biotype ne "pseudogene" && $biotype ne "rRNA" && $biotype ne "tRNA" && $biotype ne "other" && $biotype ne "noncoding_rna" && $biotype ne "NONE") { 
        print STDERR "WARNING: Found unsupported biotype $biotype for locus tag $lt; changed to \"other\".\n"; 
        $biotype = "other"; 
      }  
      $genes->{$lt}->{biotype} = $biotype;
      $genes->{$lt}->{chr} = $t[0];   
      $genes->{$lt}->{beg} = $t[3];   
      $genes->{$lt}->{end} = $t[4];   
      $genes->{$lt}->{strand} = $t[6];
    } elsif ($lt ne "NONE") {
      ## get the coords for the longer feature 
      my $length_old = $genes->{$lt}->{end} - $genes->{$lt}->{beg}; 
      my $length_new = $t[4] - $t[3];
      if ($length_new > $length_old) { 
        $genes->{$lt}->{beg} = $t[3];  
        $genes->{$lt}->{end} = $t[4]; 
      }  
      ## everything else is defined and stays the same 
    } 
  } elsif (m/Parent=gene/i) {
    ## parse Level 2 features 
    my @t = split /\t+/; 
    my $id = ($t[8] =~ m/Parent=(.*?);/i) ? $1 : "NONE"; 
    my $product = ($t[8] =~ m/product=(.*?);/i || $t[8] =~ m/product=(.*)$/i) ? $1 : "NONE"; 
    my $note = ($t[8] =~ m/note=(.*?);/i || $t[8] =~ m/note=(.*)$/i) ? $1 : "NONE";
     
    ## gene0 can be found on chr and plasmid sometimes, so ref via chr
    ## coordinates for CDS to replace gene, in case gene has UTR defined with it (see LT2)  
    if (! defined $prod->{$t[0]}->{$id}) { 
      $prod->{$t[0]}->{$id}->{product} = $product; 
      $prod->{$t[0]}->{$id}->{note} = $note; 
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
    ## let's keep track of all things we're reading in 
    if ($t[2] eq "CDS") {
      $in_cds++; 
      $prod->{$t[0]}->{$id}->{type} = "CDS"; 
    } elsif ($t[2] eq "tRNA") {
      $in_trna++;  
      $prod->{$t[0]}->{$id}->{type} = "tRNA"; 
    } elsif ($t[2] eq "rRNA") { 
      $in_rrna++;  
      $prod->{$t[0]}->{$id}->{type} = "rRNA"; 
    } elsif ($t[2] =~ m/rna/i || $t[2] =~ m/transcript/i) { 
      $in_ncrna++;
      $prod->{$t[0]}->{$id}->{type} = "ncRNA"; 
    } else {
      $in_other++; 
      $prod->{$t[0]}->{$id}->{type} = "other"; 
      print STDERR "WARNING: found level 2 feature with unsupported type: $t[2], gene id $id\n; changing feature type to \"other\"."; 
    }  
  }
} 

printf STDERR "Input GFF %s, level 1: %d genes, %d pseudogenes, %d overall.\n",$gff,$in_gene,$in_pseudo,$in_gene+$in_pseudo; 
printf STDERR "Input GFF %s, level 2: %d CDS, %d ncRNA, %d tRNA, %d rRNA, %d others.\n",$gff,$in_cds,$in_ncrna,$in_trna,$in_rrna,$in_other; 

my @keys = sort { $genes->{$a}->{chr} cmp $genes->{$b}->{chr} || $genes->{$a}->{beg} <=> $genes->{$b}->{beg} } keys %{$genes};

printf STDERR "Number of unique locus tags: %d\n",scalar @keys; 

foreach my $lt (@keys) {
  if ($lt ne "NONE") {
    my $chr = $genes->{$lt}->{chr}; 
    my $id = $genes->{$lt}->{id}; 
    my $strand = $genes->{$lt}->{strand}; 
    my $name = $genes->{$lt}->{name};
    my $biotype = $genes->{$lt}->{biotype};
    my $beg = $genes->{$lt}->{beg}; 
    my $end = $genes->{$lt}->{end}; 
    my $type = "NONE";
    my $product = "NONE"; 
    my $note = "NONE"; 

    ## if we have a child (level2) defined for this $lt, let's process the values
    if (defined $prod->{$chr}->{$id}) {  
      $type = $prod->{$chr}->{$id}->{type};  
      $beg = $prod->{$chr}->{$id}->{beg} if ($type eq "CDS");
      $end = $prod->{$chr}->{$id}->{end} if ($type eq "CDS");
      $product = $prod->{$chr}->{$id}->{product}; 
      $note = $prod->{$chr}->{$id}->{note}; 
    } 
 
    ## now let's sort biotype and type; we allow the following types for each
    ## biotype: protein_coding, pseudogene, tRNA, rRNA, ncRNA, other
    ## type: CDS, tRNA, rRNA, ncRNA, other?
    if ($type eq "NONE" && $biotype eq "NONE") {
      print STDERR "DEBUG: feature $lt does not have a defined type OR biotype!\n"; 
    } elsif ($type eq "NONE") {
      $type = $biotype; 
      $type = "CDS" if ($biotype eq "protein_coding"); 
      $type = "CDS" if ($biotype eq "pseudogene"); 
      $type = "ncRNA" if ($biotype eq "noncoding_rna"); 
    } elsif ($biotype eq "NONE") { 
      $biotype = $type; 
      $biotype = "protein_coding" if ($type eq "CDS"); 
      $biotype = "noncoding_rna" if ($type eq "ncRNA"); 
    } 
     
    my $out = sprintf "%s\tBacpipe\t%s\t%d\t%d\t.\t%s\t.\t",$chr,$type,$beg,$end,$strand;
    $out = join '',$out,"ID=",$lt,";";
    ## we skip Name= if it's not defined - no Name=locus tag here. That's for Roary.  
    $out = join '',$out,"Name=",$name,";" if ($name ne "NONE" && $name ne $lt);
 
    if ($biotype eq "protein_coding") { 
      $out_cds++; 
    } elsif ($biotype eq "pseudogene") { 
      $out_pseudo++;
    } elsif ($biotype eq "noncoding_rna") { 
      $out_ncrna++; 
    } elsif ($biotype eq "rRNA") { 
      $out_rrna++; 
    } elsif ($biotype eq "tRNA") { 
      $out_trna++; 
    } else { 
      $out_other++; 
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

printf STDERR "Output GFF: %d CDS, %d pseudogenes (as CDS), %d tRNA, %d rRNA, %d ncRNA, and %d others.\n",$out_cds,$out_pseudo,$out_trna,$out_rrna,$out_ncrna,$out_other;
printf STDERR "Total number of features in the output: %d.\n",$out_cds+$out_pseudo+$out_trna+$out_rrna+$out_ncrna+$out_other; 

close GFF; 
