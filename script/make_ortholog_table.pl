#!/usr/bin/env perl 

## this script replaces annotate_CDS.pl and annotate_ncRNA.pl in v0.6 and above
## here we make the final table that includes 
## 1) unique name (disambiguated with _*, e.g. citA/citA_2/citA_3 etc,
## 2) type - CDS or ncRNA 
## 3) location (chromosome/prophage/plasmid),
## 4-n) all IDs (locus tags) for X study + Y reference strains. 

## roary gene_presence_absence.csv file is pre-parsed and converted Dos->Unix using Perl one-liner 

use strict; 
use warnings; 
use Data::Dumper; 

if (scalar @ARGV != 4) {
  print STDERR "Usage: ./make_ortholog_table.pl <full_wdir> <roary_unix_tsv> <bacpipe_config> <modified_ref_fa>\n";
  exit 1
}

my $wdir = shift @ARGV; 
my $roary_unix_csv = shift @ARGV; 
my $config = shift @ARGV;
my $ref_fa = shift @ARGV;  


open FA,"<",$ref_fa or die "$!"; 
open ROARY,"<",$roary_unix_csv or die "$!"; 
open CONFIG,"<",$config or die "$!"; 

my $blast_type = {}; ## this is an abhorrent hack to deal with pseudogene problem in ref strains. Please think of something more elegant ffs. 
my $strains = {};   ## for config parsing  
my $genes = {};     ## strain- and lt-based data rec
my $names = {};     ## key = unique name, output hash 
my @ref_strains; 
my @study_strains; 

while (<FA>) {
  if (m/^>(.*)\.(.*?)$/) { 
    my $name = $1;
    my $type = $2;
    die "ERROR: You can't have blast types other than CDS/ncRNA/misc!\n" if ($type ne "CDS" && $type ne "ncRNA" && $type ne "misc");
    #print STDERR "$name\t$type\n";
    $blast_type->{$name} = $type;  
  } 
} 

## define strains using config file 
while (<CONFIG>) {
  chomp; 
  my @t = split /\t+/; 
  if ($t[0] =~ m/Reference/i && !defined $strains->{$t[1]}) { 
    $strains->{$t[1]} = "ref"; 
    push @ref_strains,$t[1]; 
  } elsif (!defined $strains->{$t[1]}) { 
    $strains->{$t[1]} = "study";
    push @study_strains,$t[1];
  } 
}
## sort strains within each group to make output reproducible   
@ref_strains = sort { $a cmp $b } @ref_strains; 
@study_strains = sort { $a cmp $b } @study_strains;
my @all_strains = @study_strains; 
push @all_strains,@ref_strains; 

foreach my $strain (@study_strains) { 
  my $prop_overlap = join "",$wdir,"/study_strains/",$strain,"/",$strain,".prophage_overlap.tsv";
  my $fasta_index = join "",$wdir,"/study_strains/",$strain,"/",$strain,".genome.fa.fai";
  my $united_gff = join "",$wdir,"/study_strains/",$strain,"/",$strain,".united.gff";
  my $match_tsv = join "",$wdir,"/study_strains/",$strain,"/",$strain,".match.tsv";

  ## define and remember chromosome name for each strain
  my $chr_name = `cat $fasta_index | sort -k2,2nr | head -n 1 | cut -f 1`;
  chomp $chr_name;
  $genes->{$strain}->{chr_name} = $chr_name;
  
  ## parse the GFF 
  open GFF,"<",$united_gff or die "$!"; 
  while (<GFF>) { 
    chomp; 
    my @t = split /\t+/; 
    $t[8] =~ m/ID=(.*?);/; 
    my $lt = $1;
    ## all locus_tags are unique - after processing with unify_study_strain.pl
    ## you can also preprocess GFF files your way but warranty is void in this case :) 
    $genes->{$strain}->{$lt}->{chr} = $t[0]; 
    $genes->{$strain}->{$lt}->{type} = $t[2]; 
    $genes->{$strain}->{$lt}->{beg} = $t[3]; 
    $genes->{$strain}->{$lt}->{end} = $t[4]; 
    $genes->{$strain}->{$lt}->{strand} = $t[6]; 
    if ($t[0] eq $genes->{$strain}->{chr_name}) { 
      $genes->{$strain}->{$lt}->{loc} = "chromosome";
    } else { 
      $genes->{$strain}->{$lt}->{loc} = "plasmid";
    }  
  } 
  close GFF;  
  
  ## use prophage_overlap file obtained with bedtools 
  ## 14 cols - 1-9 is GFF, 10-13 is prophage bed, 14 is overlap  
  open PROP,"<",$prop_overlap or die "$!"; 
  while (<PROP>) { 
    chomp;
    my @t = split /\t+/;
    $t[8] =~ m/ID=(.*?);/; 
    my $lt = $1; 
    $genes->{$strain}->{$lt}->{loc} = "prophage";
  }  
  close PROP; 
 
  open MATCH,"<",$match_tsv or die "$!"; 
  while (<MATCH>) { 
    chomp; 
    my $lt = (split /\t+/)[0]; 
    my $bname = (split /\t+/)[1];  ## blast name 
    $genes->{$strain}->{$lt}->{bname} = $bname;
    my $type = $genes->{$strain}->{$lt}->{type};
    my $loc = $genes->{$strain}->{$lt}->{loc};
 
    if ($type ne "CDS" && $type ne "other" && defined $names->{$bname}->{blast}->{1}->{$strain}) {
      ## if we have seen this name for this strain already, keep adding entries
      ## no paralog entries for CDS - let Roary handle this 
      print STDERR "WARNING: more than 1 defined match for blast name $bname, strain $strain!\n"; 
      my $i = 1; 
      while (defined $names->{$bname}->{blast}->{$i}->{$strain}) { 
        $i++; 
      } 
      $names->{$bname}->{blast}->{$i}->{$strain}->{lt} = $lt;  
      $names->{$bname}->{blast}->{$i}->{$strain}->{type} = $type;  
      $names->{$bname}->{blast}->{$i}->{$strain}->{loc} = $loc; 
    } elsif ($type ne "CDS" && $type ne "other") { 
      $names->{$bname}->{blast}->{1}->{$strain}->{lt} = $lt;  
      $names->{$bname}->{blast}->{1}->{$strain}->{type} = $type;  
      $names->{$bname}->{blast}->{1}->{$strain}->{loc} = $loc; 
    }  
  } 
  close MATCH; 
} 

foreach my $strain (@ref_strains) { 
  ## ref strains don't have defined prophages/intervals 
  my $fasta_index = join "",$wdir,"/ref_strains/",$strain,"/",$strain,".genome.fa.fai";
  my $clean_gff = join "",$wdir,"/ref_strains/",$strain,"/",$strain,".clean.gff";
  my $match_tsv = join "",$wdir,"/ref_strains/",$strain,"/",$strain,".match.tsv";

  ## define and remember chromosome name for each strain
  my $chr_name = `cat $fasta_index | sort -k2,2nr | head -n 1 | cut -f 1`;
  chomp $chr_name;
  $genes->{$strain}->{chr_name} = $chr_name;
  
  ## parse the GFF 
  open GFF,"<",$clean_gff or die "$!"; 
  while (<GFF>) { 
    chomp; 
    my @t = split /\t+/; 
    $t[8] =~ m/ID=(.*?);/; 
    my $lt = $1;
    ## all locus_tags are unique - after processing with unify_study_strain.pl
    ## you can also preprocess GFF files your way but warranty is void in this case :) 
    $genes->{$strain}->{$lt}->{chr} = $t[0]; 
    $genes->{$strain}->{$lt}->{type} = $t[2]; 
    $genes->{$strain}->{$lt}->{beg} = $t[3]; 
    $genes->{$strain}->{$lt}->{end} = $t[4]; 
    $genes->{$strain}->{$lt}->{strand} = $t[6]; 
  } 
  close GFF;  
  
  open MATCH,"<",$match_tsv or die "$!"; 
  while (<MATCH>) { 
    chomp; 
    my $lt = (split /\t+/)[0]; 
    my $bname = (split /\t+/)[1]; 
    $genes->{$strain}->{$lt}->{bname} = $bname;
    my $type = $blast_type->{$bname};  ## this is to take care of pseudogenes/misc  
 
    if ($type ne "CDS" && defined $names->{$bname}->{blast}->{1}->{$strain}) {
      ## if we have seen this name for this strain already, keep adding entries
      ## no location is defined for reference strains 
      print STDERR "WARNING: more than 1 defined match for blast name $bname, strain $strain!\n"; 
      my $i = 1; 
      while (defined $names->{$bname}->{blast}->{$i}->{$strain}) { 
        $i++; 
      } 
      $names->{$bname}->{blast}->{$i}->{$strain}->{lt} = $lt;  
      $names->{$bname}->{blast}->{$i}->{$strain}->{type} = $type;  
    } elsif ($type ne "CDS" && $type ne "other") { 
      $names->{$bname}->{blast}->{1}->{$strain}->{lt} = $lt;  
      $names->{$bname}->{blast}->{1}->{$strain}->{type} = $type;  
    }  
  } 
  close MATCH; 
} 

my $header_line = <ROARY>;
chomp $header_line;  
my @header = split /\t/,$header_line;
my $new_header = "Gene_name\tType\tLocation";  
for (my $i = 0; $i < scalar @header; $i++) {
  my $strain = $header[$i];
  ## if we recognize tag of ref/study strain, we record which column of Roary output is it in 
  if (defined $genes->{$strain}) { 
    $genes->{$strain}->{index} = $i;
    $new_header = join "\t",$new_header,$strain; 
  }
}

## order of strains is as follows: 1) unix-alphabetical study; 2) unix-alphabetical reference
   
OUTER: while (<ROARY>) {
  chomp;
  ## assuming pre-formatted gene presence-absence file with tabs; get all the empty fields right  
  my @t = split /\t/,$_,-1;
  my $rname = ($t[1] eq "") ? $t[0] : $t[1];
  ## everything with underscore is to be stripped to original name
  ## ie tnpA_1a, tnpA and tnpA_2 will all become tnpA, and then get a unique new name (tnpA, tnpA_2, tnpA_3, etc) 
  ## dashes and numbers remain intact, so tnpA2 or RyhB-1 remain the same 
  $rname =~ s/_(.*?)$//g if ($rname !~ m/group_/); 

  foreach my $strain (@all_strains) { 
    my $index = $genes->{$strain}->{index};
    my $lt = ($t[$index] eq "") ? "NONE" : $t[$index]; 

    $genes->{$strain}->{$lt}->{rname} = $rname if ($lt ne "NONE");
    ## we want to skip misc - if there's already blast entry associated with at least 1 lt, skip the whole roary line 
    if (defined $genes->{$strain}->{$lt}->{bname} && $blast_type->{$genes->{$strain}->{$lt}->{bname}} eq "misc") { 
      delete $names->{$rname}->{roary}; 
      next OUTER; 
    } 
 
    if (defined $names->{$rname}->{roary}->{1}->{$strain}) {
      ## if we have seen this name for this strain already, keep adding entries
      ## no location is defined for reference strains 
      my $i = 1; 
      while (defined $names->{$rname}->{roary}->{$i}->{$strain}) {
        $i++; 
      } 
      $names->{$rname}->{roary}->{$i}->{$strain}->{lt} = $lt;  
      $names->{$rname}->{roary}->{$i}->{$strain}->{type} = $genes->{$strain}->{$lt}->{type};  
      $names->{$rname}->{roary}->{$i}->{$strain}->{loc} = $genes->{$strain}->{$lt}->{loc} if (defined $genes->{$strain}->{$lt}->{loc}); 
    } else { 
      $names->{$rname}->{roary}->{1}->{$strain}->{lt} = $lt;  
      $names->{$rname}->{roary}->{1}->{$strain}->{type} = $genes->{$strain}->{$lt}->{type};  
      $names->{$rname}->{roary}->{1}->{$strain}->{loc} = $genes->{$strain}->{$lt}->{loc} if (defined $genes->{$strain}->{$lt}->{loc}); 
    } 
  }  
} 

## you got all orthology records in $names now; it's just down to printing them 

#print STDERR Dumper $names;
print "$new_header\n";    

foreach my $name (keys %{$names}) { 
  if (!defined $names->{$name}->{blast}) { 
    ## no blast hits with this name, just 1 or more roary hits 
    ## case 1-2: roary 1 or 2+, blast 0
    my $i = 1; 
    while (defined $names->{$name}->{roary}->{$i}) {  
      my $type = ""; 
      my $loc = ""; 
      my $out = ""; 
      foreach my $strain (@all_strains) {
        $type = $names->{$name}->{roary}->{$i}->{$strain}->{type} if defined $names->{$name}->{roary}->{$i}->{$strain}->{type}; 
        $loc = $names->{$name}->{roary}->{$i}->{$strain}->{loc} if defined $names->{$name}->{roary}->{$i}->{$strain}->{loc}; 
        my $lt = (defined $names->{$name}->{roary}->{$i}->{$strain}->{lt}) ? $names->{$name}->{roary}->{$i}->{$strain}->{lt} : "NONE";
        $out = join "\t",$out,$lt;   
      }
      my $appended_name = ($i > 1) ? join '_',$name,$i : $name; ## when there's tnpA x2, output tnpA and tnpA_2
      $loc = "chromosome" if ($loc eq ""); 
      $out = join "",$appended_name,"\t",$type,"\t",$loc,$out;
      print "$out\n";
      $i++; 
    }  
  } elsif (!defined $names->{$name}->{roary}) {
    ## no roary hits with this name, just 1 or more blast hits 
    ## case 3-4: roary 0, blast 1 or 2+
    my $i = 1; 
    while (defined $names->{$name}->{blast}->{$i}) {  
      my $type = ""; 
      my $loc = ""; 
      my $out = ""; 
      foreach my $strain (@all_strains) {
        $type = $names->{$name}->{blast}->{$i}->{$strain}->{type} if defined $names->{$name}->{blast}->{$i}->{$strain}->{type}; 
        $loc = $names->{$name}->{blast}->{$i}->{$strain}->{loc} if defined $names->{$name}->{blast}->{$i}->{$strain}->{loc}; 
        my $lt = (defined $names->{$name}->{blast}->{$i}->{$strain}->{lt}) ? $names->{$name}->{blast}->{$i}->{$strain}->{lt} : "NONE";
        $out = join "\t",$out,$lt;   
      }
      my $appended_name = ($i > 1) ? join '_',$name,$i : $name; ## when there's tnpA x2, output tnpA and tnpA_2
      $loc = "chromosome" if ($loc eq ""); 
      $out = join "",$appended_name,"\t",$type,"\t",$loc,$out;
      print "$out\n" if ($type ne "CDS");
      $i++; 
    }  
  } else { 
    ## both roary and blast are defined for this bee yotch  
    ## (sorry I'm really exhausted at this point) 
    ## if blast_type is CDS we ignore blast part of hash completely 
    ## otherwise (if it's misc or ncRNA) we ignore roary 
    ## basically it should only happen for misc
    if ($blast_type->{$name} eq "CDS") { 
      my $i = 1; 
      while (defined $names->{$name}->{roary}->{$i}) {  
        my $type = ""; 
        my $loc = ""; 
        my $out = ""; 
        foreach my $strain (@all_strains) {
          $type = $names->{$name}->{roary}->{$i}->{$strain}->{type} if defined $names->{$name}->{roary}->{$i}->{$strain}->{type}; 
          $loc = $names->{$name}->{roary}->{$i}->{$strain}->{loc} if defined $names->{$name}->{roary}->{$i}->{$strain}->{loc}; 
          my $lt = (defined $names->{$name}->{roary}->{$i}->{$strain}->{lt}) ? $names->{$name}->{roary}->{$i}->{$strain}->{lt} : "NONE";
          $out = join "\t",$out,$lt;   
        }
        my $appended_name = ($i > 1) ? join '_',$name,$i : $name; ## when there's tnpA x2, output tnpA and tnpA_2
        $loc = "chromosome" if ($loc eq ""); 
        $out = join "",$appended_name,"\t",$type,"\t",$loc,$out;
        print "$out\n";
        $i++; 
      }  
    } else {   
      my $i = 1; 
      while (defined $names->{$name}->{blast}->{$i}) {  
        my $type = ""; 
        my $loc = ""; 
        my $out = ""; 
        foreach my $strain (@all_strains) {
          $type = $names->{$name}->{blast}->{$i}->{$strain}->{type} if defined $names->{$name}->{blast}->{$i}->{$strain}->{type}; 
          $loc = $names->{$name}->{blast}->{$i}->{$strain}->{loc} if defined $names->{$name}->{blast}->{$i}->{$strain}->{loc}; 
          my $lt = (defined $names->{$name}->{blast}->{$i}->{$strain}->{lt}) ? $names->{$name}->{blast}->{$i}->{$strain}->{lt} : "NONE";
          $out = join "\t",$out,$lt;   
        }
        my $appended_name = ($i > 1) ? join '_',$name,$i : $name; ## when there's tnpA x2, output tnpA and tnpA_2
        $loc = "chromosome" if ($loc eq ""); 
        $out = join "",$appended_name,"\t",$type,"\t",$loc,$out;
        print "$out\n" if ($type ne "CDS");
        $i++;
      } 
    }  
  }  
}

close ROARY; 
close CONFIG; 
close FA; 
