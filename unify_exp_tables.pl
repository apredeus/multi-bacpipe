#!/usr/bin/env perl 

use strict; 
use warnings; 
use Data::Dumper; 

## take annotated tables of CDS and ncRNA 
## and expression tables of all individual strains 
## and get the master expression table 
## gene names that have a non - group_* name with underscore (e.g. argC_1) 
## are stripped or the _* part 
## non-unique genes are then assigned new unique numbers (e.g. 2 x argC becomes argC_1, argC_2)
## group_* genes are replaced with STM* whenever possible. 


my $wdir = shift @ARGV; 
my $config = shift @ARGV; 
my $ann_cds = shift @ARGV; 
my $ann_nc = shift @ARGV;
my $ref_strain = shift @ARGV;  

my $sample_counts = {};
my $strain_idx = {};  
my $ann = {}; 
my $counts = {};
my $tpms = {}; 
my @ref_strains;  

open CONFIG,"<",$config or die "$!"; 
open ANN_CDS,"<",$ann_cds or die "$!"; 
open ANN_NC,"<",$ann_nc or die "$!"; 

## parse config, get list of study strains
## get how many samples are there for each strain 
while (<CONFIG>) { 
  chomp;
  my @tt = split /\t+/;  
  my $strain = $tt[1]; 
  if (! m/^Reference/ && defined $sample_counts->{$strain}) {
    $sample_counts->{$strain}++; 
  } elsif (! m/^Reference/ && ! defined $sample_counts->{$strain}) { 
    $sample_counts->{$strain}=1; 
  } elsif (m/^Reference/) {
    ## we assume that each ref strain is only listed once 
    push @ref_strains,$tt[1];
  } 
}

my @study_strains = keys %{$sample_counts}; 
@study_strains = sort { $a cmp $b } @study_strains;
@ref_strains   = sort { $a cmp $b } @ref_strains;
my @strains = @study_strains; 
push @strains,@ref_strains; 

my $ref_index; 
#print Dumper($sample_counts); 

## parse annotation files, populate hash
my $cds_header = <ANN_CDS>; 
chomp $cds_header; 
my $nc_header = <ANN_NC>; 
chomp $nc_header; 
my @str_names = split /\t+/,$cds_header;

$cds_header =~ s/Gene_name\tGene_loc\t//g;  
$nc_header  =~ s/ncRNA_name\tncRNA_loc\t//g;  
die "ERROR: CDS and ncRNA annotations have different strain order!\n" if ($cds_header ne $nc_header); 

## get the column numbers 
for (my $i=0; $i < scalar @str_names; $i++) { 
  foreach my $strain (@strains) {
    $strain_idx->{$strain} = $i-2 if ($str_names[$i] eq $strain);
  }
} 

#print Dumper($strain_idx); 

## populate annotation hash 
while (<ANN_CDS>) { 
  chomp;
  my $name_count = 2;  
  my @tt   = split /\t+/;
  my $name = shift @tt; 
  my $loc  = shift @tt;
  my $ref_idx = $strain_idx->{$ref_strain};
  #print "DEBUG1: $ref_strain $ref_idx\n";  
  my $ref_lt = $tt[$ref_idx]; 
  
  $name =~ s/_*//g if ($name !~ m/^group_/);
  $name = $ref_lt if ($name =~ m/^group_/ && $ref_lt ne "NONE"); 
  
  if (! defined $ann->{$name}) { 
    foreach my $strain (@strains) { 
      my $idx = $strain_idx->{$strain}; 
      my $locus_tag = $tt[$idx];                    ## account for two shifts above  
      $ann->{$name}->{$strain}->{lt} = $locus_tag;
      $ann->{$name}->{loc} = $loc;                    ## chr/prophage/plasmid  
    } 
  } else { 
    $name = join "_",$name,$name_count; 
    $name_count++; 
    foreach my $strain (@strains) {
      my $idx = $strain_idx->{$strain};
      my $locus_tag = $tt[$idx];                    ## account for two shifts above  
      $ann->{$name}->{$strain}->{lt} = $locus_tag;
      $ann->{$name}->{loc} = $loc;                    ## chr/prophage/plasmid  
    }
  }
}

while (<ANN_NC>) { 
  chomp;
  my $name_count = 2;  
  my @tt   = split /\t+/;
  my $name = shift @tt; 
  my $loc  = shift @tt;
  $name    =~ s/_*//g if ($name !~ m/^group_/); 
  
  if (! defined $ann->{$name}) { 
    foreach my $strain (@strains) { 
      my $idx = $strain_idx->{$strain}; 
      my $locus_tag = $tt[$idx];                    ## account for two shifts above  
      $ann->{$name}->{$strain}->{lt} = $locus_tag;
      $ann->{$name}->{loc} = $loc;                  ## chr/prophage/plasmid  
    } 
  } else { 
    $name = join "_",$name,$name_count; 
    $name_count++; 
    foreach my $strain (@strains) {
      my $idx = $strain_idx->{$strain};
      my $locus_tag = $tt[$idx];                    ## account for two shifts above  
      $ann->{$name}->{$strain}->{lt} = $locus_tag;
      $ann->{$name}->{loc} = $loc;                  ## chr/prophage/plasmid  
    }
  }
}

#print Dumper($ann); 

## finally, let's print the full annotated table
## for this, we first populate the counts hash
foreach my $strain (@study_strains) { 
  my $counts_file = join "/",$wdir,"exp_tables",$strain; 
  my $tpm_file = join "/",$wdir,"exp_tables",$strain; 
  $counts_file = join ".",$counts_file,"counts.tsv";
  $tpm_file = join ".",$tpm_file,"TPM.tsv";
  open COUNTS,"<",$counts_file or die "$!";  
  open TPMS,"<",$tpm_file or die "$!";  
  
  my $samples = <TPMS>; 
  $samples = <COUNTS>;
  $samples =~ s/$strain\t//g; 
  $counts->{$strain}->{samples}=$samples; 
  $tpms->{$strain}->{samples}=$samples; 
  
  while (<COUNTS>) { 
    chomp; 
    m/^(.*?)\t(.*)/;
    my $locus_tag = $1;
    my $count_str = $2; 
    $counts->{$strain}->{$locus_tag}=$count_str; 
  } 
  while (<TPMS>) { 
    chomp; 
    m/^(.*?)\t(.*)/;
    my $locus_tag = $1;
    my $tpm_str = $2; 
    $tpms->{$strain}->{$locus_tag}=$tpm_str; 
  }
  
  close COUNTS; 
  close TPMS; 
}

#print Dumper($tpms);   

## and now, the actual table. Serously, fuck this thing sideways. 
foreach my $name (keys %{$ann}) { 
  my $flag = 1; 
  my $output = join "\t",$name,$ann->{$name}->{loc};
  foreach my $strain (@strains) { 
    $output = join "\t",$output,$ann->{$name}->{$strain}->{lt};
  }
  foreach my $strain (@study_strains) { 
    my $locus_tag = $ann->{$name}->{$strain}->{lt}; 
    my $counts = $counts->{$strain}->{$locus_tag}; 
    $output = join "\t",$output,$counts; 
    #my $tpms = $tpms->{$strain}->{$locus_tag}; 
  } 
  print "$output\n"; 
}

close CONFIG; 
close ANN_CDS; 
close ANN_NC; 
