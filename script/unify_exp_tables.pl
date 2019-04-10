#!/usr/bin/env perl 

## TODO: make sure this works correctly without the reference strain

use strict; 
use warnings; 
#use Data::Dumper; 

## take annotated tables of CDS and ncRNA 
## and expression tables of all individual strains (in exp_tables)  
## and get the master expression table in counts and TPMs 
## gene names that have a non - group_* name with underscore (e.g. argC_1) 
## are stripped or the _* part 
## non-unique genes are then assigned new unique numbers (e.g. 2 x argC becomes argC_1, argC_2)
## group_* genes are replaced with reference strain locus tags (eg STM for Salmonella) whenever possible. 


my $wdir = shift @ARGV; 
my $config = shift @ARGV; 
my $ortho = shift @ARGV; 
my $refstr = shift @ARGV;  

my $sample_counts = {};
my $strain_idx = {};  
my $ann = {}; 
my $counts = {};
my $tpms = {}; 
my @ref_strains;  

open CONFIG,"<",$config or die "$!"; 
open ORTHO,"<",$ortho or die "$!"; 

## parse config, get list of study strains
## get how many samples are there for each strain 
while (<CONFIG>) { 
  chomp;
  my @t = split /\t+/;  
  my $strain = $t[1]; 
  if (! m/^Reference/ && defined $sample_counts->{$strain}) {
    $sample_counts->{$strain}++; 
  } elsif (! m/^Reference/ && ! defined $sample_counts->{$strain}) { 
    $sample_counts->{$strain}=1; 
  } elsif (m/^Reference/) {
    ## we assume that each ref strain is only listed once 
    push @ref_strains,$t[1];
  } 
}

my @study_strains = keys %{$sample_counts}; 
@study_strains = sort { $a cmp $b } @study_strains;
@ref_strains   = sort { $a cmp $b } @ref_strains;
my @all_strains = @study_strains; 
push @all_strains,@ref_strains; 

my $ref_index; 

## parse annotation files, populate hash
my $ortho_header = <ORTHO>; 
chomp $ortho_header; 
my @str_names = split /\t+/,$ortho_header;

$ortho_header =~ s/Gene_name\tType\tLocation\t//g;  

## get the column numbers 
for (my $i=0; $i < scalar @str_names; $i++) { 
  foreach my $strain (@all_strains) {
    ## -2 accounts for two shifts, see below 
    $strain_idx->{$strain} = $i-3 if ($str_names[$i] eq $strain);
  }
} 

#print Dumper($strain_idx); 

## populate annotation hash 
while (<ORTHO>) { 
  chomp;
  my @t   = split /\t+/;
  my $name = shift @t; 
  my $type = shift @t; 
  my $loc  = shift @t;
  
  my $present_in_study = 0; 
  for (my $i = 0; $i < scalar @study_strains; $i++) { 
    $present_in_study = 1 if ($t[$i] ne "NONE"); 
  } 
  next if (! $present_in_study); 

  my $ref_idx = $strain_idx->{$refstr};
  my $ref_lt = $t[$ref_idx]; 
  $name = $ref_lt if ($name =~ m/^group_/ && $ref_lt ne "NONE"); 
  
  foreach my $strain (@all_strains) {
    my $idx = $strain_idx->{$strain};
    $ann->{$name}->{$strain}->{lt} = $t[$idx];
  }
  $ann->{$name}->{type} = $type; 
  $ann->{$name}->{loc} = $loc; 
}

#print Dumper $ann;

foreach my $strain (@study_strains) { 
  my $counts_file = join "/",$wdir,"exp_tables",$strain; 
  my $tpm_file = join "/",$wdir,"exp_tables",$strain; 
  $counts_file = join ".",$counts_file,"counts.tsv";
  $tpm_file = join ".",$tpm_file,"TPM.tsv";
  open COUNTS,"<",$counts_file or die "$!";  
  open TPMS,"<",$tpm_file or die "$!";  
  
  my $samples = <TPMS>; 
  $samples = <COUNTS>;
  chomp $samples;
  $samples =~ s/$strain\t//g; 
  $counts->{$strain}->{samples} = $samples; 
  $tpms->{$strain}->{samples} = $samples; 
  
  while (<COUNTS>) { 
    chomp; 
    m/^(.*?)\t(.*)/;
    my $lt = $1;
    my $count_str = $2; 
    $counts->{$strain}->{$lt} = $count_str; 
  } 
  while (<TPMS>) { 
    chomp; 
    m/^(.*?)\t(.*)/;
    my $lt = $1;
    my $tpm_str = $2; 
    $tpms->{$strain}->{$lt} = $tpm_str; 
  }
  
  close COUNTS; 
  close TPMS; 
}

#print Dumper($tpms);   

open ALL_COUNTS,">","Master_table.counts.tsv" or die "$!";
open   ALL_TPMS,">","Master_table.TPM.tsv"    or die "$!";

## and now, the actual table. Serously, fuck this thing sideways. 
my $print_header = "Gene_name\tType\tLocation"; 
foreach my $strain (@all_strains) { 
  $print_header = join "\t",$print_header,$strain;
}
foreach my $strain (@study_strains) { 
  ## @study_strains and @ref_strains are always Unix-sorted
  $print_header = join "\t",$print_header,$counts->{$strain}->{samples};  
}

print ALL_COUNTS "$print_header\n";
print ALL_TPMS   "$print_header\n";

my @names = keys %{$ann}; 
@names   = sort { $a cmp $b } @names;

foreach my $name (@names) { 
  ## flag defines if this gene is present in either of @study_strains
  ## genes absent in every study strain are not to be printed
  my $flag = 0;
                                         
  my $cnt_output = join "\t",$name,$ann->{$name}->{type},$ann->{$name}->{loc};

  foreach my $strain (@all_strains) { 
    $cnt_output = join "\t",$cnt_output,$ann->{$name}->{$strain}->{lt};
  }
  my $tpm_output = $cnt_output; 
 
  foreach my $strain (@study_strains) {
    my $counts_str = 0;
    my $tpm_str   = 0.000;  
    my $lt = $ann->{$name}->{$strain}->{lt};
    if ($lt ne "NONE") { 
      $flag       = 1;  
      $counts_str = $counts->{$strain}->{$lt}; 
      $tpm_str    = $tpms->{$strain}->{$lt}; 
      $cnt_output = join "\t",$cnt_output,$counts_str;
      $tpm_output = join "\t",$tpm_output,$tpm_str;
    } else {
      ## generate a string of appropriate number of zeroes (counts) or 0.0 for TPMs 
      for (my $j=0; $j < $sample_counts->{$strain} - 1; $j++ ) {
        $counts_str = join "\t",$counts_str,"0";
        $tpm_str    = join "\t",$tpm_str,"0.000";  
      } 
      $cnt_output = join "\t",$cnt_output,$counts_str;
      $tpm_output = join "\t",$tpm_output,$tpm_str;
    }  
  } 
  print ALL_COUNTS "$cnt_output\n" if ($flag);
  print ALL_TPMS   "$tpm_output\n" if ($flag);
}

close ALL_COUNTS; 
close ALL_TPMS; 

close CONFIG; 
close ORTHO; 
