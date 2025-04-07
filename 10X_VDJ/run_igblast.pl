#!/usr/bin/perl -w
use strict;

system "mkdir -p IgBLAST_AIRR/clonotype";
system "mkdir -p IgBLAST_AIRR_KimDB/clonotype";
system "ln -sf /data/runs/tools/vdj_references/rhesus/internal_data_Cirelli internal_data";

my $igblastn = "/data/runs/tools/igblast/ncbi-igblast-1.21.0/bin/igblastn";
my $db = "/data/runs/tools/vdj_references/rhesus/internal_data_Cirelli/rhesus_monkey";
my $constant = "$db/IG_constant.fa";
my $aux_data = "/data/runs/tools/vdj_references/rhesus/internal_data_Cirelli/optional_file/rhesus_monkey_gl.aux";
my $organism = "rhesus_monkey";
my $kimdb = "/data/runs/tools/vdj_references/rhesus/KimDB_v1.1";


my @files=<*.fasta>;
foreach my $f(@files)
{
  $f =~ /(.*)_filtered_contig.fasta/;
  my $pre = $1;
  print "$igblastn -germline_db_V $db/rhesus_monkey_V -germline_db_J $db/rhesus_monkey_J -germline_db_D $db/rhesus_monkey_D -organism $organism -num_threads 8 -domain_system imgt -ig_seqtype Ig -query $f -auxiliary_data $aux_data -c_region_db $constant  -outfmt 19 -out IgBLAST_AIRR/$pre\.blastout -clonotype_out IgBLAST_AIRR/clonotype/$pre\.clonotype\n"; 
 
  print "$igblastn -germline_db_V $kimdb/V.fasta -germline_db_J $kimdb/J.fasta -germline_db_D $kimdb/D.fasta -organism $organism -num_threads 8 -domain_system imgt -ig_seqtype Ig -query $f -auxiliary_data $aux_data -c_region_db $constant  -outfmt 19 -out IgBLAST_AIRR_KimDB/$pre\.blastout -clonotype_out IgBLAST_AIRR_KimDB/clonotype/$pre\.clonotype\n";

}


