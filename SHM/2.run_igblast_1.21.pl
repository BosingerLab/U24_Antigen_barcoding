#!/usr/bin/perl -w

my $igblast = "/data/runs/tools/igblast/ncbi-igblast-1.21.0/";
my $nproc = 32;

my $input = $ARGV[0];
$input=~/(.*).fa/;
my $output = "$1";

system "mkdir -p IgBLAST/clonotype";

system "ln -sf $igblast/internal_data .";

my $db = "/data/runs/Analysis_BLab/U24_Chris/Analysis/clonal_analysis/ag_sp/SHM/DB/fasta_incohort";
my $v_file = "V.fasta";
my $d_file = "D.fasta";
my $j_file = "J.fasta";

my $aux = "$db/J.aux";

my $cmd = "$igblast/bin/igblastn -germline_db_V $db/$v_file -germline_db_J $db/$j_file -germline_db_D $db/$d_file -organism rhesus_monkey -domain_system imgt -ig_seqtype Ig -query $input -auxiliary_data $aux -outfmt '7 std qseq sseq btop' -out IgBLAST/$output\.blastout -clonotype_out IgBLAST/clonotype/$output\.clonotype -num_threads $nproc";
print "$cmd\n";

  
