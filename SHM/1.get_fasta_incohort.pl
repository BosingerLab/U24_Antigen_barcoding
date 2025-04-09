#!/usr/bin/perl -w
use strict;

open IN,"MUSA.txt" or die $!;

my %done;

while(my $ln = <IN>)
{
	chomp $ln;

	if($ln=~/^allele/) {next;}

	my @arr = split("\t",$ln);

	my $gene_type = $arr[2];
	my $new_allele = $arr[0];
	my $seq = $arr[4];
	my $sample_genomic = $arr[22];
	my $sample_rep = $arr[24];
	
	if($done{$new_allele}){next;}


	my @genomic = split(", ",$sample_genomic);
	my @airr = split(", ",$sample_rep);

	my %seen = map { $_ => 1 } @genomic;

	# Count common elements
	my $common_count = 0;
	my %counted;

	foreach my $aid (@airr) {
	    if ($seen{$aid} && !$counted{$aid}) {
	        $common_count++;
	        $counted{$aid} = 1;  # To avoid counting duplicates
	    }
	}

	if($common_count>0)
	{
		open OUT1,">>$gene_type.fasta";
		print OUT1 ">$new_allele\n$seq\n";
		$done{$new_allele} = 1;
	}
}
