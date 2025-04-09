#!/usr/bin/perl -w
use strict;
use Math::Round;

my %count;
my %total;
my %freq;
my %clone_aid;
my %public;
my %positions;
my %overlap;


# Rounding affects the totals for karyotype file. This was fixed manually

open OUT1,">karyotype.txt" or die $!; 
open OUT2,">links.txt" or die $!; 

# Get counts and total for each clone

open IN,"Public_clones_0.7.txt" or die $!;
while(my $ln = <IN>)
{
	chomp $ln;
	my @arr = split("\t",$ln);

	my $cloneid = $arr[0];
	my $cellid = $arr[2];
	$cellid=~/\-1_(.*)/;
	my $aid = $1;

	$count{$aid}{$cloneid}++;
	$total{$aid}++;
	$clone_aid{$cloneid}{$aid} = 1;
}
close IN;


# Calculate frequencies for clones

foreach my $aid (keys %count)
{
	foreach my $cloneid(keys %{$count{$aid}})
	{
		$freq{$aid}{$cloneid} = round(($count{$aid}{$cloneid}/$total{$aid}) * 100000);
		#print "$aid\t$cloneid\t$freq{$aid}{$cloneid}\n";
	}
}


# find public clones

foreach my $c(keys %clone_aid)
{
  my $num = scalar(keys %{$clone_aid{$c}});
  {
    if($num > 1)
    {
      $public{$c} = 1;
    }
  }
}


# Get positions

foreach my $aid(sort keys %freq)
{
	print OUT1 "chr\t-\t$aid\t$aid\t0\t100000\tvvdred\n";
}

foreach my $aid(sort keys %freq)
{
	my $sum = 0;
	my $start = 0;
	my $old = "";

	foreach my $clone(sort {$freq{$aid}{$b} <=> $freq{$aid}{$a}} keys %{$freq{$aid}})
	{
		my $end = $start + $freq{$aid}{$clone};
		print OUT1 "band\t$aid\t$clone\t$clone\t$start\t$end\tgneg\n";
		$positions{$aid}{$clone} = "$start\t$end\t";
		$start = $end;
	}
}

my %done;
foreach my $clone(keys %public)
{
	foreach my $aid1(sort keys %{$clone_aid{$clone}})
	{
		foreach my $aid2(sort keys %{$clone_aid{$clone}})
		{
			if($aid1 ne $aid2)
			{
				if(!($done{$clone}))
				{
					# VH1	0	20	JH5	0	20	color=vvdred
  					print OUT2 "$clone\t$aid1\t$positions{$aid1}{$clone}\t$aid2\t$positions{$aid2}{$clone}\tcolor=vvdgreen\n";		
				}
				$done{$clone} =1;
			}
		}
	}   
}

#system "circos -conf vp.conf";
#system "mv circos.svg $subject\.svg";

