#!/usr/bin/perl -w
use strict;

use Set::Scalar;
use Graph::UnionFind;

system "mkdir -p Clones_MB";

my %cluster;

my $input = $ARGV[0];
my $thr = $ARGV[1]; # =< 1

$input=~/(.*)_Ag_VDJ_combined.txt/;
my $id = $1;
my $output = "$id\_clones_$thr.txt";

sub compare_arr
{
	my ($str_1,$str_2) = @_;
	my $match = 0;

	my @arr1 = split(",",$str_1);
	my @arr2 = split(",",$str_2);

	foreach my $a1(@arr1)
	{
		foreach my $a2(@arr2)
		{
			$a1=~s/\*.*//;
			$a2=~s/\*.*//;
			if($a1 eq $a2)
			{
				$match = 1;
				return $match;
			}
		}
	}
	return $match;
}


sub compare_cdr3
{
	my ($cdr3_1,$cdr3_2) = @_;
	my $match = 0;
	if(length($cdr3_1) == length($cdr3_2))
	{

		my $len =  length($cdr3_1);
		my $count = ( $cdr3_1 ^ $cdr3_2 ) =~ tr/\0//c;
		#print "$cdr3_1 $cdr3_2 $len $count\n";

		my $per = 1-($count/$len);

		if($per >= $thr)
		{
			$match = 1;
		}
	}
	return $match;
}


my %db;
my %clonotype;



#####################################################################################################################


# Single cell
open IN,"$input" or die $!;
while(my $ln = <IN>)
{
	chomp $ln;
	
	$ln=~s/\r|\n//g;
	my @arr = split("\t",$ln);	

	my $new_id = $arr[0]; 
	
	$db{"MB\_$new_id"} = "MB\_$new_id\t$ln";
	
	my $hcdr3_len = length($arr[56]);

	if($arr[56] && $arr[56]=~/[A-Z]/ && $arr[167] && $arr[167]=~/[A-Z]/)
	{
		$clonotype{$hcdr3_len}{"MB\_$new_id"}{'HV'} = $arr[12];
		$clonotype{$hcdr3_len}{"MB\_$new_id"}{'HJ'} = $arr[14];
		$clonotype{$hcdr3_len}{"MB\_$new_id"}{'HCDR3'} = $arr[56];
		$clonotype{$hcdr3_len}{"MB\_$new_id"}{'LV'} = $arr[123];
		$clonotype{$hcdr3_len}{"MB\_$new_id"}{'LJ'} = $arr[125];
		$clonotype{$hcdr3_len}{"MB\_$new_id"}{'LCDR3'} = $arr[167];

		# print "$arr[0]\t$arr[12]\t$arr[14]\t$arr[56]\t$arr[123]\t$arr[125]\t$arr[167]\n";
	}
}
close IN;


my $uf100 = Graph::UnionFind->new;
my %v100;

my $ctr1 = 0;
foreach my $len(sort keys %clonotype)
{
	if($len == 0){next;}
	foreach my $s1 (sort keys %{$clonotype{$len}})
	{
		foreach my $s2 (sort keys %{$clonotype{$len}})
		{
			# print "$s1 $s2\n";
        	my $hv_check = compare_arr($clonotype{$len}{$s1}{'HV'},$clonotype{$len}{$s2}{'HV'});
			if($hv_check)
			{
				# print "HV match\n";
				my $hj_check = compare_arr($clonotype{$len}{$s1}{'HJ'},$clonotype{$len}{$s2}{'HJ'});
				if($hj_check)
				{
					# print "HJ match\n";
					my $hcdr3_check = compare_cdr3($clonotype{$len}{$s1}{'HCDR3'},$clonotype{$len}{$s2}{'HCDR3'});
					if($hcdr3_check)
					{		
						# print "HCDR3 match\n";
						my $lv_check = compare_arr($clonotype{$len}{$s1}{'LV'},$clonotype{$len}{$s2}{'LV'});
						if($lv_check)
						{
							# print "LV match\n";
							my $lj_check = compare_arr($clonotype{$len}{$s1}{'LJ'},$clonotype{$len}{$s2}{'LJ'});
							if($lj_check)
							{
								# print "LJ match\n";
								my $lcdr3_check = compare_cdr3($clonotype{$len}{$s1}{'LCDR3'},$clonotype{$len}{$s2}{'LCDR3'});
								if($lcdr3_check)
								{	
									# print "Found match\n";	
									++$v100{$_} for $s1,$s2;
			    					$uf100->union($s1,$s2);
								}
							}
						}			
			
					}
				}
			}
			
		}
	}
}

my %c100;
foreach my $v (keys %v100)
{
	#print "Vertex $v\n";
    my $b = $uf100->find($v);
    die "$0: no block for $v" unless defined $b;
    # print "$b\t$v\n";	
    push @{$c100{$b}}, $v;
}


my %clones;
my $ctr = 0;
foreach my $a(keys %c100)
{
	$ctr++;
	foreach my $seq(@{$c100{$a}})
	{
		$clones{$seq} = "$id\_$ctr";
	}
	
}

open OUT,">Clones_MB/$output" or die $!;
foreach my $seq(keys %db)
{
	# print "$db{$seq}\t$clones{$seq}\n";
	if($clones{$seq})
	{
		print OUT "$clones{$seq}\t$db{$seq}\n";
	}
}
