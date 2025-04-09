#!/usr/bin/perl -w
use strict;

my %c_aid;

open IN,"$ARGV[0]" or die $!; #Public_clones_0.7.txt 
while(my $ln = <IN>)
{
  my @arr = split("\t",$ln);
  $arr[1] =~/\-1_(.*)/;
  my $aid = $1;
  $c_aid{$arr[0]}{$aid}=1; 
}
close IN;

my %public;
foreach my $c(keys %c_aid)
{
  my $num = scalar(keys %{$c_aid{$c}});
  {
    if($num > 1)
    {
      $public{$c} = 1;
    }
  }
}

open IN,"$ARGV[0]" or die $!;
open OUT,">$ARGV[1]" or die $!;

while(my $ln = <IN>)
{
  my @arr = split("\t",$ln);
  if($public{$arr[0]})
  {
    my $new = join ",", keys %{$c_aid{$arr[0]}};
    print OUT "$new\t",$ln;
  }
}
close IN;

