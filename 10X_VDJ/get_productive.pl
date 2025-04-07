#!/usr/bin/perl -w
use strict;

system "mkdir -p productive";

my @files = <*.blastout>;
foreach my $f(@files)
{
  open OUT,">productive/$f\.productive" or die $!;

  open IN,"$f" or die $!;
  while(my $ln = <IN>)
  {
    if($ln=~/sequence_id/){print OUT "$ln";}
    else
    {
      my @arr = split("\t",$ln);
      my $prod = $arr[7];
      my $cdr3 = $arr[55];
      if($prod eq "T" && $cdr3=~/[A-W]/)
      {
	    print OUT "$ln";
	  }
    }
  }
  close IN;
  close OUT;
}
