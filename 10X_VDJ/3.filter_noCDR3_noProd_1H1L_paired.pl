#!/usr/bin/perl -w
use strict;

system "mkdir -p filtered_paired";

my @files = <*.blastout.productive>;
foreach my $f(@files)
{
  my %num_chains;

  $f=~/(.*).blastout/;
  my $capture = $1; 

  open OUT,">filtered_paired/$capture\.filtered_paired" or die $!;
  
  open IN,"$f" or die $!;
  while(my $ln = <IN>)
  {
    if($ln=~/sequence_id/){chomp $ln;print OUT "Cell\t$ln\t$ln\n";}
    else
    {
      my @arr = split("\t",$ln);
      my $prod = $arr[7];
      my $cdr3 = $arr[55];
      $arr[0]=~/(.*)_contig/;
      my $cell = $1;

      # print "$cell\t$prod\t$cdr3\n";


      if($prod eq "T" && $cdr3=~/[A-W]/)
      {
        if($arr[3] eq "IGH"){$num_chains{$cell}{'Heavy'}++;}
        if($arr[3] eq "IGL" || $arr[3] eq "IGK"){$num_chains{$cell}{'Light'}++;}
      }
    }
  }
  close IN;
 
  my %shortlist;

  open OUT1,">$capture\_num_chains" or die $!;
  for my $cell(keys %num_chains)
  {
    if(!($num_chains{$cell}{'Heavy'})){$num_chains{$cell}{'Heavy'}=0;}
    if(!($num_chains{$cell}{'Light'})){$num_chains{$cell}{'Light'}=0;}

    print OUT1 "$cell\t$num_chains{$cell}{'Heavy'}\t$num_chains{$cell}{'Light'}\n";

    if($num_chains{$cell}{'Heavy'} == 1 && $num_chains{$cell}{'Light'} == 1){$shortlist{$cell} = 1;}
  }
  close OUT1;

  my %chains;
  open IN,"$f" or die $!;
  while(my $ln = <IN>)
  {
    chomp $ln;
    my @arr = split("\t",$ln);
    $arr[0]=~/(.*)_contig/;
    my $cell = $1;
    if(!($shortlist{$cell})){next;}
    my $type = "";
    if($arr[3] eq "IGH"){$type = "Heavy";}
    elsif($arr[3] eq "IGK" || $arr[3] eq "IGL"){$type = "Light";}

    $chains{$cell}{$type} = $ln;
  }
  close IN;


  foreach my $c(keys %chains)
  {
    if(!($chains{$c}{'Heavy'})){print "No heavy chains for $c in $capture\n";}
    if(!($chains{$c}{'Light'})){print "No light chains for $c in $capture\n";}
  
    if($chains{$c}{'Heavy'} && $chains{$c}{'Light'}) 
    {
      print OUT "$c\t$chains{$c}{'Heavy'}\t$chains{$c}{'Light'}\n";
    }
  }
  close OUT;
}
