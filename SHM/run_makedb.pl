#!/usr/bin/perl -w
use strict;
use Cwd;

my $path = getcwd;

system "mkdir -p makedb";

my $image = "docker.io/immcantation/suite:4.4.0";

my $input = $ARGV[0];
$input=~/(.*?).blastout/;
my $output = "$1";

my $cmd_makedb = "podman run --rm -w $path -v /data:/data $image /usr/local/bin/MakeDb.py igblast -i $ARGV[0] -s $ARGV[1] -r imgt_ref --failed --extended --format changeo --outdir makedb --outname $output\.DB --log $output\.log > makedb/$output\.DB.out";


print "$cmd_makedb\n";
system "$cmd_makedb";
