#!/usr/bin/perl

use strict;
use warnings;

my $infile=shift;
open (my $fh ,'<',$infile) or die $!;

my $header = <$fh>;
my $depth_thres = shift;

while (<$fh>) {
    chomp;
    my @row = split "\t";
    my @depths = @row[2..$#row];
    my @n_minDP = grep $_ >= $depth_thres,@depths;
    print $row[0],"\t",$row[1],"\t",scalar(@n_minDP), "\n";
}
