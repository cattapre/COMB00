#!/usr/bin/perl
# ChIPseq_wrapper.pl by Rinaldo Catta-Preta @UCDavis

use strict;
use warnings 'FATAL' => 'all';

my ($bayes_dir, $edit_file) = @ARGV;
$edit_file = $bayes_dir . '/' . $edit_file;
my $outfile = $edit_file . ".bed";

open (my $in, '<' , $edit_file) or die "Can't open $edit_file for reading\n";
open (my $out, '>' , $outfile) or die "Can't open $outfile for writing\n";

my $number = 0;
my @editfile = split /\//, $edit_file;

while (<$in>)
{
	chomp;
	next if /^\s*space/;
	
	my @data = split / /;
	
	my $chrom = $data[1];
	my $chromStart = $data[2];
	my $chromEnd = $data[3];
	my $peak = "$editfile[-1]" . '_' . $number;
	my $score = $data[5] * 1000;
	my $strand = '.';
	my $tStart = $chromStart;
	my $tEnd = $chromStart;
	my $RGB = 0;
	my $bCount = 1;
	my $bSizes = '0';
	my $bStarts = '0';
	
	print $out "$chrom\t$chromStart\t$chromEnd\t$peak\t$score\t$strand\t$tStart\t$tEnd\t$RGB\t$bCount\t$bSizes\t$bStarts\n";
	$number+= 1;
}

close $out;
close $in;

print "ready\n";

