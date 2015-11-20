#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

use vars qw($opt_s $opt_e $opt_c $opt_o $opt_m $opt_k $opt_n $opt_r $opt_t $opt_v);

# Usage
my $usage = "
ReadJellyfish.pl - Finds Differentially Present Kmers in Jellyfish Output
Copyright (C) 2014 by Jacob A Tennessen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

This program reads output from Jellyfish 1.0.2.
Reference:
Marçais G, Kingsford C. 2011.
A fast, lock-free approach for efficient parallel counting of occurrences of k-mers.
Bioinformatics. 27(6):764-70. doi: 10.1093/bioinformatics/btr011.

This program assumes a set of files in which the first column is a kmer sequence and the second column is a count
Files of this type can be generated with the command 'jellyfish count' folowed by 'jellyfish stats'
One of these files is designated the 'target' sample from which candidate kmers will be extracted
Some of these files are designated 'experimental samples' and others 'control samples'
This program finds candidate kmers that are significantly more common in the experimental samples than in the control samples

Usage: perl ReadJellyfish.pl options
 required:
  -s	the target sample name
  -e	a file listing names of experimental samples
  -c	a file listing names of control samples
  -o	output file name, will be used in tab-delimited column file and FASTA file
 optional:
  -m	maximum size of hash containing candidate kmers [default = 10000000]
  -k	kmer size used [default = 31]
  -n	number of times a kmer must be seen in target sample [default = 2]
  -r	minimum ratio of kmer count in experimental samples to control samples [default = 2]
  -t	minimum t-statistic between count in experimental samples and count in control samples [default = 5]
  -v	minimum proportion of experimental samples which must contain kmer [default = 0.5]
";

#############

# command line processing.
getopts('s:e:c:o:m:k:n:r:t:v');
die $usage unless ($opt_s);
die $usage unless ($opt_e);
die $usage unless ($opt_c);
die $usage unless ($opt_o);

my $hashsize;
my $kmersize;
my $indcountthresh;
my $ratiothresh;
my $tstatthresh;
my $minrefinds;

$hashsize = $opt_m ? $opt_m : 10000000;
$kmersize = $opt_k ? $opt_k : 31;
$indcountthresh = $opt_n ? $opt_n : 2;
$ratiothresh = $opt_r ? $opt_r : 2;
$tstatthresh = $opt_t ? $opt_t : 5;
$minrefinds = $opt_v ? $opt_v : 0.5;

my %ref_inds;

open(EXP, $opt_e) || die "can't open $opt_e\n";

while (<EXP>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    $ref_inds{$line} = 1;
}

close (EXP);

my %check_inds;

open(CON, $opt_c) || die "can't open $opt_c\n";

while (<CON>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    $check_inds{$line} = 1;
}

close (CON);

unless ( open(OUT, ">$opt_o") ) {
    print "Cannot open file \"$opt_o\" to write to!!\n\n";
    exit;
}
print OUT "Kmer\tCount\tExpInds\tExpDepth\tControlInds\tControlDepth\tRatio\tTstat\n";

my $bookmarkkmer = "n";

my %out;

my $seeninds = 0;

my @ref_inds = (keys %ref_inds);

my $reftotal = scalar (@ref_inds);

my @check_inds = (keys %check_inds);

my $checktotal = scalar (@check_inds);

until ($seeninds >= 1) {
    
    my %kmers;

    open(IN, $opt_s) || die "can't open $opt_s\n";
    while (<IN>) {
	my $line = $_;
	$line =~ s/\r|\n//g;
	unless ($bookmarkkmer =~ /^n$/) {
	    if ($line =~ /^$bookmarkkmer/) {
		print "Bookmark re-found at $line\n";
		$bookmarkkmer = "n";
	    } else {
		next;
	    }
	}
	my @data = split /\s+/, $line;
	unless (defined $out{$data[0]}) {
	    if ($data[1] >= $indcountthresh) {
		$kmers{$data[0]} = "0\n0";
	    }
	    if ((scalar (keys %kmers)) > $hashsize) {
		$bookmarkkmer = $data[0];
		last;
	    }
	}
    }
    
    close (IN);
    
    if ($bookmarkkmer =~ /^n$/) {
	$seeninds += 1;
    } else  {
	print "Max hash size met at $opt_s $bookmarkkmer\n";
    }

    print "Kmers uploaded.\n";
    
    foreach my $rind (@ref_inds) {
	
	open(RIN, $rind) || die "can't open $rind\n";
	while (<RIN>) {
	    my $line = $_;
	    $line =~ s/\r|\n//g;
	    my @data = split /\s+/, $line;
	    if (defined $kmers{$data[0]}) {
		my @kmerdata = split "\n", $kmers{$data[0]};
		if ($kmerdata[0] =~ /^0/) {
		    $kmerdata[0] = $data[1];
		} else {
		    $kmerdata[0] = "$kmerdata[0]\t$data[1]";
		}
		$kmers{$data[0]} = join "\n", @kmerdata;
	    }
	}
	
	close (RIN);
	
    }
    
    print "Kmer counts uploaded.\n";
    
    foreach my $cind (@check_inds) {
	
	open(CIN, $cind) || die "can't open $cind\n";
	while (<CIN>) {
	    my $line = $_;
	    $line =~ s/\r|\n//g;
	    my @data = split /\s+/, $line;
	    if (defined $kmers{$data[0]}) {
		my @kmerdata = split "\n", $kmers{$data[0]};
		if ($kmerdata[1] =~ /^0/) {
		    $kmerdata[1] = $data[1];
		} else {
		    $kmerdata[1] = "$kmerdata[1]\t$data[1]";
		}
		$kmers{$data[0]} = join "\n", @kmerdata;
	    }
	}
	
	close (CIN);
	
    }
    
    my $outcount = 0;
    
    foreach my $k (keys %kmers) {
	my @kmerdata = split "\n", $kmers{$k};
	my @refdepths = split "\t", $kmerdata[0];
	my @checkdepths = split "\t", $kmerdata[1];
	my $refindcount = 0;
	my $checkindcount = 0;
	my $refdepthtotal = 0;
	my $checkdepthtotal = 0;
	unless ($refdepths[0] =~ /^0/) {
	    $refindcount = scalar (@refdepths);
	    foreach my $r (@refdepths) {
		$refdepthtotal += $r;
	    }
	}
	unless ($checkdepths[0] =~ /^0/) {
	    $checkindcount = scalar (@checkdepths);
	    foreach my $c (@checkdepths) {
		$checkdepthtotal += $c;
	    }
	}
	if ($refindcount >= ($reftotal*$minrefinds)) {
	    until ((scalar (@refdepths)) >= $reftotal) {
		push @refdepths, 0;
	    }
	    my $refmean = $refdepthtotal/$reftotal;
	    my $refstd = stdev (@refdepths);
	    until ((scalar (@checkdepths)) >= $checktotal) {
		push @checkdepths, 0;
	    }
	    my $checkmean = $checkdepthtotal/$checktotal;
	    my $checkstd = stdev (@checkdepths);
	    my $sxx = sqrt((($refstd*$refstd)/$reftotal)+(($checkstd*$checkstd)/$checktotal));
	    my $tstat = 1000;
	    if ($sxx > 0) {
		$tstat = sprintf "%.2f", (abs($refmean-$checkmean))/$sxx;
	    }
	    my $ratio = 1000;
	    if ($checkdepthtotal > 0) {
		$ratio = ($refdepthtotal*$checktotal)/($checkdepthtotal*$reftotal);
	    }
	    if (($ratio >= $ratiothresh)&&($tstat > $tstatthresh)) {
		$outcount +=1;
		$out{$k} = "$outcount\t$refindcount\t$refdepthtotal\t$checkindcount\t$checkdepthtotal\t$ratio\t$tstat";
	    }
	}
    }
    
    my $outsize = (scalar(keys %out));
    
    print "$outsize kmers now in output hash.\n";
    
}

foreach my $o (keys %out) {
    print OUT "$o\t$out{$o}\n";
}

close (OUT);

my $fastafile = "$opt_o.fas";
unless ( open(FASTA, ">$fastafile") ) {
    print "Cannot open file \"$fastafile\" to write to!!\n\n";
    exit;
}

foreach my $o (keys %out) {
    print FASTA ">K$out{$o}\n$o\n";
}

close (FASTA);

#############################

sub average{
        my(@data) = @_;
        if (!defined $data[0]) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach my $d (@data) {
                $total += $d;
        }
        my $average = $total / (scalar(@data));
        return $average;
}

sub stdev{
    my(@data) = @_;
    if(!defined $data[1]){
	    return 0;
    }
    my $average = average(@data);
    my $sqtotal = 0;
    foreach my $d (@data) {
	    $sqtotal += ($average-$d) ** 2;
    }
    my $std = ($sqtotal / ((scalar(@data))-1)) ** 0.5;
    return $std;
}

