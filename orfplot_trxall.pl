#!/usr/bin/perl

use strict;
use warnings;

# get trx subsets of genes
my $greater_file = "../trxfreq_greater10perhr.txt";
open( GREATER, $greater_file ) || die "Couldn't open file: $greater_file\n";

my %genelist;
while ( <GREATER> )
{
	chomp $_;
	my @fields = split /\t/, $_;
	$genelist{$fields[0]} = $fields[1];
}

close( GREATER );

my $mid_file = "../trxfreq_1to10perhr.txt";
open( MID, $mid_file ) || die "Couldn't open file: $mid_file\n";

while ( <MID> )
{
        chomp $_;
        my @fields = split /\t/, $_;
        $genelist{$fields[0]} = $fields[1];
}

close( MID );

my $lesser_file = "../trxfreq_less1perhr.txt";
open( LESS, $lesser_file ) || die "Couldn't open file: $lesser_file\n";

while ( <LESS> )
{
        chomp $_;
        my @fields = split /\t/, $_;
        $genelist{$fields[0]} = $fields[1];
}

close( LESS );

print STDERR "Please enter filename for data matrix\n";
my $filename = <STDIN>;
chomp $filename;

my $fileout = $filename;
$fileout =~ s/\.txt/_trxsorted\.txt/ || die "filename of data matrix is wrong type!\n";
open ( OUT, ">$fileout") || die "couldn't open file\n";

open( FILE, $filename ) || die "Couldn't open file: $filename\n";

if ( $filename =~ /_centered/ || $filename =~ /_center/ || $filename =~ /_norm/ )
{
	;
}
else
{
	# remove top header because (normally removed by logtransform_center.pl program)
	my $topline = <FILE>;
}

my $header = <FILE>;
my %matrix;
while ( my $line = <FILE> )
{
	chomp $line;
	my @fields = split /\t/, $line;
	$matrix{$fields[0]} = $line;
}

#print less list
print OUT $header;
foreach my $acc ( sort {$genelist{$b} <=> $genelist{$a} or $a cmp $b } keys %genelist )
{
        if ( exists $matrix{$acc} )
        {
		print OUT $matrix{$acc} . "\n";
	}
}

