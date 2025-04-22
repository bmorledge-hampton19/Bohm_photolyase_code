#!/usr/bin/perl

use strict;
use warnings;

print "Please enter filename of matrix to split ts and nts strand data\n";

my $file = <STDIN>;
chomp $file;

my $tsfile = $file;
$tsfile =~ s/\.txt/_ts\.txt/ || die "Wrong file type!\n";
open (TS, ">$tsfile") || die "couldn't open file $tsfile\n";

my $ntsfile = $file;
$ntsfile =~ s/\.txt/_nts\.txt/ || die "Wrong file type!\n";
open (NTS, ">$ntsfile") || die "couldn't open file $ntsfile\n";

open (DATA, "$file") || die "Couldn't open file $file\n";

my $firstline = 1;
my $prevmid = 0;
while ( my $line = <DATA> )
{
	chomp $line;

	if ( $line =~ /Bin size/ )
	{
		next;
	}

	my @fields = split /\t/, $line, -1;

	my $numfields = scalar @fields;

	if ( $numfields % 2 == 0 )
	{
		die "Wrong number of columns!\n";
	}

	my $mid = int ( $numfields / 2 );
	my $ts_end = $mid;
	my $nts_start = $mid + 1;

	if ( $firstline )
	{
		if ( $fields[$ts_end] =~ /\(TS\)/ && $fields[$nts_start] =~ /\(NTS\)/ )
		{
			;
		}
		else
		{
			die "error in finding column midpoint\n";
		}
	}
	elsif ( $mid != $prevmid )
	{
		die "Error finding mid in one of columns\n";
	}

	#print ts
	print TS "$fields[0]";
	for ( my $i = 1; $i <= $ts_end; $i++ )
	{
		print TS "\t$fields[$i]";
	}
	print TS "\n";

	#print nts
	print NTS "$fields[0]";
	for ( my $i = $nts_start; $i < $numfields; $i++ )
	{
                print NTS "\t$fields[$i]";
        }
        print NTS "\n";

	$prevmid = $mid;
	$firstline = 0;
}
