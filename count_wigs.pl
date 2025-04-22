#!/usr/bin/perl

use strict;
use warnings;

my $count = 0;
my $vals = 0;

while ( my $line = <STDIN> )
{
	chomp $line;
	if ( $line =~ /^variableStep/ )
	{
		next;
	}

	my @fields = split /\t/, $line;
	$count += $fields[1];
	$vals++;	
}
my $avg = 1.0 * $count / $vals;
print "Total reads: $count\nTotal positions: $vals\nAvg read per position: $avg\n";

print STDERR "reads=$count\n";
