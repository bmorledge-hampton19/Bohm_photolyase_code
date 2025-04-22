#!/usr/bin/perl

use strict;
use warnings;

# get trx subsets of genes
my $greater_file = "../greater10perhr.txt";
open( GREATER, $greater_file ) || die "Couldn't open file: $greater_file\n";

my %greaterlist;
while ( <GREATER> )
{
	chomp $_;
	$greaterlist{$_} = 1;
}

close( GREATER );

my $mid_file = "../1to10perhr.txt";
my %midlist;
open( MID, $mid_file ) || die "Couldn't open file: $mid_file\n";

while ( <MID> )
{
        chomp $_;
        $midlist{$_} = 1;
}

close( MID );

my $lesser_file = "../less1perhr.txt";
my %lesslist;
open( LESS, $lesser_file ) || die "Couldn't open file: $lesser_file\n";

while ( <LESS> )
{
        chomp $_;
        $lesslist{$_} = 1;
}

close( LESS );


# data matrix
print STDERR "Please enter filename for data matrix\n";
my $filename = <STDIN>;
chomp $filename;

open( FILE, $filename ) || die "Couldn't open file: $filename\n";

# number of bins to compute
my $numgenebins = 6;
my $numflankbins = 3;
my $interval = 167; # size of flanking interval bins
my $flank_offset = $numflankbins * $interval;
my $totalbins = $numgenebins + 2 * $numflankbins;

my %nmp_sum;
my %purine_sum;
my $list = "";
while ( my $line = <FILE> )
{
	chomp $line;
	if ( $line =~ /^(Y[A-P][LR][0-9]{3}[CW]\-?[A-H]?)/ )
	{
		my $acc = $1;
		my @fields = split /\t/, $line;	
		my $prevlist = $list;
		if ( scalar @fields != ( 2 * $totalbins + 2 ) )
		{	
			die "Wrong number of bins for gene: $acc\n";
		}
	
		if ( $greaterlist{$acc} )
		{
			$list = "greater";
		}
		elsif ( $midlist{$acc} )
		{
			$list = "mid";
		}
		elsif ( $lesslist{$acc} )
		{
			$list = "less";
		}
		else
		{
			$list = "other";
		}
		
		if ( $fields[1] eq "CPDs" )
		{
			for( my $i = 2; $i < scalar @fields; $i++ )
			{
				$nmp_sum{$list}[$i - 2] += $fields[$i];
			}
		}
		elsif ( $fields[1] eq "Dipyrimidines" )
		{
			if ( $list ne $prevlist )
			{
				die "Gene list assignment for purines different than assignment for CPDs for acc: $acc\n";
			}
                        for( my $i = 2; $i < scalar @fields; $i++ )
                        {
                                $purine_sum{$list}[$i - 2] += $fields[$i];
                        }
                }
		else
		{
			die "Count type for acc $acc is $fields[1]\n";
		}
	}		
}
my @list_types = ("greater", "mid", "less", "other");
print "From file: $filename\t for transcription frequency subsets\n";

# print results:
print "Data Type\tGene List";
for (my $i = 0; $i < $numflankbins; $i++ )
{
        my $start = 0 - $flank_offset + ( $i * $interval);
        my $end = 0 - $flank_offset + (($i + 1) * $interval) - 1;
	print "\tTS Promoter ($start to $end)";
}
for ( my $i = 1; $i <= $numgenebins; $i++ )
{
	print "\tTS Coding bin $i";
}
for ( my $i = 0; $i < $numflankbins; $i++ )
{
	my $start = ( $i * $interval) + 1; 
        my $end = ($i + 1) * $interval;
	print "\tTS Terminator ($start to $end)";
}
foreach my $type (@list_types)
{
	print "\nCPDs\t$type";
	my $midway = (scalar @{$nmp_sum{$type}}) / 2;
	for (my $i = 0; $i < $midway; $i++)
	{
		print "\t$nmp_sum{$type}[$i]";
	}
	print "\nDipyrimidines\t$type";
	if ( $midway != scalar (@{$purine_sum{$type}}) / 2 )
	{	
		die "Error NMP and Purine arrays are of different sizes!\n";
	}
	elsif ( $midway != ($numgenebins + 2 * $numflankbins ) )
	{
		die "Mismatch in number of bins in input matrix and sum_all file\n";
	}
	for (my $i = 0; $i < $midway; $i++)
	{
	        print "\t$purine_sum{$type}[$i]";
	}
	print "\nNormalized CPDs\t$type";
	for (my $i = 0; $i < $midway; $i++)
	{
		my $avg = 1.0 * $nmp_sum{$type}[$i]/$purine_sum{$type}[$i];
	        print "\t$avg";
	}
	print "\n\n";
}
print "Data Type\tGene List";
for (my $i = 0; $i < $numflankbins; $i++ )
{
        my $start = 0 - $flank_offset + ( $i * $interval);
        my $end = 0 - $flank_offset + (($i + 1) * $interval) - 1;
        print "\tNTS Promoter ($start to $end)";
}
for ( my $i = 1; $i <= $numgenebins; $i++ )
{
        print "\tNTS Coding bin $i";
}
for ( my $i = 0; $i < $numflankbins; $i++ )
{
        my $start = ( $i * $interval) + 1; 
        my $end = ($i + 1) * $interval;
        print "\tNTS Terminator ($start to $end)";
}
foreach my $type (@list_types)
{
	print "\nCPDs\t$type";
        my $midway = (scalar @{$nmp_sum{$type}}) / 2;
        if ( $midway != scalar (@{$purine_sum{$type}}) / 2 )
        {
                die "Error NMP and Purine arrays are of different sizes!\n";
        }
        elsif ( $midway != ($numgenebins + 2 * $numflankbins ) )
        {
                die "Mismatch in number of bins in input matrix and sum_all file\n";
        }

	for (my $i = $midway; $i < scalar @{$nmp_sum{$type}}; $i++)
	{
	        print "\t$nmp_sum{$type}[$i]";
	}
	print "\nDipyrimidines\t$type";
	for (my $i = $midway; $i < scalar @{$purine_sum{$type}}; $i++)
	{
	        print "\t$purine_sum{$type}[$i]";
	}
	print "\nNormalized CPDs\t$type";
	for (my $i = $midway; $i < scalar @{$nmp_sum{$type}}; $i++)
	{
	        my $avg = 1.0 * $nmp_sum{$type}[$i]/$purine_sum{$type}[$i];
	        print "\t$avg";
	}
	print "\n\n";
}
