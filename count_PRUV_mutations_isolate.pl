#!/usr/bin/perl -w

use strict;
use warnings;

my $header = <STDIN>;
print STDERR "Counting isolates - Header: $header\n";

my %isolate = ();
my $skippedlines = 0;
while( my $line = <STDIN> )
{
	chomp $line;
	if ( $line eq "")
	{
		next;
	}
	my @fields = split /\t/, $line;

	my $id = $fields[0];
	my $genotype = "";
	my $isolatename = "";
	if ( $id =~ /([a-zA-Z0-9]+_[a-zA-Z0-9]+)_([0-9]+)$/ )
	{
		$isolatename = $2;
		$genotype = $1;
	
		my $uv = "_UVC";
		$genotype .= $uv; 
	}
	else
	{
		die "error with id: $id\n";
	}
=pod
	# remove D15 Ogg1d isolate (low sequencing depth)
	if ( $fields[0] eq "D15 OGG UV" )
	{
		print STDERR "skipped line: $fields[0]\n";
		$skippedlines++;
		next;
	}
=cut
        my $chr = $fields[2];

        $chr =~ s/\smapping$// || die "weird chr: $chr\n";
        $chr = "chr" . $chr;	

	if ( $chr eq "chrMito" )
	{
        	print STDERR "Skipped line: ${chr}\n";
	        $skippedlines++;
        	next;
	}
		
	$isolate{$genotype}{$isolatename}{$fields[4]}++;
}

my $strain = "";
my $basename = "_Mutcount.txt";
my $filename = "";
foreach my $gene ( sort keys %isolate )
{
        if ( $gene ne $strain )
        {
                if ( $gene =~ /[A-Za-z0-9]+/ )
                {
                        if ( $strain ne "" )
                        {
                                close ( OUT );
				close ( ALL );
                        }
                        $strain = $gene;
                        $filename = $strain . $basename;
                        open ( OUT, ">$filename" ) or die "cannot open file: $filename\n";
			print OUT "Isolate count for $gene\n";
			print OUT "Isolate number\tIsolate name\tMutation type\tMutation count\n";

			my $allfile = $strain. "_all" . $basename;
                        open ( ALL, ">$allfile" ) or die "cannot open file: $allfile\n";
                        print ALL "Isolate count for $gene\n";
                        print ALL "Isolate number\tIsolate name\tMutation count\n";
                }
        }
	my $isonum = 0;
	foreach my $iso ( sort keys %{$isolate{$gene}} )
	{
		my $totalmut = 0;
		foreach my $type ( sort keys %{$isolate{$gene}{$iso}} )
		{
			print OUT "$isonum\t$iso\t$type\t$isolate{$gene}{$iso}{$type}\n";
			$totalmut += $isolate{$gene}{$iso}{$type};
		}
		print OUT "$isonum\t$iso\tAll\t$totalmut\n";
		print ALL "$isonum\t$iso\t$totalmut\n";	
		$isonum++;
	}
	
}

print STDERR "Skipped lines [just rad16 and rad26 isolates: $skippedlines\n\n";
