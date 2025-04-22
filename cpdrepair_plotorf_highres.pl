#!/usr/bin/perl

use strict;
use warnings;

use lib '../';
use GeneCoord;
use CPDReadValues;

# ask for probe filename to analyze
print STDERR "Enter filename of plus strand reads for repair (2hr) timepoint\n";
my $repairplusfile = <STDIN>;
chomp($repairplusfile);

print STDERR "Enter filename of minus strand reads for repair (2hr) timepoint\n";
my $repairminusfile = <STDIN>;
chomp($repairminusfile);

print STDERR "Enter filename of plus strand reads for 0hr timepoint\n";
my $plusfile = <STDIN>;
chomp($plusfile);

print STDERR "Enter filename of minus strand reads for 0hr timepoint\n";
my $minusfile = <STDIN>;
chomp($minusfile);

print STDERR "Loading Gene coordinates\n";
my $genes = GeneCoord->new();
print STDERR "Loading 0hr files\n";
my $zeroreads = CPDReadValues->new($plusfile, $minusfile);

print STDERR "Loading repair (2hr) files\n";
my $repaireads = CPDReadValues->new($repairplusfile, $repairminusfile);

# location offsets
my $upstream_offset = -500;
my $downstream_offset = 640;

# 15 bp bins
my $binsize = 15;
#my $binwindow = ( $binsize - 1 )/ 2;
my $binwindow = $binsize;

my %chromosomes = $genes->get_chromosomes();
my %trxstart = $genes->get_tss();
my %trxend = $genes->get_tts();
my %strand = $genes->get_strand();

#print header
print "Fraction Lesions Remaining\tBin size = $binsize\tData from files - 0hr: $plusfile\t$minusfile\t - repair: $repairplusfile\t$repairminusfile\n";
print "YORF";

for (my $i = $upstream_offset; $i <= $downstream_offset; $i += $binsize )
{
	print "\t$i (TS)";
}
for (my $i = $upstream_offset; $i <= $downstream_offset; $i += $binsize )
{
        print "\t$i (NTS)";
}
print "\n";

foreach my $chr (sort keys %chromosomes)
{
	print STDERR "Starting $chr\n";
	my %zeroplusreads = $zeroreads->get_plus_reads_for_chromosome($chr);
	my $num_zeroplusreads = scalar keys %zeroplusreads;
	my %zerominusreads = $zeroreads->get_minus_reads_for_chromosome($chr);
	my $num_zerominusreads = scalar keys %zerominusreads;

        my %repairplusreads = $repaireads->get_plus_reads_for_chromosome($chr);
        my $num_repairplusreads = scalar keys %repairplusreads;
        my %repairminusreads = $repaireads->get_minus_reads_for_chromosome($chr);
        my $num_repairminusreads = scalar keys %repairminusreads;
	print STDERR "$chr 0hr reads: $num_zeroplusreads plus reads and $num_zerominusreads minus reads; repair: $num_repairplusreads plus reads and $num_repairminusreads minus reads\n";
	foreach my $acc ( @{$chromosomes{$chr}} )
	{
		my $tss = $trxstart{$acc};
		my $tts = $trxend{$acc};
		
		my @fracremain_cpd_ts = ();
		my @fracremain_cpd_nts = ();
		# calculate read sums (CPD and DIPY bkgd) for gene
		for ( my $i = $upstream_offset; $i <= $downstream_offset; $i += $binsize)
		{
			my $pos;
			if ( $strand{$acc} eq "+" )
			{
				$pos = $tss + $i;
			}
			elsif ( $strand{$acc} eq "-" )
			{
				$pos = $tss - $i;
			}
			else
			{
				die "No strand information for gene: $acc\n";
			}

			my $zero_plus_cpds = 0;
			my $repair_plus_cpds = 0;
			my $zero_minus_cpds = 0;
			my $repair_minus_cpds = 0;

                        for ( my $j = $pos - $binwindow; $j <= $pos + $binwindow; $j++ )
                        {
				if ( exists $zeroplusreads{$j} )
				{
					$zero_plus_cpds += $zeroplusreads{$j};
				}
				if ( exists $zerominusreads{$j} )
				{
					$zero_minus_cpds += $zerominusreads{$j};
				}

                                if ( exists $repairplusreads{$j} )
                                {
                                        $repair_plus_cpds += $repairplusreads{$j};
                                }
                                if ( exists $repairminusreads{$j} )
                                {
                                        $repair_minus_cpds += $repairminusreads{$j};
                                }      
			}

                        my $plus_remaining = "";
			if ( $zero_plus_cpds > 0 )
			{
				$plus_remaining = 1.0 * $repair_plus_cpds / $zero_plus_cpds;
			}
			my $minus_remaining = "";
			if ( $zero_minus_cpds > 0 )
			{
                        	$minus_remaining = 1.0 * $repair_minus_cpds / $zero_minus_cpds;
			}

                        if ( $strand{$acc} eq "+" )
                        {
				push @fracremain_cpd_nts, $plus_remaining;

				push @fracremain_cpd_ts, $minus_remaining;
			}
			elsif ( $strand{$acc} eq "-" )
                        {
                                push @fracremain_cpd_ts, $plus_remaining;

                                push @fracremain_cpd_nts, $minus_remaining;
			}
		}

		if ( scalar @fracremain_cpd_ts != scalar @fracremain_cpd_nts )
		{
			die "Arrays are of different sizes!\n";
		}

		# bin values

		# print Fraction of CPD's remaining ( unrepaired ) for acc
		print "$acc";

		foreach my $val (@fracremain_cpd_ts)
		{
			print "\t$val";
		}
		foreach my $val (@fracremain_cpd_nts)
		{
			print "\t$val";
		}
		print "\n";
	}

}
