#! /usr/bin/env perl
# Plot time series of tagged components of command line argument in EMEP tagged run
# Version 0: Jane Coates 24/11/2014

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions/MCM/EMEP_tagged_solvents_only_all";
my $mecca = MECCA->new("$base/boxmodel");
my $ntime = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times = $times(1:$ntime-2);

my $species = $ARGV[0];
my @tagged_species = get_all_tagged_species($species);
my %plot_data;

foreach my $tagged (@tagged_species) {
    my $conc = $mecca->tracer($tagged);
    $conc = $conc(1:$ntime-2) * 1e9;
    $plot_data{$tagged} += $conc;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
);

$R->set('Time', [map { $_ } $times->dog]);
$R->run(q` data = data.frame(Time) `);
foreach my $item (keys %plot_data) {
    $R->set('species', $item);
    $R->set('mixing.ratio', [ map { $_ } $plot_data{$item}->dog]);
    $R->run(q` data[species] = mixing.ratio `);
}
$R->run(q` data = melt(data, id.vars = c("Time"), variable.name = "Species", value.name = "Mixing.Ratio") `);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Mixing.Ratio, colour = Species, group = Species)) `,
        q` plot = plot + geom_line() `,
);

$R->set('filename', "EMEP_${species}_mixing_ratio_comparison.pdf");
$R->run(q` CairoPDF(file = filename) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_all_tagged_species {
    my ($species) = @_;
    my @tagged;
    my $spc_file = "$base/gas.spc";
    open my $in, '<:encoding(utf-8)', $spc_file or die "Can't open $spc_file : $!";
    my @lines = <$in>;
    close $in;
    foreach my $line (@lines) {
        next unless ($line =~ /${species}_/);
        chomp $line;
        (my $tagged = $line) =~ s/^${species}_(.*?)\b.*$/${species}_$1/;
        push @tagged, $tagged;
    }
    return @tagged;
}
