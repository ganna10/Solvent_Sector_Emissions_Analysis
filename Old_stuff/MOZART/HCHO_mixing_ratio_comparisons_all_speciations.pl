#! /usr/bin/env perl
# Compare CH2O mixing ratios in all speciations
# Version 0: Jane Coates 24/11/2014

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions/MOZART";
opendir DIR, $base or die "Can't open $base : $!";
my @runs = grep { $_ =~ /_Solvents_Only$/ } readdir DIR;
closedir DIR;

my $mecca = MECCA->new("$base/DE94_Solvents_Only/boxmodel");
my $ntime = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times = $times(1:$ntime-2);

my %plot_data;
foreach my $run (@runs) {
    (my $label = $run) =~ s/^(.*?)_Solvents_Only/$1/;
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $conc = $mecca->tracer("CH2O");
    $conc = $conc(1:$ntime-2) * 1e9;
    $plot_data{$label} += $conc;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
);

$R->set('Time', [map { $_ } $times->dog]);
$R->run(q` data = data.frame(Time) `);
foreach my $speciation (keys %plot_data) {
    $R->set('speciation', $speciation);
    $R->set('mixing.ratio', [map { $_ } $plot_data{$speciation}->dog]);
    $R->run(q` data[speciation] = mixing.ratio `);
}
$R->run(q` data = melt(data, id.vars = c("Time"), variable.name = "Speciation", value.name = "Mixing.Ratio") `);
$R->run(q` my.colours = c("#010101", "#f37d22", "#008c47", "#1859a9", "#ed2d2e", "#662c91", "#12b2b2", "#b33893", "#a11d20") `);
$R->run(q` plot = ggplot(data, aes(x = Time, y = Mixing.Ratio, colour = Speciation, group = Speciation)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "CH2O_mixing_ratio_comparison.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
); 

$R->stop();
