#! /usr/bin/env perl
# compare NO source between all speciationss
# Version 0: Jane Coates 23/1/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions/MOZART";
my $mecca = MECCA->new("$base/DE94_Solvents_Only/boxmodel");
my $NTIME = $mecca->time->nelem;
my $time = $mecca->time;
$time -= $time->at(0);
$time /= 86400;
$time = $time(1:$NTIME-2);

my @speciations = qw( DE94 EMEP IPCC GR05 GR95 TNO UK98 UK08 );
my %plot_data;
foreach my $speciation (@speciations) {
    my $boxmodel = "$base/${speciation}_Solvents_Only/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/${speciation}_Solvents_Only/gas.eqn";
    my $kpp = KPP->new($eqn_file); 
    my $reaction = $kpp->producing_from("NO", "UNITY");
    my $reaction_number = $kpp->reaction_number(@$reaction);
    my $rate = $mecca->rate($reaction_number);
    $plot_data{$speciation} = $rate(1:$NTIME-2);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->set('time', [map { $_ } $time->dog]);
$R->run(q` data = data.frame(time) `);
foreach my $run (sort keys %plot_data) {
    $R->set('speciations', $run);
    $R->set('rate', [map { $_ } $plot_data{$run}->dog]);
    $R->run(q` data[speciations] = rate `);
}

$R->run(q` data = gather(data, Speciation, Rate, -time) `);
$R->run(q` my.colours = c("#010101", "#f37d22", "#008c47", "#1859a9", "#ed2d2e", "#662c91", "#12b2b2", "#b33893", "#a11d20") `);

$R->run(q` plot = ggplot(data, aes(x = time, y = Rate, colour = Speciation)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "NO_source_all_speciations.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
