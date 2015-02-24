#! /usr/bin/env perl
# Compare MEK mixing ratios as emitted species and degradation product in EMEP Solvents Only run
# Version 0: Jane Coates 25/11/2014

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use Statistics::R;

my $dir = "/local/home/coates/Solvent_Emissions/MCM/EMEP_tagged_solvents_only_all";
my $mecca = MECCA->new("$dir/boxmodel");
my $ntime = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times = $times(1:$ntime-2);
my (%data, @species);

my $spc_file = "$dir/gas.spc";
open my $in, '<:encoding(utf-8)', $spc_file or die "Can't open $spc_file : $!";
my @lines = <$in>;
close $in;

foreach my $line (@lines) {
    next unless ($line =~ /^MEK\s|^MEK_/);
    (my $species = $line) =~ s/^(.*?)\s=.*\n$/$1/;
    push @species, $species;
}

foreach my $species (@species) {
    my $conc = $mecca->tracer($species);
    $conc = $conc(1:$ntime-2) * 1e9;
    $data{$species} = $conc;
    $data{"Total"} += $conc;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->set('Time', [ map { $_ } $times->dog ]);
$R->run(q` data = data.frame(Time) `);

foreach my $species (sort keys %data) {
    $R->set('species', $species);
    $R->set('mixing.ratio', [ map { $_ } $data{$species}->dog ]);
    $R->run(q` data[species] = mixing.ratio `);
}
$R->run(q` data = gather(data, Species, Mixing.Ratio, -Time) `);
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` plot = ggplot(data, aes(x = Time, y = Mixing.Ratio, colour = Species, group = Species)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 1)) `,
);

$R->run(q` CairoPDF(file = "MEK_tagged_mixing_ratio_comparison.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
