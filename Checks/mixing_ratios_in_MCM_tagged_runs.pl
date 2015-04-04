#! /usr/bin/env perl
# Compare mixing ratios of a species (argument) in the different tagged model runs of the speciations
# Version 0: Jane Coates 31/3/2015

use strict;
use diagnostics;
use MECCA;
use PDL;
use Statistics::R;

my $species = $ARGV[0];
die "specify species\n" unless (defined $species);

my $base = "/local/home/coates/Solvent_Emissions/MCM";
my @speciations = qw( TNO IPCC DE94 GR95 GR05 UK98 UK08 );

my $mecca = MECCA->new("$base/GR05_tagged_solvents_only_alkanes/boxmodel");
my $times = $mecca->time;

my %data;
foreach my $speciation (@speciations) {
    opendir DIR, $base or die "Can't open $base : $!";
    my @dirs = grep { $_ =~ /${speciation}_tagged_solvents_only/ } readdir DIR;
    close DIR;
    foreach my $dir (@dirs) {
        my $mecca = MECCA->new("$base/$dir/boxmodel");
        (my $label = $dir) =~ s/${speciation}_tagged_solvents_only_(.*?)/$1/; 
        $data{$speciation}{$label} = $mecca->tracer($species);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);

$R->set('Time', [ map { $_ } $times->dog ]);
$R->run(q` data = data.frame() `);
foreach my $speciation (sort keys %data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $run (sort keys %{$data{$speciation}}) {
        $R->set('run', $run);
        $R->set('mr', [ map { $_ } $data{$speciation}{$run}->dog ]);
        $R->run(q` pre[run] = mr `);
    }
    $R->set('speciation', $speciation);
    $R->run(q` pre$Speciation = rep(speciation, length(Time)) `,
            q` pre = gather(pre, Run, MR, -Time, -Speciation) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` p = ggplot(data, aes(x = Time, y = MR, colour = Run, group = Run)) `,
        q` p = p + geom_line() `,
        q` p = p + facet_wrap( ~ Speciation ) `,
);

$R->set('name', "${species}_tagged_mixing_ratio_comparison.pdf");
$R->run(q` CairoPDF(file = name) `,
        q` print(p) `,
        q` dev.off() `,
);

$R->stop();
