#! /usr/bin/env perl
# standard deviations of each mechanism runs between the speciations
# Version 0: Jane Coates 25/6/2015

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $species = "O3";

my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw(MCM MOZART RADM2);
my @speciations = qw( DE94 EMEP GR05 GR95 IPCC TNO UK08 UK98 );
my @runs = qw( Solvents_Only mean_NO_source_Solvents_Only );
my %data;

my $mecca = MECCA->new("$base/RADM2/DE94_Solvents_Only/boxmodel");
my $ntime = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $n_per_day = 86400 / $dt;
my $n_days = int $ntime / $n_per_day;

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        foreach my $run (@runs) {
            my $dir = "$base/$mechanism/${speciation}_$run";
            my $mecca = MECCA->new("$dir/boxmodel");
            my $mr = $mecca->tracer($species); 
            $mr = $mr->reshape($n_per_day, $n_days);
            $mr = $mr->sumover;
            $data{$run}{$mechanism}{$speciation} = $mr;
        }
    }
}

my $R = Statistics::R->new();
$R->run(q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(ggplot2) `,
        q` library(Cairo) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` d = data.frame() `);
foreach my $run (sort keys %data) {
    $R->set('run', $run);
    foreach my $mechanism (sort keys %{$data{$run}}) {
        $R->set('mechanism', $mechanism);
        $R->run(q` pre = data.frame(Time, Run = rep(run, length(Time)), Mechanism = rep(mechanism, length(Time))) `);
        foreach my $speciation (sort keys %{$data{$run}{$mechanism}}) {
            $R->set('speciation', $speciation);
            $R->set('mixing.ratio', [  map { $_ } $data{$run}{$mechanism}{$speciation}->dog ]);
            $R->run(q` pre[speciation] = mixing.ratio `);
        }
        $R->run(q` pre = gather(pre, Speciation, Mixing.Ratio, -Time, -Run, -Mechanism) `,
                q` d = rbind(d, pre) `,
        );
    }
}

$R->run(q` mean.NO = d %>% filter(Run == "mean_NO_source_Solvents_Only") %>% select(-Run) `,
        q` mcm.mean.NO = mean.NO %>% filter(Mechanism == "MCM") %>% select(-Mechanism) `,
        q` sd.mcm.mean.NO = mcm.mean.NO %>% group_by(Time) %>% summarise(St.Dev = sd(Mixing.Ratio)) `,
        q` mozart.mean.NO = mean.NO %>% filter(Mechanism == "MOZART") %>% select(-Mechanism) `,
        q` sd.mozart.mean.NO = mozart.mean.NO %>% group_by(Time) %>% summarise(St.Dev = sd(Mixing.Ratio)) `,
        q` radm2.mean.NO = mean.NO %>% filter(Mechanism == "RADM2") %>% select(-Mechanism) `,
        q` sd.radm2.mean.NO = radm2.mean.NO %>% group_by(Time) %>% summarise(St.Dev = sd(Mixing.Ratio)) `,
);

$R->run(q` solvents = d %>% filter(Run == "Solvents_Only") %>% select(-Run) `,
        q` mcm.solvents = solvents %>% filter(Mechanism == "MCM") %>% select(-Mechanism) `,
        q` sd.mcm.solvents = mcm.solvents %>% group_by(Time) %>% summarise(St.Dev = sd(Mixing.Ratio)) `,
        q` mozart.solvents = solvents %>% filter(Mechanism == "MOZART") %>% select(-Mechanism) `,
        q` sd.mozart.solvents = mozart.solvents %>% group_by(Time) %>% summarise(St.Dev = sd(Mixing.Ratio)) `,
        q` radm2.solvents = solvents %>% filter(Mechanism == "RADM2") %>% select(-Mechanism) `,
        q` sd.radm2.solvents = radm2.solvents %>% group_by(Time) %>% summarise(St.Dev = sd(Mixing.Ratio)) `,
);

open my $out, '>:encoding(utf-8)', "standard_deviations.txt" or die $!;
print $out "Solvents Only Runs\n";
print $out "MCM standard deviation of each day\n";
my $p = $R->run(q` print(sd.mcm.solvents) `);
print $out $p, "\n";
print $out "\nMOZART standard deviation of each day\n";
$p = $R->run(q` print(sd.mozart.solvents) `);
print $out $p, "\n";
print $out "\nRADM2 standard deviation of each day\n";
$p = $R->run(q` print(sd.radm2.solvents) `);
print $out $p, "\n";
print $out "\n\nMean NO Source Solvents Only Runs\n";
print $out "MCM standard deviation of each day\n";
$p = $R->run(q` print(sd.mcm.mean.NO) `);
print $out $p, "\n";
print $out "\nMOZART standard deviation of each day\n";
$p = $R->run(q` print(sd.mozart.mean.NO) `);
print $out $p, "\n";
print $out "\nRADM2 standard deviation of each day\n";
$p = $R->run(q` print(sd.radm2.mean.NO) `);
print $out $p, "\n";
close $out;

$R->run(q` d = d %>% group_by(Time, Mechanism, Run) %>% summarise(St.Dev = sd(Mixing.Ratio)) `);
$R->run(q` p = ggplot(d, aes(x = Time, y = St.Dev, colour = Mechanism)) `,
        q` p = p + geom_point() `,
        q` p = p + facet_wrap( ~ Run) `,
);

$R->run(q` CairoPDF(file = "std_devs.pdf", width = 10, height = 7) `,
        q` print(p) `,
        q` dev.off() `,
);

$R->stop();
