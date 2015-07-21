#! /usr/bin/env perl
# plot time series of argument species mixing ratios in each specification, each mechanism and each model run condition
# Version 0: Jane Coates 21/7/2015

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $species = "O3";

my $base = "/work/users/jco/Solvent_Emissions";
#my $base = "/local/home/coates/Solvent_Emissions";
#my @mechanisms = qw(MOZART );
my @mechanisms = qw(MCM MOZART RADM2);
my @speciations = qw( DE94 EMEP GR05 GR95 IPCC TNO UK08 UK98 );
#my @runs = qw( Solvents_Only all_sectors mean_NO_source_all_sectors mean_NO_source_Solvents_only );
my @runs = qw( Solvents_Only mean_NO_source_Solvents_Only );
my %data;

my $mecca = MECCA->new("$base/RADM2/DE94_Solvents_Only/boxmodel");
my $ntime = $mecca->time->nelem;
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;
$times =  $times(1:$ntime-2);

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        foreach my $run (@runs) {
            my $dir = "$base/$mechanism/${speciation}_$run";
            my $mecca = MECCA->new("$dir/boxmodel");
            my $mr = $mecca->tracer($species); 
            $data{$run}{$mechanism}{$speciation} = $mr(1:$ntime-2) * 1e9;
        }
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
        q` library(scales) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(grid) `,
);

$R->set('Time', [ map { $_ } $times->dog ]);
$R->run(q` my.colours = c("TNO" = "#000000", "IPCC" = "#2b9eb3", "EMEP" = "#ef6638", "DE94" = "#f9c500", "GR95" = "#6c254f", "GR05" = "#0352cb", "UK98" = "#0e5c28", "UK08" = "#b569b3") `);
$R->run(q` data = data.frame() `);
foreach my $run (sort keys %data) {
    $R->set('run', $run);
    foreach my $mechanism (sort keys %{$data{$run}}) {
        $R->run(q` pre = data.frame(Time, Run = rep(run, length(Time))) `);
        foreach my $speciation (sort keys %{$data{$run}{$mechanism}}) {
            $R->set('speciation', $speciation);
            $R->set('mixing.ratio', [ map { $_ } $data{$run}{$mechanism}{$speciation}->dog ]);
            $R->run(q` pre[speciation] = mixing.ratio `);
        }
        $R->run(q` pre = gather(pre, Speciation, Mixing.Ratio, -Time, -Run) `);
        $R->set('mechanism', $mechanism) ;
        $R->run(q` pre$Mechanism = rep(mechanism, length(pre$Mixing.Ratio)) `,
                q` data = rbind(data, pre) `,
        );
    }
}

$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` data = mutate(data, mechanism = factor(Mechanism, labels = c("MCM v3.2", "MOZART-4", "RADM2"))) `);
$R->run(q` data = mutate(data, run = factor(Run, labels = c("Solvents Only", "Mean NO Source"))) `);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Mixing.Ratio, colour = Speciation, group = Speciation)) `,
        q` plot = plot + facet_grid(run ~ mechanism) `,
        q` plot = plot + geom_vline(xintercept = 0:7, colour = "grey") `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + geom_line() `,
        q` plot = plot + scale_x_continuous(breaks = seq(0, 7, 0.5), expand = c(0, 0.01), labels = c("06:00", "18:00", "06:00", "18:00", "06:00", "18:00", "06:00", "18:00", "06:00", "18:00", "06:00", "18:00", "06:00", "18:00", "06:00")) `,
        q` plot = plot + ylab("Mixing Ratio (ppbv)") `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(strip.text.y = element_text(angle = 0)) `,
        q` plot = plot + theme(plot.title = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(panel.margin = unit(5, "mm")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
        q` plot = plot + theme(legend.position = "top") `,
);

$R->run(q` CairoPDF(file = "O3_mixing_ratio_time_series_Solvents_Only_mean_NO_source_and_tuned_NO_source.pdf", width = 11, height = 5) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
