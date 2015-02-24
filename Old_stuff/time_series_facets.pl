#! /usr/bin/env perl
# Mixing Ratio time series plots, facet by solvent speciation. ARGV is species name
# Version 0 : Jane Coates 19/1/2015

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use Statistics::R;
use MECCA;

my $base = "/local/home/coates/Solvent_Emissions";
#Create x-axis for plot in hours
my $mecca = MECCA->new("$base/MCM/TNO_Solvents_Only/boxmodel"); 
my $ntime = $mecca->time->nelem; #number of time points
my $times = $mecca->time;
$times -= $times->at(0);
$times = $times(1:$ntime-2);
$times /= 86400;

my %data;
my $species = $ARGV[0];
my @mechanisms = qw( MCM MOZART RADM2 );
my @speciations = qw( DE94 EMEP GR05 GR95 IPCC TNO UK08 UK98 );

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $boxmodel = "$base/$mechanism/${speciation}_Solvents_Only/boxmodel";
        my $mecca = MECCA->new($boxmodel);
        my $tracer = $mecca->tracer($species);
        $data{$mechanism}{$speciation} = $tracer(1:$ntime-2);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(grid) `,
        q` library(Cairo) `,
);

$R->set('Time', [map { $_ } $times->dog]);
$R->set('title', "$species Mixing Ratios");
$R->set('file.name', "${species}_mixing_ratio_comparison.pdf");
$R->run(q` data = data.frame() `);

foreach my $mechanism (sort keys %data) {
    $R->run(q` pre = data.frame(Time) `);
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$data{$mechanism}} ) {
        $R->set('speciation', $speciation);
        $R->set('mixing.ratio', [ map { $_ } $data{$mechanism}{$speciation}->dog ]);
        $R->run(q` pre[speciation] = mixing.ratio * 1e9 `);
    }
    $R->run(q` pre$Mechanism = rep(mechanism, length(time)) `,
            q` pre = gather(pre, Speciation, Mixing.Ratio, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
} 
#my $p = $R->run(q` print(pre) `);
#print "$p\n";

$R->run(q` my.colours = c("MCM" = "#6c254f", "MOZART" = "#ef6638", "RADM2" = "#0e5c28") `);
$R->run(q` plot = ggplot(data, aes(x = Time, y = Mixing.Ratio, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + facet_wrap( ~ Speciation) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + ylab("Mixing Ratio (ppbv)") `,
        q` plot = plot + xlab("Time (days)") `,
        q` plot = plot + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), expand = c(0, 0)) `,
        q` plot = plot + ggtitle(title) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(panel.margin.x = unit(0.3, "cm")) `,
        q` plot = plot + theme(panel.border = element_rect(colour  = "black")) `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.key = element_blank()) `,
        q` plot = plot + theme(legend.position = c(0.85, 0.1)) `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = file.name) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
