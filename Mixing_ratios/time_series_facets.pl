#! /usr/bin/env perl
# Mixing Ratio time series plots, facet by solvent speciation. ARGV is species name
# Version 0 : Jane Coates 19/1/2015
# Version 1: Jane Coates 25/2/2015 plotting aesthetic changes
# Version 2: Jane Coates 29/4/2015 plot changes

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use Statistics::R;
use MECCA;

my $species = $ARGV[0];
die "Need to specify species!" unless defined $species;

my $base = "/local/home/coates/Solvent_Emissions";
#Create x-axis for plot in hours
my $mecca = MECCA->new("$base/MCM/TNO_Solvents_Only/boxmodel"); 
my $ntime = $mecca->time->nelem; #number of time points
my $times = $mecca->time;
$times -= $times->at(0);
$times = $times(1:$ntime-2);
$times /= 86400;

my %data;
my @mechanisms = qw( MCM MOZART RADM2 );
my @speciations = qw( DE94 EMEP GR05 GR95 IPCC TNO UK08 UK98 );

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $boxmodel = "$base/$mechanism/${speciation}_Solvents_Only/boxmodel";
        my $mecca = MECCA->new($boxmodel);
            my $lookup;
            if ($mechanism eq "MOZART" and $species eq "HCHO") {
                $lookup = "CH2O";
            } else {
                $lookup = $species;
            }
        my $tracer = $mecca->tracer($lookup);
        $data{$mechanism}{$speciation} = $tracer(1:$ntime-2);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(grid) `,
        q` library(ggthemes) `,
        q` library(Cairo) `,
);

$R->set('Time', [map { $_ } $times->dog]);
$R->set('title', "$species Mixing Ratios");
$R->set('file.name', "${species}_mixing_ratio_comparison_by_speciation.pdf");
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
#my $p = $R->run(q` print(data) `);
#print "$p\n";

$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` my.colours = c("MCM" = "#6c254f", "MOZART" = "#ef6638", "RADM2" = "#0e5c28") `);
$R->run(q` plot = ggplot(data, aes(x = Time, y = Mixing.Ratio, colour = Mechanism, group = Mechanism)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + ylab("Mixing Ratio (ppbv)") `,
        q` plot = plot + xlab("Time (days)") `,
        q` plot = plot + scale_x_continuous(limits = c(0, 7), breaks = seq(0, 7, 1), expand = c(0, 0)) `,
        q` plot = plot + ggtitle(title) `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(plot.title = element_text(face = "bold")) `,
        q` plot = plot + theme(panel.margin.x = unit(0.3, "cm")) `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(legend.position = "top") `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = file.name, width = 10, height = 7) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
