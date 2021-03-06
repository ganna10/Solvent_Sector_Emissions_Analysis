#! /usr/bin/env perl
# Mixing Ratio time series plots, facet by solvent speciation. ARGV is species name
# Version 0 : Jane Coates 19/1/2015
# Version 1: Jane Coates 25/2/2015 plotting aesthetic changes
# Version 2: Jane Coates 29/4/2015 plot changes
# Version 3: Jane Coates 13/5/2015 adding mean_NO_source runs to script

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
my @runs = qw( Solvents_Only );
#my @runs = qw( Solvents_Only mean_NO_source );

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        foreach my $run (@runs) {
            my $boxmodel = "$base/$mechanism/${speciation}_Solvents_Only/boxmodel";
            my $mecca = MECCA->new($boxmodel);
                my $lookup;
                if ($mechanism eq "MOZART" and $species eq "HCHO") {
                    $lookup = "CH2O";
                } else {
                    $lookup = $species;
                }
            my $tracer = $mecca->tracer($lookup);
            $data{$run}{$mechanism}{$speciation} = $tracer(1:$ntime-2);
        }
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
foreach my $run (keys %data) {
    $R->set('title', "$species Mixing Ratios in $run");
    $R->set('file.name', "${species}_mixing_ratio_comparison_by_speciation_${run}.pdf");
    $R->run(q` data = data.frame() `);

    foreach my $mechanism (sort keys %{$data{$run}}) {
        $R->run(q` pre = data.frame(Time) `);
        $R->set('mechanism', $mechanism);
        foreach my $speciation (sort keys %{$data{$run}{$mechanism}} ) {
            $R->set('speciation', $speciation);
            $R->set('mixing.ratio', [ map { $_ } $data{$run}{$mechanism}{$speciation}->dog ]);
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
            q` plot = plot + geom_vline(xintercept = 0:7, colour = "grey") `,
            q` plot = plot + geom_line() `,
            q` plot = plot + facet_wrap( ~ Speciation, nrow = 2) `,
            q` plot = plot + theme_tufte() `,
            #q` plot = plot + theme_hc() `,
            q` plot = plot + ylab("Mixing Ratio (ppbv)") `,
            q` plot = plot + scale_x_continuous(breaks = seq(0, 7, 0.5), expand = c(0, 0.01), labels = c("06:00", "18:00", "06:00", "18:00", "06:00", "18:00", "06:00", "18:00", "06:00", "18:00", "06:00", "18:00", "06:00", "18:00", "06:00")) `,
            q` plot = plot + scale_y_continuous(limits = c(35, 110)) `,
            q` plot = plot + theme(strip.background = element_blank()) `,
            q` plot = plot + theme(axis.text = element_text(colour = "black")) `,
            q` plot = plot + theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) `,
            q` plot = plot + theme(axis.ticks = element_line(colour = "black")) `,
            q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
            q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
            q` plot = plot + theme(panel.margin.x = unit(0.3, "cm")) `,
            q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
            q` plot = plot + theme(axis.title.x = element_blank()) `,
            q` plot = plot + theme(legend.title = element_blank()) `,
            q` plot = plot + theme(legend.position = "top") `,
            q` plot = plot + scale_colour_manual(values = my.colours) `,
            #q` plot = plot + theme(panel.grid.major.y = element_line(colour ="black")) `,
    );

    $R->run(q` CairoPDF(file = file.name, width = 11, height = 7) `,
            q` print(plot) `,
            q` dev.off() `,
    );
}

$R->stop();
