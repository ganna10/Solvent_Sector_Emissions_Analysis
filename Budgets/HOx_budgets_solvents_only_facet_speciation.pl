#! /usr/bin/env perl
# Plot cumulative HOx production budgets, facet by speciation 
# Version 0: Jane Coates 24/3/2015

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );

my $mecca = MECCA->new("$base/RADM2/TNO_Solvents_Only/boxmodel");
my $ntime = $mecca->time->nelem;
my (%families, %weights, %data);
$families{"HOx"} = [ qw( OH HO2 HO2NO2 ) ];

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $boxmodel = "$base/$mechanism/${speciation}_Solvents_Only/boxmodel";
        my $mecca = MECCA->new($boxmodel);
        my $eqn = "$base/$mechanism/${speciation}_Solvents_Only/gas.eqn";
        my $kpp = KPP->new($eqn);

        my (%production_rates, %consumption_rates);
        $kpp->family({
                name    => "HOx",
                members => $families{"HOx"},
                weights => $weights{"HOx"},
        });
        my $producers = $kpp->producing("HOx");
        my $producer_yields = $kpp->effect_on("HOx", $producers);
        my $consumers = $kpp->consuming("HOx");
        my $consumer_yields = $kpp->effect_on("HOx", $consumers);
        print "No producers for HOx in $mechanism and $speciation\n" if (@$producers == 0);
        print "No consumers for HOx in $mechanism and $speciation\n" if (@$consumers == 0);

        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            next if ($rate->sum == 0);
            my $reaction_string = $kpp->reaction_string($reaction);
            my ($reactants, $products) = split / = /, $reaction_string;
            $production_rates{$reactants} += $rate(1:$ntime-2);
        }

        for (0..$#$consumers) {
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            next if ($rate->sum == 0);
            my $reaction_string = $kpp->reaction_string($reaction);
            my ($reactants, $products) = split / = /, $reaction_string;
            $consumption_rates{$reactants} += $rate(1:$ntime-2);
        }
        remove_common_processes(\%production_rates, \%consumption_rates);
        my $others = 2e8;
        foreach my $reaction (keys %production_rates) {
            my $sum = $production_rates{$reaction}->sumover;
            if ($sum < $others) {
                $production_rates{"Others"} += $sum;
                delete $production_rates{$reaction};
            } else {
                $production_rates{$reaction} = $sum;
            }
        }
        $data{$mechanism}{$speciation} = \%production_rates;
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
        q` library(grid) `,
);

$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        $R->run(q` pre = data.frame(Mechanism = mechanism) `);
        $R->set('speciation', $speciation);
        $R->run(q` pre$Speciation = speciation `);
        foreach my $process (sort keys %{$data{$mechanism}{$speciation}}) {
            $R->set('process', $process);
            $R->set('HOx.production', $data{$mechanism}{$speciation}{$process}->at(0));
            $R->run(q` pre[process] = HOx.production `);
        }
        $R->run(q` pre = gather(pre, Process, HOx.Production, -Mechanism, -Speciation) `,
                q` data = rbind(data, pre) `,
        );
        #my $p = $R->run(q` print(pre) `);
        #print $p, "\n";
    }
}
$R->set('filename', "Cumulative_HOx_budget_facet_speciation_solvents_only.pdf");
$R->set('title', "Cumulative HOx Production Budget");
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = HOx.Production, fill = Process)) `,
        q` plot = plot + geom_bar(stat = "identity", position = "stack") `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2 ) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + ggtitle(title) `,
        q` plot = plot + ylab("HOx Production (molecules cm-3)") `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + theme(axis.ticks.x = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
        q` plot = plot + theme(plot.title = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.y = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 0.7, vjust = 1.0)) `,
        q` plot = plot + theme(panel.margin = unit(5, "mm")) `,
);

$R->run(q` CairoPDF(file = filename, width = 8.6, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub remove_common_processes {
    my ($production, $consumption) = @_;
    my %common_processes;
    $common_processes{$_} = 1 for (grep { defined $production->{$_} } keys %$consumption) ;

    foreach my $process (keys %common_processes) {
        my $net_effect = $production->{$process} + $consumption->{$process};
        #print $process, $net_effect->nelem, "\n";
        if ($net_effect->sum > 0) { #if net production remove consumption processes, retain net production
            if (which($net_effect < 0)->nelem > 0) { #which gets indices of non-0 values, nelem get nr of elements
                #print "which if $process $net_effect\n";
                $production->{$process} .= $net_effect;
                $production->{$process}->where($net_effect < 0) .= 0; #where returns output pdl corresponding to input pdl and meeting the condition
                $consumption->{$process} .= $net_effect;
                $consumption->{$process}->where($net_effect > 0) .= 0; #where returns output pdl corresponding to input pdl and meeting the condition
                next;
            } 
            $production->{$process} = $net_effect;
            delete $consumption->{$process};
        } else { #net consumption
            if (which($net_effect > 0)->nelem > 0) {
                $production->{$process} .= $net_effect;
                $production->{$process}->where($net_effect < 0) .= 0;
                $consumption->{$process} .= $net_effect;
                $consumption->{$process}->where($net_effect > 0) .= 0;
                next;
            }
            $consumption->{$process} = $net_effect;
            delete $production->{$process};
        }
    }
} 
