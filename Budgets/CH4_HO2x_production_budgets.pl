#! /usr/bin/env perl
# Plot cumulative HO2x budgets from CH4 degradation-> not Ox as Ox production from CH4 is from HO2x, minor contributions in MCM from CH3NO3
# Version 0: Jane Coates 6/5/2015

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use MECCA;
use KPP;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw( MOZART RADM2 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my (%families, %weights, %data);
$families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        if ($mechanism eq "MCM") {
        } else {
            my $dir = "$base/$mechanism/${speciation}_tagged_solvents_only";
            my $mecca = MECCA->new("$dir/boxmodel");
            my $kpp = KPP->new("$dir/gas.eqn");
            $data{$mechanism}{$speciation} = get_data($mecca, $kpp, $mechanism, $speciation);
        }
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(grid) `,
        q` library(scales) `,
        q` library(ggthemes) `,
);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->run(q` pre = data.frame(Mechanism = mechanism, Speciation = speciation) `);
        foreach my $reactant (sort keys %{$data{$mechanism}{$speciation}}) {
            $R->set('reactants', $reactant);
            $R->set('ho2x', $data{$mechanism}{$speciation}{$reactant});
            $R->run(q` pre[reactants] = ho2x `);
        }
        $R->run(q` pre = gather(pre, Reactants, HO2x, -Mechanism, -Speciation) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = HO2x, fill = Reactants)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(axis.text.x = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
);

$R->run(q` CairoPDF(file = "CH4_HO2x_production_reactions_facet_speciation.pdf", width = 10, height = 7) `,
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

sub get_data {
    my ($mecca, $kpp, $mechanism, $speciation, $no_dirs) = @_;
    $no_dirs = 1 unless (defined $no_dirs);

    my (%production_rates, %consumption_rates);
    foreach my $species (qw( HO2x )) {
        $kpp->family({
                name    => $species,
                members => $families{$species},
                weights => $weights{$species},
            });
        my $producers = $kpp->producing($species);
        my $producer_yields = $kpp->effect_on($species, $producers);
        my $consumers = $kpp->consuming($species);
        my $consumer_yields = $kpp->effect_on($species, $consumers);
        print "No consumers for $species\n" if (@$consumers == 0);
        print "No producers for $species\n" if (@$producers == 0);

        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            my ($number, $parent) = split /_/, $reaction;
            next unless (defined $parent and $parent eq "CH4");
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            next if ($rate->sum == 0);
            my ($reactants, $products) = split / = /, $kpp->reaction_string($reaction);
            $reactants =~ s/_$parent//g;
            $reactants =~ s/MO2/CH3O2/;
            $reactants =~ s/CH2O/HCHO/g;
            $reactants =~ s/OP1/HCOOH/g;
            $production_rates{$species}{$reactants} += $rate;
        }

        for (0..$#$consumers) {
            my $reaction = $consumers->[$_];
            my ($number, $parent) = split /_/, $reaction;
            next unless (defined $parent and $parent eq "CH4");
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            next if ($rate->sum == 0);
            my ($reactants, $products) = split / = /, $kpp->reaction_string($reaction);
            $reactants =~ s/_$parent//g;
            $reactants =~ s/MO2/CH3O2/;
            $reactants =~ s/CH2O/HCHO/g;
            $reactants =~ s/OP1/HCOOH/g;
            $consumption_rates{$species}{$reactants} += $rate;
        }
    }
    remove_common_processes($production_rates{"HO2x"}, $consumption_rates{"HO2x"}); 
    $production_rates{"HO2x"}{$_} = $production_rates{"HO2x"}{$_}->sum foreach (keys %{$production_rates{"HO2x"}});
    return $production_rates{"HO2x"};
}
