#! /usr/bin/env perl
# Cumulative Ox production budgets, faceted by mechanism full Ox family
# Version 0: Jane Coates 24/3/2015
# Version 1: Jane Coates 27/3/2015 including both facetting in script
# Version 2: Jane Coates 29/4/2015 removing HO2x from analysis 

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions";
#my $base = "/work/users/jco/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my (%families, %weights, %data);
$families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];

my $mecca = MECCA->new("$base/RADM2/TNO_Solvents_Only/boxmodel");
my $ntime = $mecca->time->nelem;

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my (%production_rates, %consumption_rates);
        my $boxmodel = "$base/$mechanism/${speciation}_Solvents_Only/boxmodel";
        my $mecca = MECCA->new($boxmodel);
        my $eqn_file = "$base/$mechanism/${speciation}_Solvents_Only/gas.eqn";
        my $kpp = KPP->new($eqn_file); 
        my $ro2_file = "$base/$mechanism/${speciation}_Solvents_Only/RO2_species.txt";
        my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
        $families{"Ox"} = [ qw( O3 NO2 O NO3 N2O5 HO2NO2 O1D ), @no2_reservoirs ];
        $weights{"Ox"} = { NO3 => 2, N2O5 => 3 };

        foreach my $species (qw( Ox )) {
            $kpp->family({
                    name    => $species,
                    members => $families{$species},
                    weights => $weights{$species},
            });
            my $producers = $kpp->producing($species);
            my $producer_yields = $kpp->effect_on($species, $producers);
            my $consumers = $kpp->consuming($species);
            my $consumer_yields = $kpp->effect_on($species, $consumers);

            print "No producers for $species in $mechanism, $speciation\n" if (@$producers == 0);
            print "No consumers for $species in $mechanism, $speciation\n" if (@$consumers == 0);

            for (0..$#$producers) {
                my $reaction = $producers->[$_];
                my $reaction_number = $kpp->reaction_number($reaction);
                my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
                next if ($rate->sum == 0);
                my ($number, $parent) = split /_/, $reaction;
                if (defined $parent) {
                    $production_rates{$species}{$parent} += $rate(1:$ntime-2);
                } else {
                    my $reaction_string = $kpp->reaction_string($reaction);
                    my ($reactants, $products) = split / = /, $reaction_string;
                    $production_rates{$species}{$reactants} += $rate(1:$ntime-2);
                }
            }

            for (0..$#$consumers) {
                my $reaction = $consumers->[$_];
                my $reaction_number = $kpp->reaction_number($reaction);
                my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
                next if ($rate->sum == 0);
                my ($number, $parent) = split /_/, $reaction;
                if (defined $parent) {
                    $consumption_rates{$species}{$parent} += $rate(1:$ntime-2);
                } else {
                    my $reaction_string = $kpp->reaction_string($reaction);
                    my ($reactants, $products) = split / = /, $reaction_string;
                    $consumption_rates{$species}{$reactants} += $rate(1:$ntime-2);
                }
            }
        }

#        remove_common_processes($production_rates{"HO2x"}, $consumption_rates{"HO2x"});
#        my $total_ho2x_production;
#        $total_ho2x_production += $production_rates{"HO2x"}{$_} foreach (keys %{$production_rates{"HO2x"}});
#        
#        foreach my $reaction (keys %{$production_rates{'HO2x'}}) {
#            $production_rates{"Ox"}{$reaction} += $production_rates{"Ox"}{'HO2 + NO'} * $production_rates{'HO2x'}{$reaction} / $total_ho2x_production;
#            $consumption_rates{"Ox"}{$reaction} += $consumption_rates{"Ox"}{'HO2 + O3'} * $consumption_rates{'HO2x'}{$reaction} / $total_ho2x_production; 
#        }
#        delete $production_rates{"Ox"}{'HO2 + NO'};
#        delete $consumption_rates{"Ox"}{'HO2 + O3'};
        remove_common_processes($production_rates{"Ox"}, $consumption_rates{"Ox"});

        my $others = 4e8;
        foreach my $process (keys %{$production_rates{"Ox"}}) {
            if ($production_rates{"Ox"}{$process}->sum < $others) {
                $production_rates{"Ox"}{"Others"} += $production_rates{"Ox"}{$process}->sum;
                delete $production_rates{"Ox"}{$process};
            } else {
                $production_rates{"Ox"}{$process} = $production_rates{"Ox"}{$process}->sum;
            }
        }
        
        my $sort_function = sub { $_[0] };
        my @sorted_prod = sort { &$sort_function($production_rates{"Ox"}{$b}) <=> &$sort_function($production_rates{"Ox"}{$a}) } keys %{$production_rates{"Ox"}};
        my @final_sorted;
        foreach (@sorted_prod) {
            next if ($_ =~ /Others/);
            push @final_sorted, { $_ => $production_rates{"Ox"}{$_} };
        }
        push @final_sorted, { "Others" => $production_rates{"Ox"}{"Others"} } if (defined $production_rates{"Ox"}{"Others"});
        $data{$speciation}{$mechanism} = \@final_sorted;
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
foreach my $speciation (sort keys %data) {
    $R->set('speciation', $speciation);
    foreach my $mechanism (sort keys %{$data{$speciation}}) {
        $R->run(q` pre = data.frame(Speciation = speciation) `);
        $R->set('mechanism', $mechanism);
        $R->run(q` pre$Mechanism = mechanism `);
        foreach my $ref (@{$data{$speciation}{$mechanism}}) {
            foreach my $process (sort keys %$ref) {
                $R->set('process', $process);
                $R->set('Ox.production', $ref->{$process});
                $R->run(q` pre[process] = Ox.production `);
            }
        }
        $R->run(q` pre = gather(pre, Process, Ox.Production, -Mechanism, -Speciation) `,
                q` data = rbind(data, pre) `,
        );
        #my $p = $R->run(q` print(pre) `);
        #print $p, "\n";
    }
}

$R->run(q` plot.lines = function () { list( theme_tufte() ,
                                            ggtitle("Cumulative Ox Production Budget") ,
                                            ylab("Ox Production (molecules cm-3)") ,
                                            scale_y_continuous(expand = c(0, 0)) ,
                                            scale_x_discrete(expand = c(0, 0)) ,
                                            theme(axis.ticks.x = element_blank()) ,
                                            theme(legend.title = element_blank()) ,
                                            theme(axis.line = element_line(colour = "black")) ,
                                            theme(plot.title = element_text(face = "bold")) ,
                                            theme(axis.title.y = element_text(face = "bold")) ,
                                            theme(axis.title.x = element_blank()) ,
                                            theme(strip.text = element_text(face = "bold")) ,
                                            theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 0.7, vjust = 1.0)) ,
                                            theme(panel.margin = unit(5, "mm")) ) } `);

$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);

$R->run(q` plot = ggplot(data, aes(x = Speciation, y = Ox.Production, fill = Process)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism, nrow = 1 ) `,
        q` plot = plot + plot.lines() `,
);

$R->run(q` CairoPDF(file = "Cumulative_Ox_budget_facet_mechanism.pdf", width = 8.6, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = Ox.Production, fill = Process)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2 ) `,
        q` plot = plot + plot.lines() `,
);

$R->run(q` CairoPDF(file = "Cumulative_Ox_budget_facet_speciation.pdf", width = 8.6, height = 6) `,
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

sub get_no2_reservoirs { #get species that are produced when radical species react with NO2
    my ($kpp, $file) = @_; 
    open my $in, '<:encoding(utf-8)', $file or die $!; 
    my @ro2;
    for (<$in>) {
        push @ro2, split /\s+/, $_; 
    }
    close $in;
    my @no2_reservoirs;
    foreach my $ro2 (@ro2) {
        my ($reactions) = $kpp->reacting_with($ro2, 'NO2');
        foreach my $reaction (@$reactions) {
            my ($products) = $kpp->products($reaction);
            if (@$products == 1) {
                push @no2_reservoirs, $products->[0];
            }   
        }   
    }   
    return @no2_reservoirs;
} 
