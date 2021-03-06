#! /usr/bin/env perl
# HCHO production budgets in Solvents Only runs and plot both facetting variables separately
# Version 0: Jane Coates 25/3/2015

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/work/users/jco/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my %data;

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $boxmodel = "$base/$mechanism/${speciation}_Solvents_Only/boxmodel";
        my $mecca = MECCA->new($boxmodel);
        my $eqn = "$base/$mechanism/${speciation}_Solvents_Only/gas.eqn";
        my $kpp = KPP->new($eqn);
        my $ntime = $mecca->time->nelem;

        my (%production_rates, %consumption_rates);
        my $lookup;
        if ($mechanism eq "MOZART") {
            $lookup = "CH2O";
        } else {
            $lookup = "HCHO";
        }
        my $producers = $kpp->producing($lookup);
        my $producer_yields = $kpp->effect_on($lookup, $producers);
        my $consumers = $kpp->consuming($lookup);
        my $consumer_yields = $kpp->effect_on($lookup, $consumers);
        print "No producers for $lookup in $mechanism and $speciation\n" if (@$producers == 0);
        print "No producers for $lookup in $mechanism and $speciation\n" if (@$producers == 0);

        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            next if ($rate->sum == 0);
            my $reaction_string = $kpp->reaction_string($reaction);
            my ($reactants, $producers) = split / = /, $reaction_string;
            $production_rates{$reactants} += $rate(1:$ntime-2);
        }

        for (0..$#$consumers) {
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            next if ($rate->sum == 0);
            my $reaction_string = $kpp->reaction_string($reaction);
            my ($reactants, $producers) = split / = /, $reaction_string;
            $consumption_rates{$reactants} += $rate(1:$ntime-2);
        }
        remove_common_processes(\%production_rates, \%consumption_rates); 
        my $others = 1e8;
        foreach my $reaction (sort keys %production_rates) {
            if ($production_rates{$reaction}->sum < $others) {
                $production_rates{"Others"} += $production_rates{$reaction}->sum;
                delete $production_rates{$reaction};
            } else {
                $production_rates{$reaction} = $production_rates{$reaction}->sum;
            }
        }
        my $sort_function = sub { $_[0] };
        my @sorted_prod = sort { &$sort_function($production_rates{$b}) <=> &$sort_function($production_rates{$a}) } keys %production_rates;
        my @sorted;
        foreach (@sorted_prod) {
            next if ($_ =~ /Others/);
            push @sorted, { $_ => $production_rates{$_} };
        }
        push @sorted, { "Others" => $production_rates{"Others"} } if (defined $production_rates{"Others"});
        $data{$mechanism}{$speciation} = \@sorted;
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
        foreach my $ref (@{$data{$mechanism}{$speciation}}) {
            foreach my $process (sort keys %$ref) {
                $R->set('process', $process);
                $R->set('HCHO.production', $ref->{$process});
                $R->run(q` pre[process] = HCHO.production `);
            }
        }
#my $p = $R->run(q` print(pre) `);
#print $p, "\n";
        $R->run(q` pre = gather(pre, Process, HCHO.Production, -Mechanism, -Speciation) `,
                q` data = rbind(data, pre) `,
        );
    }
}

$R->run(q` plot.lines = function () { list( theme_tufte() ,
                                            ggtitle("Cumulative HCHO Production Budget") ,
                                            ylab("HCHO Production (molecules cm-3)") ,
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
$R->run(q` plot = ggplot(data, aes(x = Speciation, y = HCHO.Production, fill = Process)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism, nrow = 1 ) `,
        q` plot = plot + plot.lines() `,
);

$R->run(q` CairoPDF(file = "Cumulative_HCHO_budget_full_HCHO_facet_mechanism_solvents_only.pdf", width = 8.6, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = HCHO.Production, fill = Process)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2 ) `,
        q` plot = plot + plot.lines() `,
);

$R->run(q` CairoPDF(file = "Cumulative_HCHO_budget_full_HCHO_facet_speciation_solvents_only.pdf", width = 8.6, height = 6) `,
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
