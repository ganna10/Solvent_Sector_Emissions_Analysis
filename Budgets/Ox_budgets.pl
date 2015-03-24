#! /usr/bin/env perl
# Plot Ox production and consumption budgets for each mechanism and speciation, conditions
# Version 0: Jane Coates 26/2/2015
# Version 1: Jane Coates 3/3/2015 plotting solvents only Ox production and all mechanisms on one plot by facettingnd cumulative Ox production over whole model run
#
use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my (%families, %weights, %data);

my $mecca = MECCA->new("$base/RADM2/TNO_Solvents_Only/boxmodel");
my $ntime = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $n_per_day = 43200 / $dt;
my $n_days = int $ntime / $n_per_day;

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $boxmodel = "$base/$mechanism/${speciation}_Solvents_Only/boxmodel";
        my $mecca = MECCA->new($boxmodel);
        my $eqn_file = "$base/$mechanism/${speciation}_Solvents_Only/gas.eqn";
        my $kpp = KPP->new($eqn_file);
        my $ro2_file = "$base/$mechanism/${speciation}_Solvents_Only/RO2_species.txt";
        my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
        $families{"Ox_${mechanism}_$speciation"} = [ qw( O3 O1D O HO2NO2 NO2 NO3 N2O5 ), @no2_reservoirs ];
        $weights{"Ox_${mechanism}_$speciation"} = { NO3 => 2, N2O5 => 3 };
        $data{$mechanism}{$speciation} = get_data($mecca, $kpp, $mechanism, $speciation);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
);

$R->set('Time', [ ("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7") ]);
$R->run(q` my.colours = c(  "Production Others" = "#623812",
                            "Consumption Others" = "#0352cb",
                            "O3" = "#f9c500",
                            "HO2 + NO" = "#6c254f",
                            "O1D" = "#ae4901",
                            "NO2 + OH" = "#0e5c28",
                            "HC3P + NO" = "#86b650",
                            "CH3O2 + NO" = "#2b9eb3", "MO2 + NO" = "#2b9eb3", 
                            "CH3CO3 + NO" = "#ef6638", "ACO3 + NO" = "#ef6638" ) `);

foreach my $mechanism (sort keys %data) {
    $R->run(q` data = data.frame() `);
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        $R->run(q` pre = data.frame(Time) `); 
        foreach my $ref (@{$data{$mechanism}{$speciation}}) {
            foreach my $reaction (sort keys %$ref) {
                next if ($reaction eq "N2O5");
                $R->set('reaction', $reaction);
                $R->set('rates', [ map { $_ } $ref->{$reaction}->dog ]);
                $R->run(q` pre[reaction] = rates `);
            }
        }
        $R->set('speciation', $speciation);
        $R->run(q` pre$Speciation = rep(speciation, length(Time)) `,
                q` pre = gather(pre, Reaction, Rates, -Time, -Speciation) `,
                q` data = rbind(data, pre) `,
        );
    }
    #my $p = $R->run(q` print(data) `);
    #print $p, "\n";
    $R->set('filename', "Ox_budgets_all_speciations_${mechanism}.pdf");
    $R->set('title', "$mechanism : Ox Budgets");
    $R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);

    if ($mechanism eq "RADM2") {
        $R->run(q` data$Reaction = factor(data$Reaction, levels = c("Production Others", "HC3P + NO", "ACO3 + NO", "MO2 + NO", "HO2 + NO", "O3", "NO2 + OH", "O1D", "Consumption Others")) `);
    } elsif ($mechanism eq "MCM" or $mechanism eq "MOZART") {
        $R->run(q` data$Reaction = factor(data$Reaction, levels = c("Production Others", "CH3CO3 + NO", "CH3O2 + NO", "HO2 + NO", "O3", "NO2 + OH", "O1D", "Consumption Others")) `);
    }

    $R->run(q` plot = ggplot(data, aes(x = Time, y = Rates, fill = Reaction)) `,
            q` plot = plot + geom_bar(data = subset(data, Rates < 0), stat = "identity") `,
            q` plot = plot + geom_bar(data = subset(data, Rates > 0), stat = "identity") `,
            q` plot = plot + facet_wrap(~ Speciation, scales = "free") `,
            q` plot = plot + ggtitle(title) `,
            q` plot = plot + theme_tufte() `,
            q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
            q` plot = plot + scale_y_continuous(limits = c(-1e9, 2e9), breaks = seq(-1e9, 2e9, 5e8), expand = c(0, 0)) `,
            q` plot = plot + ylab("Reaction Rates (molecules cm-3 s-1)") `,
            q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
            q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
            q` plot = plot + theme(plot.title = element_text(face = "bold")) `,
            q` plot = plot + theme(axis.title.x = element_blank()) `,
            q` plot = plot + theme(axis.title.y = element_text(face = "bold")) `,
            q` plot = plot + theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 0.7, vjust = 0.8)) `,
            q` plot = plot + theme(axis.ticks.x = element_blank()) `,
            q` plot = plot + theme(legend.title = element_blank()) `,
            q` plot = plot + theme(legend.position = c(0.85, 0.14)) `,
            q` plot = plot + scale_fill_manual(values = my.colours, limits = levels(data$Reaction)) `,
    );

    $R->run(q` CairoPDF(file = filename, width = 11.3, height = 8) `,
            q` print(plot) `,
            q` dev.off() `,
    );
}

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism, $speciation) = @_;

    my ($producers, $producer_yields, %production_rates, $consumers, $consumer_yields, %consumption_rates);
    if (exists $families{"Ox_${mechanism}_$speciation"}) {
        $kpp->family({
                name    => "Ox_${mechanism}_$speciation",
                members => $families{"Ox_${mechanism}_$speciation"},
                weights => $weights{"Ox_${mechanism}_$speciation"},
        });
        $producers = $kpp->producing("Ox_${mechanism}_$speciation");
        $producer_yields = $kpp->effect_on("Ox_${mechanism}_$speciation", $producers);
        $consumers = $kpp->consuming("Ox_${mechanism}_$speciation");
        $consumer_yields = $kpp->effect_on("Ox_${mechanism}_$speciation", $consumers);
    } else {
        print "No Ox family for $mechanism and $speciation\n";
    }
    print "No producers found for $mechanism and $speciation\n" if (@$producers == 0);
    print "No consumers found for $mechanism and $speciation\n" if (@$consumers == 0);
    
    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        my $reaction_string = $kpp->reaction_string($reaction);
        my ($reactants, $products) = split / = /, $reaction_string;
        $production_rates{$reactants} += $rate(1:$ntime-2)->sumover;
    }

    for (0..$#$consumers) {
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
        next if ($rate->sum == 0);
        my $reaction_string = $kpp->reaction_string($reaction);
        my ($reactants, $products) = split / = /, $reaction_string;
        $consumption_rates{$reactants} += $rate(1:$ntime-2)->sumover;
    }
    remove_common_processes(\%production_rates, \%consumption_rates);

    my $others_max  = 5e8;
    foreach my $reaction (sort keys %production_rates) {
        if ($integrate->sum < $others_max) {
            $production_rates{"Production Others"} += $production_rates{$reaction};
            delete $production_rates{$reaction};
        }
    }

    foreach my $reaction (sort keys %consumption_rates) {
        my $reshape = $consumption_rates{$reaction}->reshape($n_per_day, $n_days);
        my $integrate = $reshape->sumover;
        if ($integrate->sum > -$others_max) {
            $consumption_rates{"Consumption Others"} += $integrate(0:13:2);
            delete $consumption_rates{$reaction};
        } else {
            $consumption_rates{$reaction} = $integrate(0:13:2);
        }
    }

    my $sort_function = sub { $_[0]->sum }; 
    my @final_sorted;
    my @sorted_cons = sort { &$sort_function($consumption_rates{$a}) <=> &$sort_function($consumption_rates{$b}) } keys %consumption_rates;
    foreach (@sorted_cons) {
        next if ($_ =~ /Others/);
        push @final_sorted, { $_ => $consumption_rates{$_} };
    }
    push @final_sorted, { "Consumption Others" => $consumption_rates{"Consumption Others"} } if (defined $consumption_rates{"Consumption Others"});

    my @sorted_prod = reverse sort { &$sort_function($production_rates{$a}) <=> &$sort_function($production_rates{$b}) } keys %production_rates;
    foreach (@sorted_prod) {
        next if ($_ =~ /Others/);
        push @final_sorted, { $_ => $production_rates{$_} };
    }
    push @final_sorted, { "Production Others" => $production_rates{"Production Others"} } if (defined $production_rates{"Production Others"});
    return \@final_sorted;
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
