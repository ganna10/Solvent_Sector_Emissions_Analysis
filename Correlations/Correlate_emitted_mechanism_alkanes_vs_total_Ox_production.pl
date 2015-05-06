#! /usr/bin/env perl
# Correlate total Ox production to total amount of emitted alkanes in each speciation
# Version 0: Jane Coates 6/5/2015

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use KPP;
use MECCA;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my (%families, %weights, %emissions , %Ox);

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $dir = "$base/$mechanism/${speciation}_Solvents_Only";
        my $mecca = MECCA->new("$dir/boxmodel");
        my $kpp = KPP->new("$dir/gas.eqn");
        my $ro2_file = "$dir/RO2_species.txt";
        my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
        $families{"Ox"} = [ qw( O3 O O1D NO2 NO3 N2O5 HO2NO2 ), @no2_reservoirs ];
        $weights{"Ox"} = { NO3 => 2, N2O5 => 3 };
        $Ox{$mechanism}{$speciation} = get_Ox($mecca, $kpp);
        $emissions{$mechanism}{$speciation} = get_alkane_emissions($mecca, $kpp);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
        q` library(grid) `,
        q` library(scales) `,
);
$R->run(q` data = data.frame() `);

foreach my $mechanism (sort keys %Ox) {
    $R->set('mechanism', $mechanism);
    $R->run(q` pre = data.frame(Mechanism = mechanism) `);
    foreach my $speciation (sort keys %{$Ox{$mechanism}}) {
        $R->set('speciation' ,$speciation);
        $R->set('ox', $Ox{$mechanism}{$speciation});
        $R->set('alkane.emissions', $emissions{$mechanism}{$speciation});
        $R->run(q` pre$Speciation = speciation `,
                q` pre$Ox = ox `,
                q` pre$Alkane.Emissions = alkane.emissions `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` data = mutate(data, mechanism = factor(Mechanism, labels = c("MCM v3.2", "MOZART-4", "RADM2"))) `);
$R->run(q` my.colours = c("MCM v3.2" = "#6c254f", "MOZART-4" = "#ef6638", "RADM2" = "#0e5c28") `);

$R->run(q` plot.lines = function () { list( theme_tufte() ,
                                            ylab("Total Ox Production (molecules cm-3)") ,
                                            theme(legend.title = element_blank()) ,
                                            theme(axis.line = element_line(colour = "black")) ,
                                            theme(axis.title = element_text(face = "bold")) ,
                                            theme(strip.text = element_text(face = "bold")) ,
                                            theme(panel.margin = unit(5, "mm")) ) } `);

$R->run(q` eqn = data %>% group_by(Speciation) %>% do(mod2 = cor(.$Ox, .$Alkane.Emissions, method = "pearson")) %>% mutate(Label = paste('"r = ', sprintf('%1.3f', mod2[1]), '"')) %>% select(-mod2) `);
#my $p = $R->run(q` print.data.frame(eqn) `);
#print $p, "\n"; 

$R->run(q` plot = ggplot(data, aes(x = Alkane.Emissions, y = Ox, colour = mechanism)) `,
        q` plot = plot + geom_point(size = 4) `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2) `,
        q` plot = plot + scale_x_continuous(label = percent) `,
        q` plot = plot + xlab("Percent Emissions") `,
        q` plot = plot + scale_colour_manual(values = my.colours, limits = levels(data$mechanism)) `,
        q` plot = plot + plot.lines() `,
        q` plot = plot + theme(strip.text.x = element_text(face = "bold")) `,
        q` plot = plot + theme(legend.position = "top") `,
        q` plot = plot + geom_text(data = eqn, aes(x = 0.55, y = 8.7e9, label = Label), colour = "black", inherit.aes = FALSE, parse = TRUE) + geom_smooth(method = "lm", se = FALSE, colour = "black") `,
);

$R->run(q` CairoPDF(file = "Ox_production_vs_emitted_alkanes.pdf", width = 10, height = 7) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

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

sub get_Ox {
    my ($mecca, $kpp) = @_;
    $kpp->family({
            name    => "Ox",
            members => $families{"Ox"},
            weights => $weights{"Ox"},
    });
    my $producers = $kpp->producing("Ox");
    my $producer_yields = $kpp->effect_on("Ox", $producers);
    print "No producers\n" if (@$producers == 0);

    my $total_Ox;
    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        $total_Ox += $rate->sum;
    }
    return $total_Ox;
}

sub get_alkane_emissions {
    my ($mecca, $kpp) = @_;
    my @alkanes = qw( C2H6 C3H8 NC4H10 IC4H10 HC5H12 IC5H12 NEOP NC6H14 M2PE M3PE NC7H16 M2HEX M3HEX NC8H18 NC920 NC10H22 NC11H24 NC12H26 CHEX BIGALK ETH HC3 HC5 HC8 );
    my ($total_alkane_emissions, $total_emissions);
    my $emission_reactions = $kpp->consuming("UNITY");
    for (0..$#$emission_reactions) {
        my $reaction = $emission_reactions->[$_];
        my ($unity, $product) = split  / = /, $kpp->reaction_string($reaction);
        next if ($product =~ /^NO/);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate =  $mecca->rate($reaction_number);
        next if ($rate->sum == 0);
        $total_emissions += $rate->sum;
        $total_alkane_emissions += $rate->sum if ($product ~~ @alkanes);
    }

    return $total_alkane_emissions / $total_emissions;
}
