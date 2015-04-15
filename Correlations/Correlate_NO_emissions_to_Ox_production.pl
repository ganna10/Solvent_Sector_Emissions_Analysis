#! /usr/bin/env perl
# Correlate cumulative NO emissions to Ox Production 
# Version 0: Jane Coates 14/4/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/work/users/jco/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my (%families, %weights, %NO_emissions, %Ox_production);

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $boxmodel = "$base/$mechanism/${speciation}_Solvents_Only/boxmodel";
        my $mecca = MECCA->new($boxmodel);
        my $eqn = "$base/$mechanism/${speciation}_Solvents_Only/gas.eqn";
        my $kpp = KPP->new($eqn);
        my $ro2_file = "$base/$mechanism/${speciation}_Solvents_Only/RO2_species.txt";
        my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
        $families{"Ox"} = [ qw(O3 O O1D NO2 NO3 N2O5 HO2NO2), @no2_reservoirs ];
        $weights{"Ox"} = { NO3 => 2, N2O5 => 3 };
        $NO_emissions{$mechanism}{$speciation} = get_NO_emissions($kpp, $mecca);
        $Ox_production{$mechanism}{$speciation} = get_Ox_production($kpp, $mecca);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
        q` library(ggthemes) `,
        q` library(grid) `,
);

$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %NO_emissions) {
    $R->set('mechanism', $mechanism);
    $R->run(q` pre = data.frame(Mechanism = mechanism) `);
    foreach my $speciation (sort keys %{$NO_emissions{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->set('no.emissions', $NO_emissions{$mechanism}{$speciation});
        $R->set('ox.production', $Ox_production{$mechanism}{$speciation});
        $R->run(q` pre$Speciation = speciation `,
                q` pre$NO.Emissions = no.emissions `,
                q` pre$Ox.Production = ox.production `,
        );
        $R->run(q` data = rbind(data, pre) `);
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` mechanism.colours = c("MCM" = "#6c254f", "MOZART" = "#ef6638", "RADM2" = "#0e5c28") `);
$R->run(q` speciation.colours = c("TNO" = "#000000", "IPCC" = "#2b9eb3", "EMEP" = "#ef6638", "DE94" = "#f9c500", "GR95" = "#6c254f", "GR05" = "#0352cb", "UK98" = "#0e5c28", "UK08" = "#b569b3") `);
$R->run(q` plot.lines = function () { list( theme_tufte() ,
                                            ylab("Ox Production (molecules cm-3)") ,
                                            xlab("NO Emissions (molecules cm-3)"),
                                            theme(legend.title = element_blank()) ,
                                            theme(axis.line = element_line(colour = "black")) ,
                                            theme(axis.title = element_text(face = "bold")) ,
                                            theme(strip.text = element_text(face = "bold")) ,
                                            theme(panel.margin = unit(5, "mm")) ) } `);

$R->run(q` plot = ggplot(data, aes(x = NO.Emissions, y = Ox.Production, colour = Speciation)) `,
        q` plot = plot + geom_point(size = 4) `,
        q` plot = plot + facet_wrap( ~ Mechanism, nrow = 1) `,
        q` plot = plot + scale_colour_manual(values = speciation.colours, limits = levels(data$Speciation)) `,
        q` plot = plot + plot.lines() `,
);


$R->run(q` CairoPDF(file = "NO_emissions_vs_Ox_production_facet_mechanism.pdf", width = 8.7, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->run(q` plot = ggplot(data, aes(x = NO.Emissions, y = Ox.Production, colour = Mechanism)) `,
        q` plot = plot + geom_point(size = 4) `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2) `,
        q` plot = plot + scale_colour_manual(values = mechanism.colours, limits = levels(data$Mechanism)) `,
        q` plot = plot + plot.lines() `,
);


$R->run(q` CairoPDF(file = "NO_emissions_vs_Ox_production_facet_speciation.pdf", width = 8.7, height = 6) `,
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

sub get_NO_emissions {
    my ($kpp, $mecca) = @_;
    my $emission_reaction = $kpp->producing_from("NO", "UNITY");
    my $reaction_number = $kpp->reaction_number($emission_reaction->[0]);
    my $emission_rate = $mecca->rate($reaction_number);
    return $emission_rate->sum;
}

sub get_Ox_production {
    my ($kpp, $mecca) = @_;
    $kpp->family({
            name    => "Ox",
            members => $families{"Ox"},
            weights => $weights{"Ox"},
    });
    my $producers = $kpp->producing("Ox");
    my $producer_yields = $kpp->effect_on("Ox", $producers);

    my $Ox_production;
    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        $Ox_production += $rate->sum;
    }
    return $Ox_production;
}
