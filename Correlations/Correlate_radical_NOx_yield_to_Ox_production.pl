#! /usr/bin/env perl
# Correlate net radical to NOx yield as used in NO emissions to Ox production
# Version 0: Jane Coates 15/4/15

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/work/users/jco/Solvent_Emissions";
#my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my (%families, %weights, %radical_NOx, %Ox_production);

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $dir = "$base/$mechanism/${speciation}_Solvents_Only";
        my $boxmodel = "$dir/boxmodel";
        my $mecca = MECCA->new($boxmodel);
        my $eqn = "$dir/gas.eqn";
        my $kpp = KPP->new($eqn);
        my $ro2_file = "$dir/RO2_species.txt";
        my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
        my $radical_file = "$dir/radicals.txt";
        my @radicals = get_radicals($radical_file);
        %families = (
            "NOx"       => [ qw( NO NO2 NO3 N2O5 ) ],
            "Ox"        => [ qw( O3 O O1D NO2 NO3 N2O5 HO2NO2 ), @no2_reservoirs ],
            "radicals"  => [ @radicals ],
        );
        %weights = (
            "NOx"   => { N2O5 => 2 },
            "Ox"    => { NO3 => 2, N2O5 => 3 },
        );
        $radical_NOx{$mechanism}{$speciation} = get_radical_NOx_yield($kpp, $mecca);
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
foreach my $mechanism (sort keys %radical_NOx) {
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$radical_NOx{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->set('radical.nox', $radical_NOx{$mechanism}{$speciation});
        $R->set('ox.production', $Ox_production{$mechanism}{$speciation});
        $R->run(q` pre = data.frame(Mechanism = mechanism) `,
                q` pre$Speciation = speciation `,
                q` pre$Radical.NOx = radical.nox `,
                q` pre$Ox.Production = ox.production `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` mechanism.colours = c("MCM" = "#6c254f", "MOZART" = "#ef6638", "RADM2" = "#0e5c28") `);
$R->run(q` speciation.colours = c("TNO" = "#000000", "IPCC" = "#2b9eb3", "EMEP" = "#ef6638", "DE94" = "#f9c500", "GR95" = "#6c254f", "GR05" = "#0352cb", "UK98" = "#0e5c28", "UK08" = "#b569b3") `);
$R->run(q` plot.lines = function () { list( theme_tufte() ,
                                            ylab("Ox Production (molecules cm-3)") ,
                                            xlab("Radical to NOx Yield (molecules cm-3)"),
                                            theme(legend.title = element_blank()) ,
                                            theme(axis.line = element_line(colour = "black")) ,
                                            theme(axis.title = element_text(face = "bold")) ,
                                            theme(strip.text = element_text(face = "bold")) ,
                                            theme(panel.margin = unit(5, "mm")) ) } `);

$R->run(q` plot = ggplot(data, aes(x = Radical.NOx, y = Ox.Production, colour = Speciation)) `,
        q` plot = plot + geom_point(size = 4) `,
        q` plot = plot + facet_wrap( ~ Mechanism, nrow = 1) `,
        q` plot = plot + plot.lines() `,
        q` plot = plot + scale_colour_manual(values = speciation.colours, limits = levels(data$Speciation)) `,
);

$R->run(q` CairoPDF(file = "Radical_NOx_yield_vs_Ox_Production_facet_mechanism.pdf", width = 8.7, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->run(q` plot = ggplot(data, aes(x = Radical.NOx, y = Ox.Production, colour = Mechanism)) `,
        q` plot = plot + geom_point(size = 4) `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2) `,
        q` plot = plot + plot.lines() `,
        q` plot = plot + scale_colour_manual(values = mechanism.colours, limits = levels(data$Mechanism)) `,
);

$R->run(q` CairoPDF(file = "Radical_NOx_yield_vs_Ox_Production_facet_speciation.pdf", width = 8.7, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_radicals {
    my ($file) = @_;
    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file : $!";
    my $line = <$in>;
    close $in;
    my @radicals = split /\s+/, $line;
    return @radicals;
}

sub get_radical_NOx_yield {
    my ($kpp, $mecca) = @_;
    
    $kpp->family({ #NOx family
                    name    => 'NOx',
                    members =>$families{'NOx'},
                    weights => $weights{'NOx'},
    });

    my ($radical_producers, $radical_producer_yields, $NOx_yields, $production);
    $kpp->family({                                                                                                                           
                name    => "radicals",
                members => $families{"radicals"},
                weights => $weights{"radicals"},
    }); 
    $radical_producers = $kpp->producing("radicals");
    $radical_producer_yields = $kpp->effect_on("radicals", $radical_producers);
    $NOx_yields = $kpp->effect_on('NOx', $radical_producers);
    
    for (0..$#$radical_producers) { #get rates for all radical producing reactions
        next if ($NOx_yields->[$_] < 0);
        my $reaction = $radical_producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $net_radical_yield = $radical_producer_yields->[$_] - $NOx_yields->[$_];
        next if ($net_radical_yield == 0);
        my $rate = $net_radical_yield * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        $production += $rate->sum;
    }
    return $production;
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
