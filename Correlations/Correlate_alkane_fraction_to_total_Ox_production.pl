#! /usr/bin/env perl
# Correlate fraction of alkane emissions to total Ox production
# Version 0: Jane Coates 26/4/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw( MCM );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my (%families, %weights, %Ox_production, %emissions);

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $dir = "$base/$mechanism/${speciation}_Solvents_Only";
        my $mecca = MECCA->new("$dir/boxmodel");
        my $kpp = KPP->new("$dir/gas.eqn");
        my $ro2_file = "$dir/RO2_species.txt";
        my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
        $families{"Ox"} = [ qw( O3 O O1D NO2 NO3 N2O5 HO2NO2 ), @no2_reservoirs ];
        $weights{"Ox"} = { NO3 => 2, N2O5 => 3 };
        $Ox_production{$mechanism}{$speciation} = get_Ox_production($kpp, $mecca);
        $emissions{$mechanism}{$speciation} = get_emission_fractions($kpp, $mecca);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(scales) `,
        q` library(grid) `,
        q` library(ggthemes) `,
);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %emissions) {
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$emissions{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->set('ox', $Ox_production{$mechanism}{$speciation});
        $R->run(q` pre = data.frame(Mechanism = mechanism, Speciation = speciation, Ox = ox) `);
        foreach my $type (sort keys %{$emissions{$mechanism}{$speciation}}) {
            $R->set('type', $type);
            $R->set('emission.fraction', $emissions{$mechanism}{$speciation}{$type});
            $R->run(q` pre[type] = emission.fraction `);
        }
        $R->run(q` pre = gather(pre, Type, Emission.Fraction, -Mechanism, -Speciation, -Ox) `,
                q` data = rbind(data, pre) `,
        );
    }
}
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` mechanism.colours = c("MCM" = "#6c254f", "MOZART" = "#ef6638", "RADM2" = "#0e5c28") `);
$R->run(q` speciation.colours = c("TNO" = "#000000", "IPCC" = "#2b9eb3", "EMEP" = "#ef6638", "DE94" = "#f9c500", "GR95" = "#6c254f", "GR05" = "#0352cb", "UK98" = "#0e5c28", "UK08" = "#b569b3") `);
$R->run(q`  lm_eqn = function(m) {  l <- list(a = format(coef(m)[1], digits = 2), 
                                    b = format(abs(coef(m)[2]), digits = 2), 
                                    r2 = format(summary(m)$r.squared, digits = 3));
                                    if (coef(m)[2] >= 0)  {
                                        eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
                                    } else {
                                        eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
                                    }
                                    as.character(as.expression(eq));    } `);

$R->run(q` plot.lines = function () { list( theme_tufte() ,
                                            ylab("Ox Production (molecules cm-3)") ,
                                            theme(legend.title = element_blank()) ,
                                            theme(axis.line = element_line(colour = "black")) ,
                                            theme(axis.title = element_text(face = "bold")) ,
                                            theme(strip.text = element_text(face = "bold")) ,
                                            theme(panel.margin = unit(5, "mm")) ) } `);

$R->run(q` alkanes = filter(data, Type == "Alkanes") `);
#my $p = $R->run(q` print(alkanes) `);
#print $p, "\n"; 

$R->run(q` plot = ggplot(alkanes, aes(x = Emission.Fraction, y = Ox, colour = Speciation)) `,
        q` plot = plot + geom_point(size = 4) `,
        q` plot = plot + scale_x_continuous(label = percent) `,
        q` plot = plot + xlab("Alkane Emissions") `,
        q` plot = plot + scale_colour_manual(values = speciation.colours, limits = levels(data$Speciation)) `,
        q` plot = plot + geom_text(aes(x = 0.3, y = 8.7e9, label = lm_eqn(lm(Ox ~ Emission.Fraction, alkanes))), colour = "black", parse = TRUE) + geom_smooth(method = "lm", se = FALSE, colour = "black") `,
        q` plot = plot + plot.lines() `,
);

$R->run(q` CairoPDF(file = "Ox_production_vs_alkane_emission_fraction.pdf", width = 8.7, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->run(q` oxygenated = filter(data, Type == "Oxygenated") `);
$R->run(q` plot = ggplot(oxygenated, aes(x = Emission.Fraction, y = Ox, colour = Speciation)) `,
        q` plot = plot + geom_point(size = 4) `,
        q` plot = plot + scale_x_continuous(label = percent) `,
        q` plot = plot + xlab("Oxygenated Emissions") `,
        q` plot = plot + scale_colour_manual(values = speciation.colours, limits = levels(data$Speciation)) `,
        q` plot = plot + geom_text(aes(x = 0.4, y = 8.7e9, label = lm_eqn(lm(Ox ~ Emission.Fraction, oxygenated))), colour = "black", parse = TRUE) + geom_smooth(method = "lm", se = FALSE, colour = "black") `,
        q` plot = plot + plot.lines() `,
);

$R->run(q` CairoPDF(file = "Ox_production_vs_oxygenated_emission_fraction.pdf", width = 8.7, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->run(q` aromatic = filter(data, Type == "Aromatics") `);
$R->run(q` plot = ggplot(aromatic, aes(x = Emission.Fraction, y = Ox, colour = Speciation)) `,
        q` plot = plot + geom_point(size = 4) `,
        q` plot = plot + scale_x_continuous(label = percent) `,
        q` plot = plot + xlab("Aromatic Emissions") `,
        q` plot = plot + scale_colour_manual(values = speciation.colours, limits = levels(data$Speciation)) `,
        q` plot = plot + geom_text(aes(x = 0.15, y = 8.7e9, label = lm_eqn(lm(Ox ~ Emission.Fraction, aromatic))), colour = "black", parse = TRUE) + geom_smooth(method = "lm", se = FALSE, colour = "black") `,
        q` plot = plot + plot.lines() `,
);

$R->run(q` CairoPDF(file = "Ox_production_vs_aromatic_emission_fraction.pdf", width = 8.7, height = 6) `,
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

sub get_Ox_production {
    my ($kpp, $mecca) = @_;
    $kpp->family({
            name    => "Ox",
            members => $families{"Ox"},
            weights => $weights{"Ox"},
    });
    my $producers = $kpp->producing("Ox");
    my $producer_yields = $kpp->effect_on("Ox", $producers);
    print "No producers\n" if (@$producers == 0);

    my $total_Ox = 0;
    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        $total_Ox += $rate->sum;
    }
    return $total_Ox;
}

sub get_emission_fractions {
    my ($kpp, $mecca) = @_;
    my %emissions;

    my $emitting_reactions = $kpp->consuming("UNITY");
    for (0..$#$emitting_reactions) {
        my $reaction = $emitting_reactions->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number);
        next if ($rate->sum == 0);
        my $VOC = $kpp->products($reaction)->[0];
        next if ($VOC eq "NO"); 
        if ($VOC =~ /^NC|^IC|^NEOP|^M2|^M3/ or $VOC eq "C3H8" or $VOC eq "CHEX" or $VOC eq "C2H6") {
            #print "$VOC : alkane\n";
            $emissions{"Alkanes"} += $rate->sum;
        } elsif ($VOC eq "C2H4" or $VOC eq "C3H6" or $VOC eq "C2H2" or $VOC eq "LIMONENE") {
            #print "$VOC : alkene or alkyne\n";
            $emissions{"Alkenes, Alkynes"} += $rate->sum;
        } elsif ($VOC =~ /^TOL|XYL|^TM|BENZ|^STY|ETHTOL$|^DIM/) {
            #print "$VOC : aromatic\n";
            $emissions{"Aromatics"} += $rate->sum;
        } elsif ($VOC eq "CH2CL2" or $VOC eq "CH3CCL3" or $VOC eq "TCE" or $VOC eq "TRICLETH") {
            #print "$VOC : chlorinated\n";
            $emissions{"Chlorinated"} += $rate->sum;
        } else {
            #print "$VOC : oxygenated\n";
            $emissions{"Oxygenated"} += $rate->sum;
        }
    }
    my $total_emissions = 0;
    $total_emissions += $emissions{$_} foreach (keys %emissions);
    $emissions{$_} /= $total_emissions foreach (keys %emissions);
    return \%emissions;
}
