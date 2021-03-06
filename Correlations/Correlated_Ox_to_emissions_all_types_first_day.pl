#! /usr/bin/env perl
# Correlate fraction of all types of VOC emissions to first day Ox production in all mechanisms
# Version 0: Jane Coates 22/7/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/work/users/jco/Solvent_Emissions";
#my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw( MOZART RADM2 );
#my @mechanisms = qw( MCM MOZART RADM2 );
#my @speciations = qw( EMEP UK08 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my (%families, %weights, %Ox_production, %emissions);

my $mecca = MECCA->new("$base/MOZART/DE94_Solvents_Only/boxmodel");
my $ntime = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $n_per_day = 86400 / $dt;
my $n_days = int $ntime / $n_per_day;

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
        $emissions{$mechanism}{$speciation} = get_emission_fractions($kpp, $mecca, $mechanism, $speciation);
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
        $R->run(q` pre = data.frame(Mechanisms = mechanism, Speciation = speciation, Ox = ox) `);
        foreach my $type (sort keys %{$emissions{$mechanism}{$speciation}}) {
            $R->set('type', $type);
            $R->set('emission.fraction', $emissions{$mechanism}{$speciation}{$type});
            $R->run(q` pre[type] = emission.fraction `);
        }
        $R->run(q` pre = gather(pre, Type, Emission.Fraction, -Mechanisms, -Speciation, -Ox) `,
                q` data = rbind(data, pre) `,
        );
    }
}
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` data = data %>% mutate(Mechanism = factor(Mechanisms, labels = c("MOZART-4", "RADM2")))  %>% select(-Mechanisms) `);
#$R->run(q` data = mutate(data, mechanism = factor(Mechanism, labels = c("MCM v3.2", "MOZART-4", "RADM2"))) `);

#theme(legend.title = element_blank()) ,
$R->run(q` plot.lines = function () { list( theme_tufte() ,
                                            ylab("Total Ox Production (molecules cm-3)") ,
                                            theme(axis.line = element_line(colour = "black")) ,
                                            theme(axis.title = element_text(face = "bold")) ,
                                            theme(strip.text = element_text(face = "bold")) ,
                                            theme(panel.margin = unit(5, "mm")) ) } `);

$R->run(q` eqn = data %>% group_by(Type) %>% do(mod2 = cor(.$Ox, .$Emission.Fraction, method = "pearson")) %>% mutate(Label = paste('"r = ', sprintf('%1.3f', mod2[1]), '"')) %>% select(-mod2) `);
#my $p = $R->run(q` print.data.frame(eqn) `);
#print $p, "\n"; 
$R->run(q` speciation.colours = c("TNO" = "#000000", "IPCC" = "#2b9eb3", "EMEP" = "#ef6638", "DE94" = "#f9c500", "GR95" = "#6c254f", "GR05" = "#0352cb", "UK98" = "#0e5c28", "UK08" = "#b569b3") `);
#$R->run(q` speciation.colours = c("#000000", "#000000", "#2b9eb3", "#2b9eb3", "#ef6638", "#ef6638", "#f9c500", "#f9c500", "#6c254f", "#6c254f", "#0352cb", "#0352cb", "#0e5c28", "#0e5c28", "#b569b3", "#b569b3") `);
$R->run(q` shapes.colours = c("TNO, MOZART-4", "TNO, RADM2", "IPCC, MOZART-4", "IPCC, RADM2", "EMEP, MOZART-4", "EMEP, RADM2", "DE94, MOZART-4", "DE94, RADM2", "GR95, MOZART-4", "GR95, RADM2", "GR05, MOZART-4", "GR05, RADM2", "UK98, MOZART-4", "UK98, RADM2", "UK08, MOZART-4", "UK08, RADM2") `);

$R->run(q` plot = ggplot(data, aes(x = Emission.Fraction, y = Ox)) `,
        q` plot = plot + geom_point(size = 4, aes(colour = Speciation, shape = Mechanism)) `,
        q` plot = plot + facet_wrap( ~ Type, scales = "free", ncol = 2) `,
        q` plot = plot + scale_x_continuous(label = percent, limits = c(0, 0.7)) `,
        q` plot = plot + scale_y_continuous(limits = c(1e9, 2e9)) `,
        q` plot = plot + xlab("Percent Emissions") `,
        q` plot = plot + scale_colour_manual(values = speciation.colours, limits = levels(data$Speciation)) `,
        #q` plot = plot + scale_shape_manual(values = c(16, 17, 16, 17, 16, 17, 16, 17, 16, 17, 16, 17, 16, 17, 16, 17), labels = shapes.colours) `, #16 = circle, 17 = triangle
        q` plot = plot + plot.lines() `,
        q` plot = plot + theme(legend.position = "top") `,
        q` plot = plot + geom_text(data = eqn, aes(x = 0.45, y = 2.0e9, label = Label), colour = "black", inherit.aes = FALSE, parse = TRUE) + geom_smooth(method = "lm", se = FALSE, colour = "black") `,
);

$R->run(q` CairoPDF(file = "First_day_Ox_production_vs_all_VOC_emission_fractions_all_data.pdf", width = 7, height = 10) `,
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
    print "Ox\n";

    my $total_Ox = 0;
    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        $rate = $rate->reshape($n_per_day, $n_days);
        $rate = $rate->sumover;
        $total_Ox += $rate->at(0);
    }
    return $total_Ox;
}

sub get_emission_fractions {
    my ($kpp, $mecca, $mechanism, $speciation) = @_;
    my %emissions;
    my @alkanes = qw( C2H6 C3H8 NC4H10 IC4H10 NC5H12 IC5H12 NEOP NC6H14 M2PE M3PE NC7H16 M2HEX M3HEX NC8H18 NC9H20 NC10H22 NC11H24 NC12H26 BIGALK ETH HC3 HC5 HC8 CHEX );
    my @alkenes = qw( C2H4 C3H6 BUT1ENE LIMONENE APINENE BPINENE BIGENE OL2 OLT OLI C10H16 C2H2 ISOP ISO C5H8 );
    my @aromatics = qw( BENZENE TOLUENE MXYL OXYL PXYL EBENZ PBENZ IPBENZ TM123B TM124B TM135B STYRENE METHTOL OETHTOL PETHTOL DIME35EB TOL XYL );
    my @chlorinated = qw( CH2CL2 CH3CCL3 TCE TRICLETH );
    my @oxygenated = qw( ETHACET NBUTACET IPROACET NPROACET CH3COCH3 MEK MIBK CYHEXONE CH3OH C2H5OH IPROPOL NBUTOL NPROPOL BUT2OL IBUTOL MIBKAOH C6H5CH2OH ETHGLY PROPGLY BUOX2ETOH PR2OHMOX EOX2EOL CH3OCH3 MO2EOL HCOOH CH3CO2H PROPACID HCHO CH3CHO C2H5CHO C3H7CHO IPRCHO C4H9CHO ACR MACR C4ALDB CH2O HCHO ALD KET CH3COOH ORA1 ORA2 );

    print "Emissions: $mechanism  : $speciation\n";
    my $emitting_reactions = $kpp->consuming("UNITY");
    for (0..$#$emitting_reactions) {
        my $reaction = $emitting_reactions->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number);
        next if ($rate->sum == 0);
        my $VOC = $kpp->products($reaction)->[0];
        next if ($VOC eq "NO"); 
        if ($VOC ~~ @alkanes) {
            #print "$VOC : alkane\n";
            $emissions{"Alkanes"} += $rate->sum;
        } elsif ($VOC ~~ @alkenes) {
            #print "$VOC : alkene or alkyne\n";
            $emissions{"Alkenes, Alkynes"} += $rate->sum;
        } elsif ($VOC ~~ @aromatics) {
            #print "$VOC : aromatic\n";
            $emissions{"Aromatics"} += $rate->sum;
        } elsif ($VOC ~~ @chlorinated) {
            #print "$VOC : chlorinated\n";
            $emissions{"Chlorinated"} += $rate->sum;
        } elsif ($VOC ~~ @oxygenated) {
            #print "$VOC : oxygenated\n";
            $emissions{"Oxygenated"} += $rate->sum;
        } else {
            print "No match for $VOC\n";
        }
    }

    if ($mechanism eq "RADM2") {
        my $contribution_file = "Fractional_Contributions_of_individual_RADM2_species_to_functional_groups.csv";
        open my $in, '<:encoding(utf-8)', $contribution_file or die $!;
        my @lines = <$in>;
        close $in;
        my %RADM2;
        foreach my $line (@lines) {
            next if ($line =~ /^Speci/);
            chomp $line;
            my ($speciation, $group, $contribution) = split /,/, $line;
            $RADM2{$speciation}{$group} = $contribution;
        }

        $emissions{"Chlorinated"} += $RADM2{$speciation}{"Chlorinated"} * $emissions{"Alkanes"};
        $emissions{"Oxygenated"} += $RADM2{$speciation}{"Oxygenated"} * $emissions{"Alkanes"};
        $emissions{"Alkenes, Alkynes"} += $RADM2{$speciation}{"Alkenes;Alkynes"} * $emissions{"Alkanes"};
        $emissions{"Alkanes"} += $RADM2{$speciation}{"Alkanes"} * $emissions{"Alkanes"};
    } elsif ($mechanism eq "MOZART") {
        my $contribution_file = "Fractional_Contributions_of_individual_MOZART_species_to_functional_groups.csv";
        open my $in, '<:encoding(utf-8)', $contribution_file or die $!;
        my @lines = <$in>;
        close $in;
        my %MOZART;
        foreach my $line (@lines) {
            next if ($line =~ /^Speci/);
            chomp $line;
            my ($speciation, $group, $contribution) = split /,/, $line;
            $MOZART{$speciation}{$group} = $contribution;
        }

        $emissions{"Chlorinated"} += $MOZART{$speciation}{"Chlorinated"} * $emissions{"Alkanes"};
        $emissions{"Oxygenated"} += $MOZART{$speciation}{"Oxygenated"} * $emissions{"Alkanes"};
        $emissions{"Alkanes"} += $MOZART{$speciation}{"Alkanes"} * $emissions{"Alkanes"};
    }
    my $total_emissions = 0;
    $total_emissions += $emissions{$_} foreach (keys %emissions);
    $emissions{$_} /= $total_emissions foreach (keys %emissions);
    return \%emissions;
}
