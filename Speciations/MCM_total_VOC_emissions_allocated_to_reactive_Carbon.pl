#! /usr/bin/env perl
# Get total emission rate of all VOC in speciations and allocate by reactive carbon, correlated to first day NO emissions and Ox production
# Version 0: Jane Coates 21/4/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions/MCM";
#my @speciations = qw( EMEP );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my (%data, %families, %weights, %Ox);

my $carbons_file = "$base/TNO_Solvents_Only/carbons.txt";
my $carbons = mcm_n_carbon($carbons_file);

foreach my $speciation (@speciations) {
    my $dir = "$base/${speciation}_Solvents_Only";
    my $boxmodel = "$dir/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn = "$dir/gas.eqn";
    my $kpp = KPP->new($eqn);
    my $ro2_file = "$dir/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
    $families{"Ox"} = [ qw( O3 O O1D NO2 NO3 N2O5 HO2NO2 ), @no2_reservoirs ];
    $weights{"Ox"} = { NO3 => 2, N2O5 => 3 };
    $data{$speciation} = get_voc_emissions($mecca, $kpp, $speciation);
    $Ox{$speciation} = get_Ox_production($kpp, $mecca);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
);

$R->run(q` data = data.frame() `);
foreach my $speciation (keys %data) {
    $R->set('speciation', $speciation);
    $R->run(q` pre = data.frame(Speciation = speciation) `);
    foreach my $VOC (sort keys %{$data{$speciation}}) {
        next if ($VOC eq "NO");
        print "No carbon for $VOC\n" unless (defined $carbons->{$VOC});
        $R->set('voc', $VOC);
        $R->set('emissions', $data{$speciation}{$VOC});
        $R->set('carbon', $carbons->{$VOC});
        $R->run(q` pre$VOC = voc `,
                q` pre$Emissions = emissions `,
                q` pre$Carbon = carbon `,
                q` data = rbind(data, pre) `,
        );
    }
}
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` data = data %>% group_by(Carbon, Speciation) %>% summarize(Total.Emissions = sum(Emissions)) `);
$R->run(q` data$Carbon = factor(data$Carbon, levels = rev(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))) `);
$R->run(q` my.colours = c("12" = "#ba8b01", "11" = "#4c9383", "10" = "#ae4901", "9" = "#0352cb", "8" = "#6db875", "7" = "#0c3f78", "6" = "#b569b3", "5" = "#2b9eb3", "4" = "#ef6638", "3" = "#0e5628", "2" = "#f9c500", "1" = "#6c254f") `);
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` plot = ggplot(data, aes(x = Speciation, y = Total.Emissions, order = -as.numeric(Carbon), fill = Carbon)) `,
        q` plot = plot + geom_bar(stat = "identity", position = "stack") `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + ylab("Total VOC emissions (molecules cm-3)") `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.ticks.x = element_blank()) `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.text.x = element_text(face = "bold")) `,
        q` plot = plot + scale_fill_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "MCM_VOC_emissions_vs_reactive_C.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

###Correlate total VOC emissions with total NO emissions
foreach my $speciation (keys %data) {
    my $total_C = 0;
    foreach my $VOC (keys %{$data{$speciation}}) {
        next if ($VOC eq "NO");
        $total_C += $data{$speciation}{$VOC} * $carbons->{$VOC};
        delete $data{$speciation}{$VOC};
    }
    $data{$speciation}{"Total.C"} = $total_C;
}

$R->run(q` no.data = data.frame() `);
foreach my $speciation (keys %data) {
    $R->set('speciation', $speciation);
    $R->run(q` pre = data.frame(Speciation = speciation) `);
    foreach my $item (sort keys %{$data{$speciation}}) {
        $R->set('emissions', $data{$speciation}{$item});
        if ($item eq "NO") {
            $R->run(q` pre$NO.emissions = emissions `);
        } elsif ($item =~ "Total") {
            $R->run(q` pre$Reactive.C = emissions `);
        }
        $R->set('ox', $Ox{$speciation});
        $R->run(q` pre$Ox = ox `);
    }
    $R->run(q` no.data = rbind(no.data, pre) `);
}
$R->run(q` no.data$Speciation = factor(no.data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);

#my $p = $R->run(q` print(no.data) `);
#print $p, "\n";
$R->run(q` speciation.colours = c("TNO" = "#000000", "IPCC" = "#2b9eb3", "EMEP" = "#ef6638", "DE94" = "#f9c500", "GR95" = "#6c254f", "GR05" = "#0352cb", "UK98" = "#0e5c28", "UK08" = "#b569b3") `);
$R->run(q` lm_eqn = function(m) { l <- list(a = format(coef(m)[1], digits = 2),
                                      b = format(abs(coef(m)[2]), digits = 2),
                                      r2 = format(summary(m)$r.squared, digits = 3));
                                  if (coef(m)[2] >= 0)  {
                                    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
                                  } else {
                                    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
                                  }
                                  as.character(as.expression(eq));                 
                                } `);

# Total Reactive C vs NO emissions
$R->run(q` p = ggplot(no.data, aes(x = Reactive.C, y = NO.emissions, colour = Speciation)) `,
        q` p = p + geom_point(size = 4) `,
        q` p = p + theme_tufte() `,
        q` p = p + theme(axis.line = element_line(colour = "black")) `,
        q` p = p + geom_text(aes(x = 8e8, y = 3.5e8, label = lm_eqn(lm(Reactive.C ~ NO.emissions, no.data))), colour = "black", parse = TRUE) + geom_smooth(method = "lm", se = FALSE, colour = "black") `,
        q` p = p + scale_colour_manual(values = speciation.colours) `,
);

$R->run(q` CairoPDF(file = "NO_vs_reactive_C.pdf") `,
        q` print(p) `,
        q` dev.off() `,
);

# Total reactive C vs Ox production
$R->run(q` p1 = ggplot(no.data, aes(x = Reactive.C, y = Ox, colour = Speciation)) `,
        q` p1 = p1 + geom_point(size = 4) `,
        q` p1 = p1 + theme_tufte() `,
        q` p1 = p1 + theme(axis.line = element_line(colour = "black")) `,
        q` p1 = p1 + geom_text(aes(x = 8e8, y = 1.5e9, label = lm_eqn(lm(Reactive.C ~ Ox, no.data))), colour = "black", parse = TRUE) + geom_smooth(method = "lm", se = FALSE, colour = "black") `,
        q` p1 = p1 + scale_colour_manual(values = speciation.colours) `,
);

$R->run(q` CairoPDF(file = "Ox_vs_reactive_C.pdf") `,
        q` print(p1) `,
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

sub mcm_n_carbon {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die $!;
    my @lines = (<$in>);
    close $in;
    my %carbons;
    foreach my $line (@lines) {
        chomp $line;
        my ($species, $smile) = split ' ', $line;
        my $C_number = 0;
        if ($smile =~ /\./) {
            $C_number = 8;
        } else {
            $C_number++ while ($smile =~ m/C/gi);
        }
        $carbons{$species} = $C_number;
    }
    return \%carbons;
}

sub get_voc_emissions {
    my ($mecca, $kpp, $speciation) = @_;

    my @no_match = qw( NO2 CO CH4 );
    my $emissions = $kpp->consuming("UNITY");
    my %emissions;
    for (0..$#$emissions) {
        my $reaction = $emissions->[$_];
        my $products = $kpp->products($reaction);
        my $VOC = $products->[0];
        next if ($VOC ~~ @no_match);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number);
        next if ($rate->sum == 0);
        if ($VOC eq "NO") {
            $rate = $rate(1:73);
            $emissions{$VOC} += $rate->sum;
        } else {
            $emissions{$VOC} += $rate->sum;
        }
    }
    return \%emissions;
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
        $rate = $rate(1:73);
        $Ox_production += $rate->sum;
    }
    return $Ox_production;
}
