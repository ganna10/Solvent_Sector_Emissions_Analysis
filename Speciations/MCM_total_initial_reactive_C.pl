#! /usr/bin/env perl
# Compare starting initial reactive C levels of VOC emissions
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
my (%data);

my $carbons_file = "$base/TNO_Solvents_Only/carbons.txt";
my $carbons = mcm_n_carbon($carbons_file);

foreach my $speciation (@speciations) {
    my $dir = "$base/${speciation}_Solvents_Only";
    my $boxmodel = "$dir/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn = "$dir/gas.eqn";
    my $kpp = KPP->new($eqn);
    $data{$speciation} = get_voc_emissions($mecca, $kpp, $speciation);
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
        print "No carbon for $VOC\n" unless (defined $carbons->{$VOC});
        $R->set('voc', $VOC);
        $R->set('emissions', $data{$speciation}{$VOC});
        $R->set('carbon', $carbons->{$VOC});
        $R->run(q` pre$VOC = voc `,
                q` pre$Carbon = emissions * carbon `,
                q` data = rbind(data, pre) `,
        );
    }
}
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` data = data %>% group_by(Speciation) %>% summarize(Reactive.Carbon = sum(Carbon)) `);
#$R->run(q` data$Carbon = factor(data$Reactive.Carbon, levels = rev(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))) `);
#$R->run(q` my.colours = c("12" = "#ba8b01", "11" = "#4c9383", "10" = "#ae4901", "9" = "#0352cb", "8" = "#6db875", "7" = "#0c3f78", "6" = "#b569b3", "5" = "#2b9eb3", "4" = "#ef6638", "3" = "#0e5628", "2" = "#f9c500", "1" = "#6c254f") `);
my $p = $R->run(q` print(data) `);
print $p, "\n";

$R->run(q` plot = ggplot(data, aes(x = Speciation, y = Reactive.Carbon)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        #q` plot = plot + ylab("Total VOC emissions (molecules cm-3)") `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.ticks.x = element_blank()) `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.text.x = element_text(face = "bold")) `,
);

$R->run(q` CairoPDF(file = "MCM_initial_reactive_C.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

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

    my @no_match = qw( NO NO2 CO CH4 );
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
        $emissions{$VOC} += $rate->sum;
    }
    return \%emissions;
}
