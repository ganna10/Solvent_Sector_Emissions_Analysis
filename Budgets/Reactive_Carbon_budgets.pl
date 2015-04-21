#! /usr/bin/env perl
# Daily Reactive Carbon budgets
# Version 0: Jane Coates 21/4/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions/MCM";
#my @speciations = qw( EMEP TNO );
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
    $data{$speciation} = get_data($mecca, $kpp, $carbons);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` data = data.frame(Time) `);
foreach my $speciation (keys %data) {
    $R->set('speciation', $speciation);
    $R->set('reactive.c', [ map { $_ } $data{$speciation}->dog]);
    $R->run(q` data[speciation] = reactive.c `);
}
$R->run(q` data = gather(data, Speciation, Reactive.C, -Time) `);
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
#$R->run(q` data = data %>% group_by(Speciation) %>% summarize(Reactive.Carbon = sum(Carbon)) `);
##$R->run(q` data$Carbon = factor(data$Reactive.Carbon, levels = rev(c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))) `);
##$R->run(q` my.colours = c("12" = "#ba8b01", "11" = "#4c9383", "10" = "#ae4901", "9" = "#0352cb", "8" = "#6db875", "7" = "#0c3f78", "6" = "#b569b3", "5" = "#2b9eb3", "4" = "#ef6638", "3" = "#0e5628", "2" = "#f9c500", "1" = "#6c254f") `);
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` my.colours = c("TNO" = "#000000", "IPCC" = "#2b9eb3", "EMEP" = "#ef6638", "DE94" = "#f9c500", "GR95" = "#6c254f", "GR05" = "#0352cb", "UK98" = "#0e5c28", "UK08" = "#b569b3") `);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Reactive.C, colour = Speciation, group = Speciation)) `,
        q` plot = plot + geom_point() `,
        q` plot = plot + geom_line() `,
        #q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        #q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.ticks.x = element_blank()) `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.text.x = element_text(face = "bold", angle = 45, vjust = 0.7, hjust = 0.8)) `,
        q` plot = plot + theme(legend.position = c(1, 1)) `,
        q` plot = plot + theme(legend.justification = c(1, 1)) `,
        q` plot = plot + scale_colour_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "MCM_total_reactive_C.pdf", width = 8.7, height = 6) `,
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

sub get_data {
    my ($mecca, $kpp, $carbons) = @_; 
    my $ntime = $mecca->time->nelem;
    my $dt = $mecca->dt->at(0);
    my $n_per_day = 86400 / $dt;
    my $n_days = int $ntime / $n_per_day;
    
    my $reactive_C;
    my $all_reactions = $kpp->all_reactions();
    for (0..$#$all_reactions) {
        my $reaction = $all_reactions->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number);
        next if ($rate->sum == 0);
        my $reactant_C = 0;
        my $product_C = 0;
        my $reactants = $kpp->reactants($reaction);
        foreach my $reactant (@$reactants) {
            if (defined $carbons->{$reactant}) {
                $reactant_C += $carbons->{$reactant};
            } 
        }
        my $products = $kpp->products($reaction);
        foreach my $product (@$products) {
            if (defined $carbons->{$product}) {
                $product_C += $carbons->{$product};
            } 
        }
        my $net_reactive_C = $product_C - $reactant_C;
        next if ($net_reactive_C == 0);
        $reactive_C += $net_reactive_C * $rate(1:$ntime-2);
    }
    $reactive_C = $reactive_C->reshape($n_per_day, $n_days);
    $reactive_C = $reactive_C->sumover;
    return $reactive_C;
}
