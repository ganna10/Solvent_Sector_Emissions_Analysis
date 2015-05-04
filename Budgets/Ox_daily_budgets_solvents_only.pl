#! /usr/bin/env perl
# Compare daily total Ox budget from solvent only runs between mechanisms, facet Time and speciation
# Version 0: Jane Coates 29/4/2015

use strict;
use diagnostics;
use MECCA;
use KPP;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my (%families, %weights, %data);

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $dir = "$base/$mechanism/${speciation}_Solvents_Only";
        my $mecca = MECCA->new("$dir/boxmodel");
        my $kpp = KPP->new("$dir/gas.eqn");
        my $RO2_file = "$dir/RO2_species.txt";
        my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
        $families{"Ox"} = [ qw( O3 O O1D NO2 NO3 N2O5 HO2NO2 ), @no2_reservoirs ];
        $weights{"Ox"} = { NO3 => 2, N2O5 => 3 };
        $data{$mechanism}{$speciation} = get_data($kpp, $mecca);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
);
$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    $R->run(q` pre = data.frame(Time, Mechanism = rep(mechanism, length(Time))) `);
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->set('ox', [ map { $_ } $data{$mechanism}{$speciation}->dog ]);
        $R->run(q` pre[speciation] = ox `);
    }
    $R->run(q` pre = gather(pre, Speciation, Ox, -Time, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = Ox)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_grid(Time ~ Speciation) `,
);

$R->run(q` CairoPDF(file = "Daily_Ox_budgets_Solvents_only.pdf", width = 10, height = 7) `,
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

sub get_data {
    my ($kpp, $mecca) = @_;
    my $ntime = $mecca->time->nelem;
    my $dt = $mecca->dt->at(0);
    my $n_per_day = 86400 / $dt;
    my $n_days = int $ntime / $n_per_day;

    my $production;
    $kpp->family({
            name    => "Ox",
            members => $families{"Ox"},
            weights => $weights{"Ox"},
    });
    my $producers = $kpp->producing("Ox");
    my $producer_yields = $kpp->effect_on("Ox", $producers);
    print "No producers\n" if (@$producers == 0);

    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        $production += $rate;
    }
    $production = $production->reshape($n_per_day, $n_days);
    $production = $production->sumover;
    return $production;
}
