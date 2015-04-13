#! /usr/bin/env perl
# Compare total radical production budgets in solvents only runs, facet by mechanism and speciation
# Version 0: Jane Coates 10/4/2015

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

foreach my $mechanism (@mechanisms) { 
    foreach my $speciation (@speciations) {
        my $boxmodel = "$base/$mechanism/${speciation}_Solvents_Only/boxmodel";
        my $mecca = MECCA->new($boxmodel);
        my $eqn  = "$base/$mechanism/${speciation}_Solvents_Only/gas.eqn";
        my $kpp = KPP->new($eqn);
        my $radical_file = "$base/$mechanism/${speciation}_Solvents_Only/radicals.txt";
        my @radicals = get_radicals($radical_file);
        $families{"radicals"} = [ @radicals ];
        $data{$mechanism}{$speciation} = get_data($mecca, $kpp, $mechanism, $speciation);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    $R->run(q` pre = data.frame(Mechanism = mechanism) `);
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->set('Total.radical.prod', $data{$mechanism}{$speciation});
        $R->run(q` pre[speciation] = Total.radical.prod `);
    }
    $R->run(q` pre = gather(pre, Speciation, Total.Radical.Prod, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` data$Speciation = factor(data$Speciation, levels = rev(c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08"))) `);

$R->run(q` plot = ggplot(data, aes(x = Speciation, y = Total.Radical.Prod)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism ) `,
        q` plot = plot + coord_flip() `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + theme(axis.title.y = element_blank()) `,
);

$R->run(q` CairoPDF(file = "total_radical_production_facet_mechanism.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = Total.Radical.Prod)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Speciation ) `,
        q` plot = plot + coord_flip() `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + theme(axis.title.y = element_blank()) `,
);

$R->run(q` CairoPDF(file = "total_radical_production_facet_speciation.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism, $speciation) = @_;
    my ($total_production_rates);
    $kpp->family({
            name    => "radicals",
            members => $families{"radicals"},
            weights => $weights{"radicals"},
    });
    my $producers = $kpp->producing("radicals");
    my $producer_yields = $kpp->effect_on("radicals", $producers);

    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        $total_production_rates += $rate->sum;
    }
    return $total_production_rates;
}

sub get_radicals {
    my ($file) = @_;
    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file : $!";
    local $/ = undef; 
    my @radicals = split /\s+/, <$in>;
    close $in;
    return @radicals;
}
