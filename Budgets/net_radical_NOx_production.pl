#! /usr/bin/env perl
# Plot cumulative net radical to NOx production rates used in calculating NO emissions
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
my (%data, %families, %weights);

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $dir = "$base/$mechanism/${speciation}_Solvents_Only";
        my $boxmodel = "$dir/boxmodel";
        my $mecca = MECCA->new($boxmodel);
        my $eqn = "$dir/gas.eqn";
        my $kpp = KPP->new($eqn);
        my $radical_file = "$dir/radicals.txt";
        my @radicals = get_radicals($radical_file);
        $families{"radicals"} = [ @radicals ];
        $data{$mechanism}{$speciation} = get_data($kpp, $mecca);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(grid) `,
        q` library(ggthemes) `,
);
$R->run(q` data = data.frame() `);

foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    #$R->run(q` pre = data.frame(Mechanism = mechanism) `);
    foreach my $speciation (keys %{$data{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->set('production', $data{$mechanism}{$speciation});
        $R->run(q` pre = data.frame(Mechanism = mechanism) `,
                q` pre[speciation] = production `,
        );
    $R->run(q` pre = gather(pre, Speciation, Production, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` plot = ggplot(data, aes(x = Speciation, y = Production)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism, nrow = 1) `,
);

$R->run(q` CairoPDF(file = "radical_NOx_production_facet_mechanism.pdf", width = 8.7, height = 6) `,
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

sub get_data {
    my ($kpp, $mecca) = @_;
    $families{'NOx'} = [ qw( NO NO2 NO3 N2O5 ) ];
    $weights{'NOx'} = { N2O5 => 2 };
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
