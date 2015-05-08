#! /usr/bin/env perl
# Compare cumulative HO2x production from C2H5OH in each mechanism and speciation
# Version 0: Jane Coates 8/5/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
my @speciations = qw( TNO IPCC EMEP DE94 UK98 UK08 );
my (%families, %weights, %data);
$families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];
my %mapping = ( ####check these and use csv
    MOZART  =>  {   
        TNO     => 0.676,
        IPCC    => 0.676,
        EMEP    => 1.000,
        DE94    => 0.495,
        UK98    => 0.630,
        UK08    => 0.685,
                },
    RADM2   =>  {
        TNO     => 0.183,
        IPCC    => 0.057,
        EMEP    => 0.273,
        DE94    => 0.334,
        UK98    => 0.448,
        UK08    => 0.327,
                },
);

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        if ($mechanism eq "MCM") {
            opendir DIR, "$base/$mechanism" or die $!;
            my @dirs = grep { $_ =~ /${speciation}_tagged/ } readdir DIR;
            closedir DIR;
            my $no_dirs = scalar @dirs;
            foreach my $directory (@dirs) {
                my $dir = "$base/$mechanism/$directory";
                my $mecca = MECCA->new("$dir/boxmodel");
                my $kpp = KPP->new("$dir/gas.eqn");
                $data{$mechanism}{$speciation}{$directory} = get_data($mecca, $kpp, $mechanism, $speciation, $no_dirs);
            }
        } else {
            my $dir = "$base/$mechanism/${speciation}_tagged_solvents_only";
            my $mecca = MECCA->new("$dir/boxmodel");
            my $kpp = KPP->new("$dir/gas.eqn");
            $data{$mechanism}{$speciation} = get_data($mecca, $kpp, $mechanism, $speciation);
        }
    }
}

foreach my $speciation (keys %{$data{"MCM"}}) {
    my %allocated;
    foreach my $dir (keys %{$data{"MCM"}{$speciation}}) {
        foreach my $reaction (keys %{$data{"MCM"}{$speciation}{$dir}}) {
            $allocated{$reaction} += $data{"MCM"}{$speciation}{$dir}{$reaction};
        }
    }
    $data{"MCM"}{$speciation} = \%allocated;
}

my $others_max = 1e7;
foreach my $mechanism (keys %data) {
    foreach my $speciation (keys %{$data{$mechanism}}) {
        foreach my $reaction (keys %{$data{$mechanism}{$speciation}}) {
            if ($data{$mechanism}{$speciation}{$reaction} < $others_max) {
                $data{$mechanism}{$speciation}{"Others"} += $data{$mechanism}{$speciation}{$reaction};
                delete $data{$mechanism}{$speciation}{$reaction};
            }
        }
        my @sorted = sort { $data{$mechanism}{$speciation}{$b} <=> $data{$mechanism}{$speciation}{$a} } keys %{$data{$mechanism}{$speciation}};
        my @final_sorted;
        foreach (@sorted) {
            next if ($_ =~ /Others/);
            push @final_sorted, { $_ => $data{$mechanism}{$speciation}{$_} };
        }
        push @final_sorted, { "Others" => $data{$mechanism}{$speciation}{"Others"} } if (defined $data{$mechanism}{$speciation}{"Others"});
        $data{$mechanism}{$speciation} = \@final_sorted;
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
        q` library(grid) `,
);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->run(q` pre = data.frame(Mechanism = mechanism, Speciation = speciation) `);
        foreach my $ref (@{$data{$mechanism}{$speciation}}) {
            foreach my $reaction (sort keys %$ref) {
                $R->set('reaction', $reaction);
                $R->set('ho2x', $ref->{$reaction});
                $R->run(q` pre[reaction] = ho2x `);
            }
        }
        $R->run(q` pre = gather(pre, Reaction, HO2x, -Mechanism, -Speciation) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = HO2x, fill = Reaction, order = Reaction)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2) `,
);

$R->run(q` CairoPDF(file = "C2H5OH_HO2x_production_budgets_facet_speciation.pdf", width = 10, height = 7) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism, $speciation, $no_dirs) = @_;
    $no_dirs = 1 unless (defined $no_dirs);
    my %production_rates;
    my $species = "HO2x";

    $kpp->family({
            name    => $species,
            members => $families{$species},
            weights => $weights{$species},
    });
    my $producers = $kpp->producing($species);
    my $producer_yields = $kpp->effect_on($species, $producers);
    print "No producers for $species\n" if (@$producers == 0);

    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        next unless ($reaction =~ /C2H5OH|HC3/);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        my ($reactants, $products) = split / = /, $kpp->reaction_string($reaction);
        $reactants =~ s/_(.*?)\b//g;
        $production_rates{$reactants} += $rate->sum;
    }

    if (exists $mapping{$mechanism}) {
        if (exists $mapping{$mechanism}{$speciation}) {
            $production_rates{$_} *= $mapping{$mechanism}{$speciation} foreach (keys %production_rates);
        }
    }
    $production_rates{"CO + OH"} /= $no_dirs if (defined $production_rates{"CO + OH"});
    return \%production_rates;
}
