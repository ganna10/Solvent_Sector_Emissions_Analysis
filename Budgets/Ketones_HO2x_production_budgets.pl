#! /usr/bin/env perl
# Daily Ketones HO2x production budgets
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
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my (%families, %weights, %data);
$families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];

my $mecca = MECCA->new("$base/RADM2/TNO_tagged_solvents_only/boxmodel");
my $ntime = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $n_per_day = 86400 / $dt;
my $n_days = int $ntime / $n_per_day;

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
my $sort_function = sub { $_[0]->sum } ;
foreach my $mechanism (keys %data) {
    foreach my $speciation (keys %{$data{$mechanism}}) {
        foreach my $reaction (keys %{$data{$mechanism}{$speciation}}) {
            if ($data{$mechanism}{$speciation}{$reaction}->sum < $others_max) {
                $data{$mechanism}{$speciation}{"Others"} += $data{$mechanism}{$speciation}{$reaction};
                delete $data{$mechanism}{$speciation}{$reaction};
            }
        }
        my @sorted = sort { &$sort_function($data{$mechanism}{$speciation}{$b}) <=> &$sort_function($data{$mechanism}{$speciation}{$a}) } keys %{$data{$mechanism}{$speciation}};
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
        q` library(dplyr) `,
);
$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->run(q` pre = data.frame(Time, Mechanism = rep(mechanism, length(Time)), Speciation = rep(speciation, length(Time))) `);
        foreach my $ref (@{$data{$mechanism}{$speciation}}) {
            foreach my $reaction (sort keys %$ref) {
                $R->set('reaction', $reaction);
                $R->set('ho2x', [ map { $_ } $ref->{$reaction}->dog ]);
                $R->run(q` pre[reaction] = ho2x `);
            }
        }
        $R->run(q` pre = gather(pre, Reaction, HO2x, -Mechanism, -Speciation, -Time) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` day.1 = filter(data, Time == "Day 1") `);
$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = HO2x, fill = Reaction)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_grid(Time ~ Speciation) `,
);

$R->run(q` CairoPDF(file = "Ketones_HO2x_production_budgets_facet_speciation.pdf", width = 10, height = 7) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism, $speciation, $no_dirs) = @_;
    $no_dirs = 1 unless (defined $no_dirs);
    my %production_rates;
    my $species = "HO2x";
    my @ketones = qw( CH3COCH3 MEK MIBK CYHEXONE HC5 KET );

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
        my ($number, $parent) = split /_/, $reaction;
        next unless (defined $parent and $parent ~~ @ketones);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        my ($reactants, $products) = split / = /, $kpp->reaction_string($reaction);
        $reactants =~ s/_(.*?)\b//g;
        $reactants =~ s/CH2O/HCHO/g;
        $reactants =~ s/MO2/CH3O2/g;
        $production_rates{$parent}{$reactants} += $rate;
    }

    my (%contributions, %final);
    if ($mechanism eq "RADM2" or $mechanism eq "MOZART") { #allocation of RADM2 and MOZART emitted species to Higher Alkanes
        open my $in, '<:encoding(utf-8)', "Ketones_RADM2_MOZ_fractional_contributions.csv" or die $!;
        my @lines = <$in>;
        close $in;
        foreach my $line (@lines) {
            chomp $line;
            my ($mechanism, $speciation, $species, $contribution) = split /,/, $line;
            $contributions{$mechanism}{$speciation}{$species} = $contribution;
        }

        foreach my $parent (keys %production_rates) {
            foreach my $reactants (keys %{$production_rates{$parent}}) {
                #print "$reactants : ", $production_rates{$parent}{$reactants}->sum, "\n";
                if (defined $contributions{$mechanism}{$speciation}{$parent}) {
                    #print $contributions{$mechanism}{$speciation}{$parent}, "\n";
                    $production_rates{$parent}{$reactants} *= $contributions{$mechanism}{$speciation}{$parent};
                } else {
                    delete $production_rates{$parent}{$reactants};
                } 
            }
        }
    }

    foreach my $parent (keys %production_rates) {
        foreach my $reaction (keys %{$production_rates{$parent}}) {
            my $reshape = $production_rates{$parent}{$reaction}->reshape($n_per_day, $n_days);
            $final{$reaction} += $reshape->sumover;
        }
    }
    $final{"CO + OH"} /= $no_dirs if (defined $final{"CO + OH"});
    return \%final;
}
