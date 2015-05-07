#! /usr/bin/env perl
#Compare OH production budgets, both total and only those from organic sources
#Version 0: Jane Coates 6/5/2015

use strict;
use diagnostics;
use PDL;
use PDL::NiceSlice;
use KPP;
use MECCA;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
#my @speciations = qw( TNO );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my (%data);

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        if ($mechanism eq "MCM") {
            opendir DIR, "$base/$mechanism" or die $!;
            my @dirs = grep { $_ =~ /${speciation}_tagged/ } readdir DIR;
            close DIR;
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
        foreach my $item (keys %{$data{"MCM"}{$speciation}{$dir}}) {
            $allocated{$item} += $data{"MCM"}{$speciation}{$dir}{$item};
        }
    }
    $data{"MCM"}{$speciation} = \%allocated;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(ggthemes) `,
        q` library(grid) `,
);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->run(q` pre = data.frame(Mechanism = mechanism, Speciation = speciation) `);
        foreach my $parent (sort keys %{$data{$mechanism}{$speciation}}) {
            $R->set('item', $parent);
            $R->set('oh', $data{$mechanism}{$speciation}{$parent});
            $R->run(q` pre[item] = oh `);
        }
        $R->run(q` pre = gather(pre, Item, OH, -Mechanism, -Speciation) `,
                q` data = rbind(data, pre) `,
        );
    }
}

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = OH, fill = Item)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
);

$R->run(q` CairoPDF(file = "Total_OH_production_all_items_facet_speciation.pdf", width = 10, height = 7) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->run(q` organic = data %>% filter(!grepl("=", Item)) `);
#my $p = $R->run(q` print(organic) `);
#print $p, "\n";

$R->run(q` plot = ggplot(organic, aes(x = Mechanism, y = OH, fill = Item)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
);

$R->run(q` CairoPDF(file = "Total_OH_production_organic_facet_speciation.pdf", width = 10, height = 7) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism, $speciation, $no_dirs) = @_;
    $no_dirs = 1 unless (defined $no_dirs);

    my $species = "OH";
    my %production_rate;
    my $producers = $kpp->producing($species);
    my $producer_yields = $kpp->effect_on($species, $producers);
    print "No producers for $species\n" if (@$producers == 0);

    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        my ($number, $parent) = split /_/, $reaction;
        if (defined $parent) {
            $production_rate{$parent} += $rate->sum / $no_dirs;
        } else {
            my $reaction_string = $kpp->reaction_string($reaction);
            $production_rate{$reaction_string} += $rate->sum / $no_dirs;
        }
    }
    delete $production_rate{"notag"};
    return \%production_rate;
}
