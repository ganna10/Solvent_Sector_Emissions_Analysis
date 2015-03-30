#! /usr/bin/env perl
# Compare total Ox production in each mechanism, speciation in tagged and non-tagged solvent only runs
# Version 0: Jane Coates 12/3/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

#my $base = "/work/users/jco/Solvent_Emissions";
my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my @run_type = qw( Tagged Non-Tagged );

my $mecca = MECCA->new("$base/RADM2/TNO_Solvents_Only/boxmodel");
my $ntime = $mecca->time->nelem;
my %families = (
    "Ox"    => [ qw( O3 NO2 ) ],
);
my (%weights, %data);

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        foreach my $run_type (@run_type) {
            $data{$mechanism}{$speciation}{$run_type} = get_data($mechanism, $speciation, $run_type);
        }
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
        $R->run(q` pre = data.frame(Mechanism = mechanism) `);
        $R->set('speciation', $speciation);
        foreach my $run_type (sort keys %{$data{$mechanism}{$speciation}}) {
            $R->set('run', $run_type);
            $R->set('total.prod', [$data{$mechanism}{$speciation}{$run_type}->at(0)]);
            $R->run(q` pre[run] = total.prod `);
        }
        $R->run(q` pre$Speciation = speciation `);
        $R->run(q` pre = gather(pre, Run, Total.Prod, -Mechanism, -Speciation) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` my.colours = c("Non-Tagged" = "#2b9eb3", "Tagged" = "#ef6638") `);

$R->run(q` plot = ggplot(data, aes(x = Speciation, y = Total.Prod, fill = Run)) `,
        q` plot = plot + geom_bar(stat = "identity", position = "dodge") `,
        q` plot = plot + facet_wrap( ~ Mechanism ) `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + ggtitle("Total Ox Production after 7 Days") `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + ylab("Total Ox Production (molecules cm-3)") `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
        q` plot = plot + theme(strip.text = element_text(size = 22, face = "bold")) `,
        q` plot = plot + theme(plot.title = element_text(size = 24, face = "bold")) `,
        q` plot = plot + theme(axis.title = element_text(size = 17, face = "bold")) `,
        q` plot = plot + theme(axis.text = element_text(size = 15)) `,
        q` plot = plot + theme(legend.text = element_text(size = 15)) `,
        q` plot = plot + theme(legend.key.size = unit(7, "mm")) `,
        q` plot = plot + theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 0.8, vjust = 0.9)) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.ticks.x = element_blank()) `,
        q` plot = plot + guides(guide = guide_legend(direction = "horizontal")) `,
        q` plot = plot + theme(legend.position = "bottom") `,
        q` plot = plot + theme(panel.margin = unit(5, "mm")) `,
        q` plot = plot + scale_fill_manual(values = my.colours) `,
);

$R->run(q` Cairo(file = "Ox_production_tagged_vs_non_tagged.pdf", type = "pdf", bg = "transparent", unit = "cm", width = 27, height = 20.3) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mechanism, $speciation, $run_type) = @_;
    my $total_production;
    if ($mechanism =~ /MCM|MOZ|RA/ and $run_type =~ /Non-/) {
        my $boxmodel = "$base/$mechanism/${speciation}_Solvents_Only/boxmodel";
        my $mecca = MECCA->new($boxmodel);
        my $eqn = "$base/$mechanism/${speciation}_Solvents_Only/gas.eqn";
        my $kpp = KPP->new($eqn);
        $total_production = get_Ox_production($mecca, $kpp);
    } elsif ($mechanism =~ /MOZ|RA/) { #tagged run for MOZART and RADM2 
        my $boxmodel = "$base/$mechanism/${speciation}_tagged_solvents_only/boxmodel";
        my $mecca = MECCA->new($boxmodel);
        my $eqn = "$base/$mechanism/${speciation}_tagged_solvents_only/gas.eqn";
        my $kpp = KPP->new($eqn);
        $total_production = get_Ox_production($mecca, $kpp);
    } elsif ($mechanism eq "MCM") {
        my $directory = "$base/$mechanism";
        opendir DIR, $directory or die "Can't open $directory : $!";
        my @dirs = grep { $_ =~ /${speciation}_tagged/ } readdir DIR;
        close DIR;
        $total_production = 0;
        my $number_of_runs = scalar @dirs;
        foreach my $dir (@dirs) {
            my $boxmodel = "$directory/$dir/boxmodel";
            my $mecca = MECCA->new($boxmodel);
            my $eqn = "$directory/$dir/gas.eqn";
            my $kpp = KPP->new($eqn);
            $total_production += get_Ox_production($mecca, $kpp);
        }
        $total_production /= $number_of_runs;
    } else {
        print "No data possible for $mechanism, $speciation and $run_type\n";
    }
    return $total_production;
}

sub get_Ox_production {
    my ($mecca, $kpp) = @_;
    my $production_rates = 0;
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
        $production_rates += $rate(1:$ntime-2);
    }
    return $production_rates->sumover; 
}
