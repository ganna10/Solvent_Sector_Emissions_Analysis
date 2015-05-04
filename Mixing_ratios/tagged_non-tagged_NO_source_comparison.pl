#! /usr/bin/env perl
#Compare NO emission time series between tagged and non-tagged runs in each mechanism and speciation
#Verion 0: Jane Coates 29/4/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw( MOZART MCM RADM2 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my @runs = qw( Solvents_Only tagged_solvents_only );
my %data;

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        foreach my $run (@runs) {
            if ($mechanism eq "MCM" and $run =~ /tagged/) {
                opendir DIR, "$base/$mechanism" or die $!;
                my @dirs = grep { $_ =~ /${speciation}_tagged/ } readdir DIR;
                closedir DIR;
                my $no_dirs = scalar @dirs;
                foreach my $directory (@dirs) {
                    my $dir = "$base/$mechanism/$directory";
                    my $mecca = MECCA->new("$dir/boxmodel");
                    my $kpp = KPP->new("$dir/gas.eqn");
                    $data{$mechanism}{$speciation}{$run}{$directory} += get_data($kpp, $mecca) / $no_dirs;
                }
            } else {
                my $dir = "$base/$mechanism/${speciation}_$run";
                my $mecca = MECCA->new("$dir/boxmodel");
                my $kpp = KPP->new("$dir/gas.eqn");
                $data{$mechanism}{$speciation}{$run} = get_data($kpp, $mecca);
            }
        }
    }
}

foreach my $speciation (keys %{$data{"MCM"}}) {
    my $final = 0;
    $final += $data{"MCM"}{$speciation}{"tagged_solvents_only"}{$_} foreach (keys %{$data{"MCM"}{$speciation}{"tagged_solvents_only"}});
    $data{"MCM"}{$speciation}{"tagged_solvents_only"} = $final;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
        q` library(ggthemes) `,
);

my $mecca = MECCA->new("$base/RADM2/TNO_Solvents_Only/boxmodel");
my $time = $mecca->time;
$time -= $time->at(0);
$time /= 86400;

$R->set('Time', [ map { $_ } $time->dog ]);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->run(q` pre = data.frame(Time, Mechanism = rep(mechanism, length(Time)), Speciation = rep(speciation, length(Time))) `);
        foreach my $run (sort keys %{$data{$mechanism}{$speciation}}) {
            $R->set('run', $run);
            $R->set('rate', [ map { $_ } $data{$mechanism}{$speciation}{$run}->dog ]);
            $R->run(q` pre[run] = rate `);
        }
        $R->run(q` pre = gather(pre, Run, Rate, -Time, -Mechanism, -Speciation) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, colour = Run)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + facet_grid(Mechanism ~ Speciation) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
);

$R->run(q` CairoPDF(file = "NO_source_tagged_vs_non-tagged.pdf", width = 10, height = 7) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($kpp, $mecca) = @_;
    my $reaction = $kpp->producing_from("NO", "UNITY");
    my $reaction_number = $kpp->reaction_number($reaction->[0]);
    my $rate = $mecca->rate($reaction_number);
    return $rate;
}
