#! /usr/bin/env perl
# Compare production budgets of HCHO in tagged and non-tagged runs.
# Version 0: Jane Coates 30/3/2015

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
my @runs = qw( Solvents_Only tagged_solvents_only );
my (%families, %weights, %data);

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        foreach my $run (@runs) {
            my ($boxmodel, $mecca, $eqn, $kpp, $species, $dir);
            if ($mechanism eq "MCM" and $run =~ /tagged/) {
                opendir DIR, "$base/$mechanism" or die "Can't open dir : $!";
                my @dirs = grep { $_ =~ /${speciation}_$run/ } readdir DIR;
                closedir DIR;
                my $no_dirs = scalar @dirs;
                my $data;
                foreach my $directory (@dirs) {
                    $dir = "$base/$mechanism/$directory";
                    $boxmodel = "$dir/boxmodel";
                    $mecca = MECCA->new($boxmodel);
                    $eqn = "$dir/gas.eqn";
                    $kpp = KPP->new($eqn);
                    $data += get_data($mecca, $kpp, "HCHO", $mechanism, $speciation, $run, $dir);
                }
                $data{$mechanism}{$speciation}{$run} = $data / $no_dirs;
            } else {
                $boxmodel = "$base/$mechanism/${speciation}_$run/boxmodel";
                $mecca = MECCA->new($boxmodel);
                $eqn = "$base/$mechanism/${speciation}_$run/gas.eqn";
                $kpp = KPP->new($eqn);
                if ($mechanism eq "MOZART") {
                    $species = "CH2O";
                } else {
                    $species = "HCHO";
                }
                $data{$mechanism}{$speciation}{$run} = get_data($mecca, $kpp, $species, $mechanism, $speciation, $run);
            }
        }
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
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
            $R->set('total.prod', $data{$mechanism}{$speciation}{$run_type});
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
#$R->run(q` my.colours = c("Non-Tagged" = "#2b9eb3", "Tagged" = "#ef6638") `);

$R->run(q` plot = ggplot(data, aes(x = Speciation, y = Total.Prod, fill = Run)) `,
        q` plot = plot + geom_bar(stat = "identity", position = "dodge") `,
        q` plot = plot + facet_wrap( ~ Mechanism ) `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + ggtitle("Total HCHO Production after 7 Days") `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + ylab("Total HCHO Production (molecules cm-3)") `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
        q` plot = plot + theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 0.8, vjust = 0.9)) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.ticks.x = element_blank()) `,
        q` plot = plot + guides(guide = guide_legend(direction = "horizontal")) `,
        q` plot = plot + theme(legend.position = "bottom") `,
        #q` plot = plot + scale_fill_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "HCHO_production_tagged_vs_non_tagged.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $species, $mechanism, $speciation, $run, $dir) = @_;
    my $ntime = $mecca->time->nelem;

    my ($production_rate, $consumption_rate, $producers, $producer_yields, $consumers, $consumer_yields);
    if ($run =~ /tagged/) {
        my $spc;
        if ($mechanism eq "MCM") {
            $spc = "$dir/gas.spc";
        } else {
            $spc = "$base/$mechanism/${speciation}_$run/gas.spc";
        }
        my $tagged_species = get_tagged_species($spc, $species);
        $families{"tagged_$species"} = [ @$tagged_species ];
        $kpp->family({
                name    => "tagged_$species",
                members => $families{"tagged_$species"},
                weights => $weights{"tagged_$species"},
        });
        $producers = $kpp->producing("tagged_$species");
        $producer_yields = $kpp->effect_on("tagged_$species", $producers); 
        $consumers = $kpp->consuming("tagged_$species");
        $consumer_yields = $kpp->effect_on("tagged_$species", $consumers); 
    } else {
        $producers = $kpp->producing($species);
        $producer_yields = $kpp->effect_on($species, $producers); 
        $consumers = $kpp->consuming($species);
        $consumer_yields = $kpp->effect_on($species, $consumers); 
    }
    print "No producers found in $run, $speciation and $mechanism\n" if (@$producers == 0);
    print "No consumers found in $run, $speciation and $mechanism\n" if (@$consumers == 0);

    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0); 
        $production_rate += $rate(1:$ntime-2);
    }

    for (0..$#$consumers) {
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
        next if ($rate->sum == 0); 
        $consumption_rate += $rate(1:$ntime-2);
    }

    my $net = $production_rate + $consumption_rate;
    return $net->sum;
}

sub get_tagged_species {
    my ($spc, $lookup) = @_;
    my @matches;
    open my $file, '<:encoding(utf-8)', $spc or die "Can't open $spc\n";
    local $/ = undef;
    my @lines = split /\n/, <$file>;
    close $file;
    foreach my $line (@lines) {
        next unless ($line =~ /^$lookup/);
        #next if ($line =~ /_notag/);
        $line =~ s/ = IGNORE.*$//;
        push @matches, $line;
    }
    return \@matches;
}
