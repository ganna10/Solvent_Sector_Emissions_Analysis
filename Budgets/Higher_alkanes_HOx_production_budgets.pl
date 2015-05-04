#! /usr/bin/env perl
# Compare HOx production budgets from the higher alkanes between speciations and mechanisms, Solvent only runs
# Version 0: Jane Coates 29/4/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

#my $base = "/work/users/jco/Solvent_Emissions";
my $base = "/local/home/coates/Solvent_Emissions";
#my @mechanisms = qw( RADM2 );
my @mechanisms = qw( MCM MOZART RADM2 );
#my @speciations = qw( GR95 );
my @speciations = qw( TNO IPCC DE94 GR95 GR05 UK98 UK08 );
my (%families, %weights, %data);
$families{"HOx"} = [ qw( OH HO2 HO2NO2 ) ];

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) { 
        next if ($speciation eq "EMEP");
        if ($mechanism eq "MCM") { 
            opendir DIR, "$base/$mechanism" or die "Can't open dir : $!";
            my @dirs = grep { $_ =~ /${speciation}_tagged_solvents_only_alkanes/ } readdir DIR;
            closedir DIR;
            my $no_dirs = scalar @dirs;
            foreach my $directory (@dirs) {
                my $dir = "$base/$mechanism/$directory";
                my $boxmodel = "$dir/boxmodel";
                my $mecca = MECCA->new($boxmodel);
                my $eqn = "$dir/gas.eqn";
                my $kpp = KPP->new($eqn);
                $data{$mechanism}{$speciation}{$dir} = get_data($mecca, $kpp, $mechanism, $speciation, $directory, $no_dirs);
            }
        } else {
            my $dir = "$base/$mechanism/${speciation}_tagged_solvents_only";
            my $boxmodel = "$dir/boxmodel";
            my $mecca = MECCA->new($boxmodel);
            my $eqn = "$dir/gas.eqn";
            my $kpp = KPP->new($eqn);
            $data{$mechanism}{$speciation} = get_data($mecca, $kpp, $mechanism, $speciation);
        }
    }
}

foreach my $speciation (keys %{$data{"MCM"}}) {
    my %allocated;
    foreach my $dir (keys %{$data{"MCM"}{$speciation}}) {
        foreach my $item (sort keys %{$data{"MCM"}{$speciation}{$dir}}) {
            $allocated{$item} += $data{"MCM"}{$speciation}{$dir}{$item};
        }
    }
    $data{"MCM"}{$speciation} = \%allocated;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(ggthemes) `,
        q` library(Cairo) `,
        q` library(grid) `,
);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->run(q` pre = data.frame(Mechanism = mechanism, Speciation = speciation) `);
        foreach my $species (sort keys %{$data{$mechanism}{$speciation}}) {
            $R->set('species', $species);
            $R->set('hox', $data{$mechanism}{$speciation}{$species});
            $R->run(q` pre[species] = hox `);
        }
        $R->run(q` pre = gather(pre, Species, HOx, -Mechanism, -Speciation) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` data$Species = factor(data$Species, levels = c("NC6H14", "M2PE", "M3PE", "CHEX", "NC7H16", "M2HEX", "M3HEX", "NC8H18", "NC9H20", "NC10H22", "NC11H24", "NC12H26")) `);
$R->run(q` my.colours = c("NC6H14" = "#6c254f", "M2PE" = "#f9c500", "M3PE" = "#0e5c28", "CHEX" = "#ef6638", "NC7H16" = "#2b9eb3", "M2HEX" = "#b569b3", "M3HEX" = "#f7c56c", "NC8H18" = "#0352cb", "NC9H20" = "#ae4901", "NC10H22" = "#4c9383", "NC11H24" = "#8c1531", "NC12H26" = "#77aecc") `); 

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = HOx, fill = Species, order = Species)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + theme(axis.text.x = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.title.y = element_text(face = "bold")) `,
        q` plot = plot + theme(panel.margin = unit(5, "mm")) `,
        q` plot = plot + ylab("Total HOx Production (molecules cm-3)") `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + scale_fill_manual(values = my.colours, limits = rev(levels(data$Species))) `,
);

$R->run(q` CairoPDF(file = "Higher_alkanes_HOx_production_facet_speciation.pdf", width = 10, height = 7) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism, $speciation, $directory, $no_dirs) = @_;
    my $ntime = $mecca->time->nelem;
    $no_dirs = 1 unless (defined $no_dirs);
    my @higher_alkanes = qw( NC6H14 M2PE M3PE NC7H16 M2HEX M3HEX NC8H18 NC9H20 NC10H22 NC11H24 NC12H26 CHEX BIGALK HC5 HC8 );

    my (%production_rates, %consumption_rates);
    $kpp->family({
            name    => "HOx",
            members => $families{"HOx"},
            weights => $weights{"HOx"},
    });
    my $producers = $kpp->producing("HOx");
    my $producer_yields = $kpp->effect_on("HOx", $producers);
    my $consumers = $kpp->consuming("HOx");
    my $consumer_yields = $kpp->effect_on("HOx", $consumers);
    print "No consumers \n" if (@$consumers == 0);
    print "No producers \n" if (@$producers == 0);

    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0); 
        my ($number, $parent) = split /_/, $reaction;
        next unless ($parent ~~ @higher_alkanes);
        $production_rates{$parent} += $rate(1:$ntime-2);
    }

    for (0..$#$consumers) {
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
        next if ($rate->sum == 0); 
        my ($number, $parent) = split /_/, $reaction;
        next unless ($parent ~~ @higher_alkanes);
        $consumption_rates{$parent} += $rate(1:$ntime-2);
    }
    remove_common_processes(\%production_rates, \%consumption_rates);
    $production_rates{$_} = $production_rates{$_}->sum foreach (keys %production_rates);
    
    #print "before\n";
    #print "$_ : $production_rates{$_}\n" foreach keys %production_rates;
    my %contributions;
    if ($mechanism eq "RADM2" or $mechanism eq "MOZART") { #allocation of RADM2 and MOZART emitted species to MCM species that are Higher Alkanes
        open my $in, '<:encoding(utf-8)', "Higher_alkanes_contributions_of_MCM_species_in_MOZART_RADM2.csv" or die $!;
        my @lines = <$in>;
        close $in;
        foreach my $line (@lines) {
            chomp $line;
            my ($mechanism, $speciation, $mech_species, $mcm_species, $contribution) = split /,/, $line;
            $contributions{$mechanism}{$speciation}{$mech_species}{$mcm_species} = $contribution;
        }

        foreach my $species (keys %production_rates) {
            if (defined $contributions{$mechanism}{$speciation}{$species}) {
                foreach my $mcm_species (keys %{$contributions{$mechanism}{$speciation}{$species}}) {
                    $production_rates{$mcm_species} = $contributions{$mechanism}{$speciation}{$species}{$mcm_species} * $production_rates{$species};
                }
            }
            delete $production_rates{$species};
        }
    }
    #print "after\n";
    #print "$_ : $production_rates{$_}\n" foreach keys %production_rates;
#    if ($mechanism eq "RADM2" or $mechanism eq "MOZART") { #allocation of RADM2 and MOZART emitted species to Higher Alkanes
#        open my $in, '<:encoding(utf-8)', "Higher_alkanes_RADM2_MOZ_fractional_contributions.csv" or die $!;
#        my @lines = <$in>;
#        close $in;
#        foreach my $line (@lines) {
#            chomp $line;
#            my ($mechanism, $speciation, $species, $contribution) = split /,/, $line;
#            $contributions{$mechanism}{$speciation}{$species} = $contribution;
#        }
#
#        foreach my $species (keys %production_rates) {
#            if (defined $contributions{$mechanism}{$speciation}{$species}) {
#                $production_rates{$species} *= $contributions{$mechanism}{$speciation}{$species};
#            } else {
#                delete $production_rates{$species};
#            } 
#        }
#    }
    return \%production_rates;
}

sub remove_common_processes {
    my ($production, $consumption) = @_;
    my %common_processes;
    $common_processes{$_} = 1 for (grep { defined $production->{$_} } keys %$consumption) ;

    foreach my $process (keys %common_processes) {
        my $net_effect = $production->{$process} + $consumption->{$process};
        #print $process, $net_effect->nelem, "\n";
        if ($net_effect->sum > 0) { #if net production remove consumption processes, retain net production
            if (which($net_effect < 0)->nelem > 0) { #which gets indices of non-0 values, nelem get nr of elements
                #print "which if $process $net_effect\n";
                $production->{$process} .= $net_effect;
                $production->{$process}->where($net_effect < 0) .= 0; #where returns output pdl corresponding to input pdl and meeting the condition
                $consumption->{$process} .= $net_effect;
                $consumption->{$process}->where($net_effect > 0) .= 0; #where returns output pdl corresponding to input pdl and meeting the condition
                next;
            } 
            $production->{$process} = $net_effect;
            delete $consumption->{$process};
        } else { #net consumption
            if (which($net_effect > 0)->nelem > 0) {
                $production->{$process} .= $net_effect;
                $production->{$process}->where($net_effect < 0) .= 0;
                $consumption->{$process} .= $net_effect;
                $consumption->{$process}->where($net_effect > 0) .= 0;
                next;
            }
            $consumption->{$process} = $net_effect;
            delete $production->{$process};
        }
    }
} 
