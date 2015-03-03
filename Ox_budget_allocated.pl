#! /usr/bin/env perl
# Allocate Cumulative Ox production budgets back to categories used to defined speciations
# Jane Coates 1/3/2015

######################normalise by total emissions????

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
$families{"Ox"} = [ qw( O3 NO2 ) ];
$families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];

my $mecca = MECCA->new("$base/RADM2/TNO_Solvents_Only/boxmodel");
my $ntime = $mecca->time->nelem;

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my (%production_rates, %consumption_rates);
        if ($mechanism =~ /MOZART|RADM2/) { 
            my $boxmodel = "$base/$mechanism/${speciation}_tagged_solvents_only/boxmodel";
            my $mecca = MECCA->new($boxmodel);
            my $eqn_file = "$base/$mechanism/${speciation}_tagged_solvents_only/gas.eqn";
            my $kpp = KPP->new($eqn_file); 

            foreach my $species (qw( Ox HO2x )) {
                $kpp->family({
                        name    => $species,
                        members => $families{$species},
                        weights => $weights{$species},
                });
                my $producers = $kpp->producing($species);
                my $producer_yields = $kpp->effect_on($species, $producers);
                my $consumers = $kpp->consuming($species);
                my $consumer_yields = $kpp->effect_on($species, $consumers);

                print "No producers for $species in $mechanism, $speciation\n" if (@$producers == 0);
                print "No consumers for $species in $mechanism, $speciation\n" if (@$consumers == 0);

                for (0..$#$producers) {
                    my $reaction = $producers->[$_];
                    my $reaction_number = $kpp->reaction_number($reaction);
                    my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
                    next if ($rate->sum == 0);
                    my ($number, $parent) = split /_/, $reaction;
                    if (defined $parent) {
                        $production_rates{$species}{$parent} += $rate(1:$ntime-2)->sumover;
                    } else {
                        my $reaction_string = $kpp->reaction_string($reaction);
                        $production_rates{$species}{$reaction_string} += $rate(1:$ntime-2)->sumover;
                    }
                }

                for (0..$#$consumers) {
                    my $reaction = $consumers->[$_];
                    my $reaction_number = $kpp->reaction_number($reaction);
                    my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
                    next if ($rate->sum == 0);
                    my ($number, $parent) = split /_/, $reaction;
                    if (defined $parent) {
                        $consumption_rates{$species}{$parent} += $rate(1:$ntime-2)->sumover;
                    } else {
                        my $reaction_string = $kpp->reaction_string($reaction);
                        $consumption_rates{$species}{$reaction_string} += $rate(1:$ntime-2)->sumover;
                    }
                }
            }
        } elsif ($mechanism =~ /MCM/) {
            my $dir = "$base/$mechanism";
            opendir DIR, $dir or die "Can't open $dir : $!";
            my @runs = grep { $_ =~ /${speciation}_tagged/ } readdir DIR;
            closedir DIR;
            my $run_number = scalar @runs;
            foreach my $run (@runs) {
                my $boxmodel = "$base/$mechanism/$run/boxmodel";
                my $mecca = MECCA->new($boxmodel);
                my $eqn_file = "$base/$mechanism/$run/gas.eqn";
                my $kpp = KPP->new($eqn_file); 

                foreach my $species (qw( Ox HO2x )) {
                    $kpp->family({
                            name    => $species,
                            members => $families{$species},
                            weights => $weights{$species},
                    });
                    my $producers = $kpp->producing($species);
                    my $producer_yields = $kpp->effect_on($species, $producers);
                    my $consumers = $kpp->consuming($species);
                    my $consumer_yields = $kpp->effect_on($species, $consumers);

                    print "No producers for $species in $mechanism, $speciation and $run\n" if (@$producers == 0);
                    print "No consumers for $species in $mechanism, $speciation and $run\n" if (@$consumers == 0);

                    for (0..$#$producers) {
                        my $reaction = $producers->[$_];
                        next if ($reaction =~ /notag/);
                        my $reaction_number = $kpp->reaction_number($reaction);
                        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
                        next if ($rate->sum == 0);
                        my ($number, $parent) = split /_/, $reaction;
                        if (defined $parent) {
                            $production_rates{$species}{$parent} += $rate(1:$ntime-2)->sumover;
                        } else {
                            my $reaction_string = $kpp->reaction_string($reaction);
                            next if ($reaction_string =~ /notag/);
                            $production_rates{$species}{$reaction_string} += $rate(1:$ntime-2)->sumover / $run_number;
                        }
                    }

                    for (0..$#$consumers) {
                        my $reaction = $consumers->[$_];
                        next if ($reaction =~ /notag/);
                        my $reaction_number = $kpp->reaction_number($reaction);
                        my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
                        next if ($rate->sum == 0);
                        my ($number, $parent) = split /_/, $reaction;
                        if (defined $parent) {
                            $consumption_rates{$species}{$parent} += $rate(1:$ntime-2)->sumover;
                        } else {
                            my $reaction_string = $kpp->reaction_string($reaction);
                            next if ($reaction_string =~ /notag/);
                            $consumption_rates{$species}{$reaction_string} += $rate(1:$ntime-2)->sumover / $run_number;
                        }
                    }
                }
            }
        }

        remove_common_processes($production_rates{"HO2x"}, $consumption_rates{"HO2x"});
        my $total_ho2x_production;
        $total_ho2x_production += $production_rates{"HO2x"}{$_} foreach (keys %{$production_rates{"HO2x"}});
        
        foreach my $reaction (keys %{$production_rates{'HO2x'}}) {
            $production_rates{"Ox"}{$reaction} += $production_rates{"Ox"}{'HO2 + NO = NO2 + OH'} * $production_rates{'HO2x'}{$reaction} / $total_ho2x_production;
            $consumption_rates{"Ox"}{$reaction} += $consumption_rates{"Ox"}{'HO2 + O3 = OH'} * $consumption_rates{'HO2x'}{$reaction} / $total_ho2x_production; 
        }
        delete $production_rates{"Ox"}{'HO2 + NO = NO2 + OH'};
        delete $consumption_rates{"Ox"}{'HO2 + O3 = OH'};
        remove_common_processes($production_rates{"Ox"}, $consumption_rates{"Ox"});

        foreach my $process (sort keys %{$production_rates{"Ox"}}) {
            my $category = get_category($process, $mechanism, $speciation);
            $data{$mechanism}{$speciation}{$category} += $production_rates{"Ox"}{$process};
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
$R->run(q` my.colours = c("Acids" = "#cc6329", "Alcohols" = "#6c254f", "Benzene" = "#898989", "Butanes" = "#77aecc", "Chlorinated HC's" = "#f9c500", "Esters" = "#623812", "Ethane" = "#86b650", "Ethene" = "#f36a71", "Ethers" = "#ba8b01", "Ethyne" = "#dc3522", "Formaldehyde" = "#9bb18d", "Higher alkanes" = "#0e5c28", "Ketones" = "#ef6638", "Other aldehydes" = "#8ed6d2", "Other alkenes, alkynes, dienes" = "#58691b", "Other aromatics" = "#b569b3", "Others" = "#2b9eb3", "Pentanes" = "#8c1531", "Propane" = "#e7e85e", "Propene" = "#0c3f78", "Terpenes" = "#ae4901", "Toluene" = "#0352cb", "Trimethylbenzenes" = "#c9a415", "Xylenes" = "#1b695b", "CO" = "#6d6537", "Methane" = "#be2448", "Inorganic" = "#000000") `);

$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    #$R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        my $sum = 0;
        #$R->run(q` pre = data.frame(Mechanism = mechanism) `);
        #$R->set('speciation', $speciation);
        #$R->run(q` pre$Speciation = speciation `);
        foreach my $type (sort keys %{$data{$mechanism}{$speciation}}) {
            $sum += $data{$mechanism}{$speciation}{$type};
            #next if ($type eq "Inorganic"); ####removing inorganic as it is such a large term
            #next if ($type eq "CO"); ####Contribution practically the same throughout
            #next if ($type eq "Methane"); ###focusing on specified VOC
            #$R->set('type', $type);
            #$R->set('Ox.production', $data{$mechanism}{$speciation}{$type}->at(0));
            #$R->run(q` pre[type] = Ox.production `);
        }
        print "$mechanism: $speciation total Ox production = $sum\n";
        #$R->run(q` pre = gather(pre, Type, Ox.Production, -Mechanism, -Speciation) `,
        #q` data = rbind(data, pre) `,
        #);
        #my $p = $R->run(q` print(pre) `);
        #print $p, "\n";
    }
}
#$R->set('filename', "Ox_Allocated_tagged_solvents_only.pdf");
#$R->set('title', "Cumulative Ox Production Budget");
#$R->run(q` data$Type = factor(data$Type, levels = c("Acids", "Alcohols", "Benzene", "Butanes", "Chlorinated HC's", "Esters", "Ethane", "Ethene", "Ethers", "Ethyne", "Formaldehyde", "Higher alkanes", "Ketones", "Other aldehydes", "Other alkenes, alkynes, dienes", "Other aromatics", "Pentanes", "Propane", "Propene", "Terpenes", "Toluene", "Trimethylbenzenes", "Xylenes")) `);
#$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);

#$R->run(q` plot = ggplot(data, aes(x = Speciation, y = Ox.Production, fill = Type)) `,
#        q` plot = plot + geom_bar(stat = "identity", position = "stack") `,
#        q` plot = plot + facet_wrap( ~ Mechanism) `,
#        q` plot = plot + theme_tufte() `,
#        q` plot = plot + ggtitle(title) `,
#        q` plot = plot + ylab("Ox Production (molecules cm-3)") `,
#        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
#        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
#        q` plot = plot + theme(axis.ticks.x = element_blank()) `,
#        q` plot = plot + theme(legend.title = element_blank()) `,
#        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
#        q` plot = plot + theme(plot.title = element_text(face = "bold")) `,
#        q` plot = plot + theme(axis.title.y = element_text(face = "bold")) `,
#        q` plot = plot + theme(axis.title.x = element_blank()) `,
#        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
#        q` plot = plot + theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 0.7, vjust = 1.0)) `,
#        q` plot = plot + theme(panel.margin = unit(5, "mm")) `,
#        q` plot = plot + scale_fill_manual(values = my.colours, limits = rev(levels(data$Type))) `,
#);
#
#$R->run(q` CairoPDF(file = filename, width = 8.6, height = 6) `,
#        q` print(plot) `,
#        q` dev.off() `,
#);

$R->stop();

sub get_category {
    my ($process, $mechanism, $speciation) = @_;
    my $category;
    if ($process =~ /CO \+ OH/) {
        $category = "CO";
    } elsif ($process =~ /=/) {
        $category = "Inorganic";
    } elsif ($process eq "CH4") {
        $category = "Methane";
    } elsif ($process =~ /^TOL/) {
        $category = "Toluene";
    } elsif ($process =~ /XYL/) {
        $category = "Xylenes";
    } elsif ($process =~ /^TM/) {
        $category = "Trimethylbenzenes";
    } elsif ($process eq "BENZENE") {
        $category = "Benzene";
    } elsif ($process =~ /BENZ|DIME|STYRENE|TOL$/) {
        $category = "Other aromatics";
    } elsif ($process eq "ETH" or $process eq "C2H6") {
        $category = "Ethane";
    } elsif ($process eq "OL2" or $process eq "C2H4") {
        $category = "Ethene";
    } elsif ($process =~ /NC4|IC4/) {
        $category = "Butanes";
    } elsif ($process eq "C2H2") {
        $category = "Ethyne";
    } elsif ($process eq "HCHO" or $process eq "CH2O") {
        $category = "Formaldehyde";
    } elsif ($process eq "HC8" or $process =~ /NC6|NC7|NC8|NC9|NC1\d|M\dPE|M\dHEX|CHEX/) {
        $category = "Higher alkanes";
    } elsif ($process =~ /CH3OH|C2H5OH|IPROPOL|NBUTOL|NPROPOL|BUT2OL|IBUTOL|MIBKAOH|C6H5CH2OH|ETHGLY|PROPGLY/) {
        $category = "Alcohols";
    } elsif ($process =~ /CH3OCH3|BUOX2ETOH|PR2OHMOX|EOX2EOL|MO2EOL/) {
        $category = "Ethers";
    } elsif ($process =~ /ACET$/) {
        $category = "Esters";
    } elsif ($process =~ /TRICLETH|CH3CCL3|CH2CL2|TCE/) {
        $category = "Chlorinated HC's";
    } elsif ($process =~ /HCOOH|CH3COOH|PROPACID|CH3CO2H/) {
        $category = "Acids";
    } elsif ($process eq "KET" or $process =~ /MEK|CH3COCH3|MIBK|CYHEXONE/) {
        $category = "Ketones";
    } elsif ($process eq "ALD" or $process =~ /MACR|CH3CHO|C2H5CHO|C3H7CHO|IPRCHO|C4H9CHO|ACR|C4ALDB/) {
        $category = "Other aldehydes";
    } elsif ($process eq "OLI" or $process =~ /C5H8|ISO|BIGENE|BUT1ENE/) {
        $category = "Other alkenes, alkynes, dienes";
    } elsif ($process eq "HC5" or $process eq "BIGALK" or $process =~ /NC5|IC5|NEOP/) {
        $category = "Pentanes";
    } elsif ($process eq "HC3" or $process eq "C3H8") {
        $category = "Propane";
    } elsif ($process eq "OLT" or $process eq "C3H6") {
        $category = "Propene";
    } elsif ($process =~ /C10H16|LIMONENE|PINENE/) {
        $category = "Terpenes";
    } else {
        print "No category found for $process in $mechanism, $speciation\n";
    }
    return $category;
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
