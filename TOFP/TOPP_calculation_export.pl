#! /usr/bin/env perl
# Calculate cumulative TOPP of emitted VOC in each mechanism
# Version 0: Jane Coates

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
my (%data, %families, %weights, %carbons, %mapping);
$families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];
$weights{"Ox"} = { NO3 => 2, N2O5 => 3 };

my $mcm_carbons_file = "$base/MCM/EMEP_tagged_solvents_only_all/carbons.txt";
$carbons{"MCM"} = mcm_n_carbon($mcm_carbons_file);
my $mozart_carbons_file = "$base/MOZART/EMEP_tagged_solvents_only/carbons.txt";
$carbons{"MOZART"} = mozart_n_carbon($mozart_carbons_file);
my $radm2_carbons_file = "$base/RADM2/EMEP_tagged_solvents_only/carbons.txt";
$carbons{"RADM2"} = carbons_others($radm2_carbons_file);
my $mapping_file = "Mapping_MOZ_RADM2_of_MCM_species.csv";
open my $in, '<:encoding(utf-8)', $mapping_file or die "Can't open $mapping_file : $!";
my @lines = <$in>;
close $in;

foreach my $line (@lines) {
    chomp $line;
    my ($mechanism, $speciation, $mapping) = split /,/,$line;
    my ($species, $mcm_species) = split /:/, $mapping;
    my @mcm_species = split /\s/, $mcm_species;
    $mapping{$mechanism}{$speciation}{$species} = \@mcm_species;
}

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        opendir DIR, "$base/$mechanism" or die "Can't open dir : $!";
        my @dirs = grep { $_ =~ /${speciation}_tagged/ } readdir DIR;
        closedir DIR;
        foreach my $directory (@dirs) {
            my $dir = "$base/$mechanism/$directory";
            my $mecca = MECCA->new("$dir/boxmodel");
            my $kpp = KPP->new("$dir/gas.eqn");
            my $ro2_file = "$dir/RO2_species.txt";
            my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
            $families{"Ox"} = [ qw( O3 O O1D NO2 NO3 N2O5 HO2NO2 ), @no2_reservoirs ];
            $data{$mechanism}{$speciation}{$directory} = get_data($kpp, $mecca, $mechanism, $speciation);
        }
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
);

foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        $R->set('speciation', $speciation);
        foreach my $dir (sort keys %{$data{$mechanism}{$speciation}}) {
            $R->set('dir', $dir);
            $R->run(q` data = data.frame(Speciation = speciation) `);
            foreach my $VOC (sort keys %{$data{$mechanism}{$speciation}{$dir}}) {
                $R->set('voc', $VOC);
                $R->set('topp', $data{$mechanism}{$speciation}{$dir}{$VOC});
                $R->run(q` data[voc] = topp `);
            }
            $R->run(q` data = gather(data, VOC, TOPP, -Speciation) `,
                    q` data = select(data, VOC, TOPP) `,
            );
            $R->set('TOPP.file', "${mechanism}_${speciation}_${dir}_cumulative_TOPP_values.csv");
            $R->run(q` write.csv(data, file = TOPP.file, quote = FALSE, row.names = FALSE) `);

            #$R->set('TOPP.plot', "${mechanism}_${speciation}_${dir}_cumulative_TOPP_plot.pdf");
            #my $p = $R->run(q` print(data) `);
            #print $p, "\n";

#            $R->run(q` plot = ggplot(data, aes(x = VOC, y = TOPP)) `,
#                    q` plot = plot + geom_point(size = 4) `,
#                    q` plot = plot + theme_tufte() `,
#                    q` plot = plot + theme(axis.title.x = element_blank()) `,
#                    q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
#            );
#
#            $R->run(q` CairoPDF(file = TOPP.plot, width = 8.7, height = 6) `,
#                    q` print(plot) `,
#                    q` dev.off() `,
#            );
        }
    }
}

$R->stop();

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

sub get_no2_reservoirs { #get species that are produced when radical species react with NO2
    my ($kpp, $file) = @_; 
    open my $in, '<:encoding(utf-8)', $file or die $!; 
    my @ro2;
    for (<$in>) {
        push @ro2, split /\s+/, $_; 
    }
    close $in;
    my @no2_reservoirs;
    foreach my $ro2 (@ro2) {
        my ($reactions) = $kpp->reacting_with($ro2, 'NO2');
        foreach my $reaction (@$reactions) {
            my ($products) = $kpp->products($reaction);
            if (@$products == 1) {
                push @no2_reservoirs, $products->[0];
            }   
        }   
    }   
    return @no2_reservoirs;
} 

sub mcm_n_carbon {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die $!;
    my @lines = (<$in>);
    close $in;
    my %carbons;
    foreach my $line (@lines) {
        chomp $line;
        my ($species, $smile) = split ' ', $line;
        my $C_number = 0;
        if ($smile =~ /\./) {
            $C_number = 8;
        } else {
            $C_number++ while ($smile =~ m/C/gi);
        }
        $carbons{$species} = $C_number;
    }
    return \%carbons;
}

sub mozart_n_carbon {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die $!;
    my $words = join ',', (<$in>);
    close $in;
    my ($string) = $words =~ /Solution(.*?)End\sSolution/s;
    $string =~ s/^\s+,//;
    $string =~ s/,\s+$//;
    $string =~ s/\s+/ /g;
    $string =~ s/RO2(.*?)->//;
    $string =~ s/ROOH(.*?)->//;
    my @species = split ',', $string;
    my %carbons;
    foreach my $species (@species) {
        $species =~ s/^\s+|\s+$//g;
        my $C_number = 0;
        if ($species !~ /->/ and $species !~ /(C[0-9])/) {
            $C_number ++ while ($species =~ m/C/g);
            $carbons{$species} = $C_number;
        } elsif ($species !~ /->/ and $species =~ /(C[0-9])/) { 
            my ($c_nr) = $species =~ /(C[0-9]+)/s;
            $c_nr =~ s/C//; 
            $C_number = $c_nr;
            $carbons{$species} = $C_number;
        } else {
            my ($mech, $molecule) = split ' -> ', $species;
            $mech =~ s/^\s+|\s+$//g;
            if ($molecule =~ /(C[0-9]+)/) { 
                my ($c_nr) = $molecule =~ /(C[0-9]+)/s;
                $c_nr =~ s/C//; 
                $C_number = $c_nr;
                $carbons{$mech} = $C_number;
            } else {
                $C_number ++ while ($molecule =~ m/C/g);
                $carbons{$mech} = $C_number;
            }
        }
    } 
    return \%carbons;
}

sub carbons_others { #get C-number from file names that have species and C# separated by space
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die "Cannot open file $file: $!";
    my (@lines) = (<$in>);
    close $in;
    my %carbons;
    foreach my $line (@lines) {
        my ($species, $C_number) = split '\s', $line;
        $carbons{$species} = $C_number;
    }
    return \%carbons;
}

sub get_data {
    my ($kpp, $mecca, $mechanism, $speciation) = @_;
    my $ntime = $mecca->time->nelem;
    my (%TOPP, %production_rates, %consumption_rates);
    
    foreach my $species ( qw( Ox HO2x ) ) {
        $kpp->family({
                name    => $species,
                members => $families{$species},
                weights => $weights{$species},
        });
        my $producers = $kpp->producing($species);
        my $producer_yields = $kpp->effect_on($species, $producers);
        my $consumers = $kpp->consuming($species);
        my $consumer_yields = $kpp->effect_on($species, $consumers);
        print "No consumers for $species\n" if (@$consumers == 0);
        print "No producers for $species\n" if (@$producers == 0);

        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            next if ($rate->sum == 0);
            my ($number, $parent) = split /_/, $reaction;
            if (defined $parent) {
                $production_rates{$species}{$parent} += $rate(1:$ntime-2);
            } else {
                my ($reactants, $products) = split / = /, $kpp->reaction_string($reaction);
                $production_rates{$species}{$reactants} += $rate(1:$ntime-2);
            }
        }

        for (0..$#$consumers) {
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
            next if ($rate->sum == 0);
            my ($number, $parent) = split /_/, $reaction;
            if (defined $parent) {
                $consumption_rates{$species}{$parent} += $rate(1:$ntime-2);
            } else {
                my ($reactants, $products) = split / = /, $kpp->reaction_string($reaction);
                $consumption_rates{$species}{$reactants} += $rate(1:$ntime-2);
            }
        }
    }
    remove_common_processes($production_rates{"HO2x"}, $consumption_rates{"HO2x"});
    my $total_ho2x;
    $total_ho2x += $production_rates{"HO2x"}{$_} foreach (keys $production_rates{"HO2x"});

    foreach my $process (keys %{$production_rates{"HO2x"}}) {
        $production_rates{"Ox"}{$process} += $production_rates{"Ox"}{"HO2 + NO"} * $production_rates{"HO2x"}{$process} / $total_ho2x ;
        $consumption_rates{"Ox"}{$process} += $consumption_rates{"Ox"}{"HO2 + O3"} * $consumption_rates{"HO2x"}{$process} / $total_ho2x;
        $consumption_rates{"Ox"}{$process} += $consumption_rates{"Ox"}{"HO2 + NO3"} * $consumption_rates{"HO2x"}{$process} / $total_ho2x;
    }
    delete $production_rates{"Ox"}{"HO2 + NO"};
    delete $consumption_rates{"Ox"}{"HO2 + O3"};
    delete $consumption_rates{"Ox"}{"HO2 + NO3"};
    remove_common_processes($production_rates{"Ox"}, $consumption_rates{"Ox"});

    foreach my $VOC (keys %{$production_rates{"Ox"}}) {
        next if ($VOC =~ /\+/ or $VOC eq "CH4" or $VOC eq "notag");
        my $emission_reaction = $kpp->producing_from($VOC, "UNITY");
        my $reaction_number = $kpp->reaction_number($emission_reaction->[0]);
        my $emission_rates = $mecca->rate($reaction_number);
        $emission_rates = $emission_rates->sum;
        $TOPP{$VOC} = $production_rates{"Ox"}{$VOC}->sum / $emission_rates;
    }
    if ($mechanism =~ /MOZ|RADM/) {
        my @moz_match = qw( TOLUENE C2H5OH MEK MACR CH3CHO );
        my @delete = qw( CH2O C10H16 BIGENE BIGALK ISOP HC3 HC5 HC8 KET ISO TOL XYL ETH OL2 OLT OLI ALD ORA1 ORA2 );
        foreach my $VOC (keys %TOPP) {
            print "No mapping for $VOC\n" unless (defined $mapping{$mechanism}{$speciation}{$VOC});
            #print "$VOC : $mapping{$mechanism}{$speciation}{$VOC}\n";
            foreach my $mcm (@{$mapping{$mechanism}{$speciation}{$VOC}}) {
                next if ($mechanism eq "MOZART" and $mcm ~~ @moz_match);
                #print "$VOC : $mcm : $carbons{$mechanism}{$VOC} and $carbons{'MCM'}{$mcm}\n";
                $TOPP{$mcm} = $TOPP{$VOC} * $carbons{'MCM'}{$mcm} / $carbons{$mechanism}{$VOC};
            }
        }
        delete $TOPP{$_} foreach (@delete);
    }
    #print "$_ : $TOPP{$_}\n" foreach keys %TOPP;
    return \%TOPP;
}
