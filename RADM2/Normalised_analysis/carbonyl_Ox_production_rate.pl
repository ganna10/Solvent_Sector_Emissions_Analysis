#! /usr/bin/env perl
# allocate carbonyl Ox production from tagged RADM2 runs to VOC 
# Version 0: Jane Coates 1/11/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base_dir = "/local/home/coates/Solvent_Emissions/RADM2";
opendir DIR, $base_dir or die "Can't open $base_dir";
my @tagged_runs = grep { $_ =~ /_tagged_/ } readdir DIR;
closedir DIR;

my $mecca = MECCA->new("$base_dir/DE94_tagged_solvents_only/boxmodel");
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0); 
my $n_per_day = 43200 / $dt;
my $n_days = int ($NTIME / $n_per_day);

my @alkanes = qw( ETH HC3 HC5 HC8 );
my @alkenes = qw( OL2 OLT OLI ISO );
my @aromatics = qw( TOL XYL );
my @carbonyls = qw( KET HCHO ALD );

my %families = ( 'HO2x' => [ qw( HO2 HO2NO2 ) ] );
my (%weights, %plot_data);

foreach my $run (@tagged_runs) {
    my ($label) = $run =~ /^(.*?)_/;
    my $model = "$base_dir/$run";
    my $boxmodel = "$model/boxmodel";
    my $eqn_file = "$model/gas.eqn";
    my $mecca = MECCA->new($boxmodel);
    my $kpp = KPP->new($eqn_file);
    my $ro2_file = "$model/RO2_species.txt";
    my @no2_reservoirs = get_no2_reservoirs($kpp, $ro2_file);
    $families{"Ox"} = [ qw( O3 O O1D NO2 HO2NO2 NO3 N2O5), @no2_reservoirs ];
    $weights{"Ox"} = { NO3 => 2, N2O5 => 3 };

    my (%production_rates, %consumption_rates, %emissions);
    foreach my $species (qw( Ox HO2x )) { 
        my ($producers, $producer_yields, $consumers, $consumer_yields);
        if (exists $families{$species}) { #get family reaction numbers and yields
            $kpp->family({ 
                            name    => $species,
                            members => $families{$species},
                            weights => $weights{$species},
            });
            $producers = $kpp->producing($species);
            $producer_yields = $kpp->effect_on($species, $producers);  
            $consumers = $kpp->consuming($species);
            $consumer_yields = $kpp->effect_on($species, $consumers);  
        } 
        die "No producers found for $species\n" if (@$producers == 0);#check that species reactions are found
        die "No consumers found for $species\n" if (@$consumers == 0);#check that species reactions are found
        
        for (0..$#$producers) { #get rates for all producing reactions
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $producer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            my $string;
            if (defined $parent) { # for tagged reactions
                if ($parent ~~ @carbonyls) {
                    $string = $parent;
                } elsif ($parent eq "CH4" or $parent ~~ @alkanes or $parent ~~ @alkenes or $parent ~~ @aromatics ) {
                    next;
                } else {
                    print "No group found for $parent\n";
                }
            } else { # for non-tagged reactions
                $string = $kpp->reaction_string($reaction);
            }
            $production_rates{$species}{$string} += $rate(1:$NTIME-2); #attribute rates to each parent tagged species and the non-tagged reactions
        }
        
        for (0..$#$consumers) { #get rates for all consuming reactions
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            my $string;
            if (defined $parent) { # for tagged reactions
                if ($parent ~~ @carbonyls) {
                    $string = $parent;
                } elsif ($parent eq "CH4" or $parent ~~ @alkanes or $parent ~~ @alkenes or $parent ~~ @aromatics ) {
                    next;
                } else {
                    print "No group found for $parent\n";
                }
            } else { # for non-tagged reactions
                $string = $kpp->reaction_string($reaction);
            }
            $consumption_rates{$species}{$string} += $rate(1:$NTIME-2); #attribute rates to each parent tagged species and the non-tagged reactions
        }
    } 
    remove_common_processes($production_rates{'HO2x'}, $consumption_rates{'HO2x'});
    my $ho2x_total_production = zeroes(PDL::float, $NTIME-2);
    $ho2x_total_production += $production_rates{'HO2x'}{$_} for (keys %{ $production_rates{'HO2x'} });

    foreach (keys %{ $production_rates{'HO2x'} }) {
        $production_rates{"Ox"}{$_} += $production_rates{"Ox"}{'HO2 + NO = NO2 + OH'} * $production_rates{'HO2x'}{$_} / $ho2x_total_production ;
        $consumption_rates{"Ox"}{$_} += $consumption_rates{"Ox"}{'HO2 + O3 = OH'} * $consumption_rates{'HO2x'}{$_} / $ho2x_total_production;
        $consumption_rates{"Ox"}{$_} += $consumption_rates{"Ox"}{'HO2 + NO3 = NO2 + OH'} * $consumption_rates{'HO2x'}{$_} / $ho2x_total_production;
    }
    delete $production_rates{"Ox"}{'HO2 + NO = NO2 + OH'};
    delete $consumption_rates{"Ox"}{'HO2 + O3 = OH'};
    delete $consumption_rates{"Ox"}{'HO2 + NO3 = NO2 + OH'};
    remove_common_processes($production_rates{'Ox'}, $consumption_rates{'Ox'});

    foreach my $parent (@carbonyls) {
        my $emission_reaction = $kpp->producing_from($parent, "UNITY");
        next if (@$emission_reaction == 0);
        my $reaction_nr = $kpp->reaction_number($emission_reaction->[0]);
        my $emission_rate = $mecca->rate($reaction_nr) * $dt;
        $emission_rate = $emission_rate(1:$NTIME-2);
        my $total_emissions = $emission_rate->sum; #integrated total emissions of parent VOC
        $emissions{$parent} += $total_emissions;
    } 

    foreach (keys %{$production_rates{"Ox"}}) {
        my $reshape = $production_rates{"Ox"}{$_}->copy->reshape($n_per_day, $n_days);
        my $integrate = $reshape->sumover;
        $integrate /= $emissions{$_};
        $integrate = $integrate(0:13:2);
        $plot_data{$label}{$_} = $integrate;
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `);
$R->run(q` library(reshape2) `);
$R->run(q` library(Cairo) `);
$R->run(q` library(grid) `);
$R->run(q` library(scales) `);
$R->run(q` library(plyr) `);

my @days = ("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7");
$R->set('time', [@days]);
$R->run(q` data = data.frame() `);

foreach my $speciation (keys %plot_data) {
    $R->run(q` pre = data.frame(time) `);
    foreach my $carbonyl (sort keys %{$plot_data{$speciation}}) {
        next if ($carbonyl =~ /\+/);
        $R->set('carbonyl', $carbonyl);
        $R->set('rate', [map { $_ } $plot_data{$speciation}{$carbonyl}->dog]);
        $R->run(q` pre[carbonyl] = rate `);
    }
    $R->set('speciation', $speciation);
    $R->run(q` pre$Speciation = rep(speciation, length(time)) `,
            q` pre = melt(pre, id.vars = c("time", "Speciation"), variable.name = "Carbonyl", value.name = "Rate") `,
            q` data = rbind(data, pre) `,
    );
}

$R->run(q` my.colours = c( "HCHO" = "#6c254f", "ALD" = "#0e5628", "KET" = "#f9c500") `,
        q` data$Carbonyl = factor(data$Carbonyl, levels = c("HCHO", "ALD", "KET")) `,
        q` scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `, #scientific label format for y-axis 
        q` data = ddply(data, .(Carbonyl)) `,
);
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` plot = ggplot(data = data, aes(x = time, y = Rate, fill = Carbonyl)) `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2)`,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + scale_y_continuous(label = scientific_10) `,
        q` plot = plot + ylab(expression(bold(paste(O[x], " Production Rate (molecules ", cm^-3, s^-1, ")")))) `,
        q` plot = plot + xlab("\n") `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(strip.text = element_text(size = 200, face = "bold")) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(panel.grid.major = element_blank()) `,
        q` plot = plot + theme(panel.grid.minor = element_blank()) `,
        q` plot = plot + theme(axis.text.x = element_text(size = 140)) `,
        q` plot = plot + theme(axis.text.y = element_text(size = 140))`,
        q` plot = plot + theme(axis.title.y = element_text(size = 200))`,
        q` plot = plot + theme(legend.title = element_blank(), legend.key.size = unit(7, "cm"), legend.text = element_text(size = 140, face = "bold"), legend.key = element_blank()) `, 
        q` plot = plot + scale_fill_manual(values = my.colours, limits = rev(levels(data$Carbonyl))) `,
);

$R->run(q` CairoPDF(file = "RADM2_carbonyls_Ox_budget.pdf", width = 200, height = 141) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop;

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

sub remove_common_processes {
    my ($production, $consumption) = @_;
    my %common_processes;
    $common_processes{$_} = 1 for (grep { defined $production->{$_} } keys %$consumption) ;

    foreach my $process (keys %common_processes) {
        my $net_effect = $production->{$process} + $consumption->{$process};
        if ($net_effect->sum > 0) { #if net production remove consumption processes, retain net production
            if (which($net_effect < 0)->nelem > 0) { #which gets indices of non-0 values, nelem get nr of elements
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
