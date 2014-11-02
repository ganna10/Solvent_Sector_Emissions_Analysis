#! /usr/bin/env perl
# allocate alkenes Ox production from tagged MCM runs
# Version 0: Jane Coates 2/11/2014

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base_dir = "/local/home/coates/Solvent_Emissions/MCM";
opendir DIR, $base_dir or die "Can't open $base_dir";
my @tagged_runs = grep { $_ =~ /_tagged_/ } readdir DIR;
closedir DIR;
#my @tagged_runs = qw( EMEP_tagged_solvents_only_all );

my $mecca = MECCA->new("$base_dir/EMEP_tagged_solvents_only_all/boxmodel");
my $NTIME = $mecca->time->nelem;
my $dt = $mecca->dt->at(0); 
my $n_per_day = 43200 / $dt;
my $n_days = int ($NTIME / $n_per_day);

my @alkanes = qw(C2H6 C3H8 CHEX IC4H10 IC5H12 M2HEX M2PE M3HEX M3PE NC10H22 NC11H24 NC12H26 NC4H10 NC5H12 NC6H14 NC7H16 NC8H18 NC9H20 NEOP );
my @alkenes = qw(APINENE BPINENE BUT1ENE C2H4 C3H6 C5H8 LIMONENE );
my @aromatics = qw(BENZENE DIME35EB EBENZ IPBENZ METHTOL MXYL OXYL PBENZ PETHTOL PXYL STYRENE TM123B TM124B TM135B TOLUENE );
my @carbonyls = qw(ACR C2H5CHO C3H7CHO C4ALDB C4H9CHO CH3CHO CH3COCH3 CYHEXONE HCHO IPRCHO MACR MEK MIBK );
my @alcohols = qw(BUT2OL C2H5OH C6H5CH2OH ETHGLY IBUTOL IPROPOL MIBKAOH NBUTOL NPROPOL PROPGLY );
my @acids = qw(CH3CO2H CH3OH HCOOH PROPACID );
my @esters = qw(ETHACET IPROACET NBUTACET NPROACET );
my @ethers = qw(BUOX2ETOH CH3OCH3 EOX2EOL MO2EOL PR2OHMOX );
my @chlorinated = qw( CH2CL2 CH3CCL3 TCE TRICLETH );
my @alkynes = qw( C2H2 ); 

my %families = ( 'HO2x' => [ qw( HO2 HO2NO2 ) ] );
my (%weights, %production_plot_rates, %production_rates, %consumption_rates, %plot_data);

foreach my $run (@tagged_runs) {
    next unless (-d "$base_dir/$run");
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
    my $ntime = $mecca->time->nelem;

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
            next if ($kpp->reaction_string($reaction) =~ /notag/);
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            my $string;
            if (defined $parent) { # for tagged reactions
                if ($parent ~~ @alkenes) {
                    $string = $parent;
                } elsif ($parent ~~ @alkanes or $parent ~~ @aromatics or $parent ~~ @alcohols or $parent ~~ @carbonyls or $parent ~~ @acids or $parent ~~ @esters or $parent ~~ @ethers or $parent ~~ @alkynes or $parent ~~ @chlorinated or $parent eq "CH4") {
                    next;
                } else {
                    print "Nothing found for $parent\n";
                }
            } else { # for non-tagged reactions
                $string = $kpp->reaction_string($reaction);
            }
            $production_rates{$label}{$species}{$string} += $rate(1:$ntime-2); #attribute rates to each parent tagged species and the non-tagged reactions
        }
        
        for (0..$#$consumers) { #get rates for all consuming reactions
            my $reaction = $consumers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $consumer_yields->[$_] * $mecca->rate($reaction_number); 
            next if ($kpp->reaction_string($reaction) =~ /notag/);
            next if ($rate->sum == 0); # do not include reactions that do not occur 
            my ($r_number, $parent) = split /_/, $reaction; #remove tag from reaction number
            my $string;
            if (defined $parent) { # for tagged reactions
                if ($parent ~~ @alkenes) {
                    $string = $parent;
                } elsif ($parent ~~ @alkanes or $parent ~~ @aromatics or $parent ~~ @alcohols or $parent ~~ @carbonyls or $parent ~~ @acids or $parent ~~ @esters or $parent ~~ @ethers or $parent ~~ @alkynes or $parent ~~ @chlorinated or $parent eq "CH4") {
                    next;
                } else {
                    print "Nothing found for $parent\n";
                }
            } else { # for non-tagged reactions
                $string = $kpp->reaction_string($reaction);
            }
            $consumption_rates{$label}{$species}{$string} += $rate(1:$ntime-2); #attribute rates to each parent tagged species and the non-tagged reactions
        }
    } 
    remove_common_processes($production_rates{$label}{'HO2x'}, $consumption_rates{$label}{'HO2x'});
    my $ho2x_total_production = zeroes(PDL::float, $ntime-2);
    $ho2x_total_production += $production_rates{$label}{'HO2x'}{$_} for (keys %{ $production_rates{$label}{'HO2x'} });

    foreach (keys %{ $production_rates{$label}{'HO2x'} }) {
        $production_rates{$label}{"Ox"}{$_} += $production_rates{$label}{"Ox"}{'HO2 + NO = NO2 + OH'} * $production_rates{$label}{'HO2x'}{$_} / $ho2x_total_production ;
        $consumption_rates{$label}{"Ox"}{$_} += $consumption_rates{$label}{"Ox"}{'HO2 + O3 = OH'} * $consumption_rates{$label}{'HO2x'}{$_} / $ho2x_total_production;
        $consumption_rates{$label}{"Ox"}{$_} += $consumption_rates{$label}{"Ox"}{'HO2 + NO3 = NO2 + OH'} * $consumption_rates{$label}{'HO2x'}{$_} / $ho2x_total_production;
    }
    delete $production_rates{$label}{"Ox"}{'HO2 + NO = NO2 + OH'};
    delete $consumption_rates{$label}{"Ox"}{'HO2 + O3 = OH'};
    delete $consumption_rates{$label}{"Ox"}{'HO2 + NO3 = NO2 + OH'};
    remove_common_processes($production_rates{$label}{'Ox'}, $consumption_rates{$label}{'Ox'});
}

foreach my $run (sort keys %production_rates) { 
    foreach my $reaction (sort keys %{$production_rates{$run}{'Ox'}}) { #average out data that comes from all the tagged runs such as CH4, CO etc
        if ($reaction eq 'CH4') {
            $production_rates{$run}{'Ox'}{$reaction} /= get_number_of_tagged_runs($run);
        } elsif ($reaction =~ /\+/) {
            delete $production_rates{$run}{"Ox"}{$reaction};
        }
    }
}

foreach my $run (sort keys %production_rates) {
    foreach (sort keys %{$production_rates{$run}{"Ox"}}) {
        my $reshape = $production_rates{$run}{"Ox"}{$_}->copy->reshape($n_per_day, $n_days);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $plot_data{$run}{$_} = $integrate;
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `);
$R->run(q` library(reshape2) `);
$R->run(q` library(scales) `);
$R->run(q` library(Cairo) `);
$R->run(q` library(grid) `);

my @days = ("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7");
$R->set('Time', [@days]);
$R->run(q` data = data.frame() `);
foreach my $run (keys %plot_data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $alkene (sort keys %{$plot_data{$run}}) {
        $R->set('alkene', $alkene);
        $R->set('rate', [map { $_ } $plot_data{$run}{$alkene}->dog]);
        $R->run(q` pre[alkene] = rate `);
    }
    $R->set('speciation', $run);
    $R->run(q` pre$Speciation = rep(speciation, length(Time)) `,
            q` pre = melt(pre, id.vars = c("Time", "Speciation"), variable.name = "Alkene", value.name = "Rate") `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print "$p\n";
$R->run(q` scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) } `, #scientific label format for y-axis
        q` my.colours = c( "C2H4" = "#6c254f", "C3H6" = "#f9c500", "BUT1ENE" = "#0e5628", "C5H8" = "#ef6638", "LIMONENE" = "#2b9eb3", "APINENE" = "#b569b3", "BPINENE" = "#0c3f78") `,
        q` my.names = c("C2H4", "C3H6", "BUT1ENE", "C5H8", "APINENE", "BPINENE", "LIMONENE" ) `,
);

$R->run(q` plot = ggplot(data = data, aes(x = Time, y = Rate, fill = Alkene)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2)`,
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
        q` plot = plot + scale_fill_manual(values = my.colours, limits = rev(my.names)) `,
);

$R->run(q` CairoPDF(file = "MCM_Ox_budget_alkenes.pdf", width = 200, height = 141) `,
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

sub get_number_of_tagged_runs {
    my ($speciation) = @_;
    my $number_of_runs = 0;
    foreach my $run (@tagged_runs) {
        $number_of_runs++ if ($run =~ $speciation) ;
    }
    return $number_of_runs;
}
