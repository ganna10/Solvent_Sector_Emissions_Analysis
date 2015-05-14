#! /usr/bin/env perl
# Allocate First day and Cumulative Ox production budgets back to categories used to defined speciations, both absolute and percent contribution
# Version 0: Jane Coates 13/5/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

#my $base = "/local/home/coates/Solvent_Emissions";
my $base = "/work/users/jco/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
#my @speciations = qw( IPCC );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my (%families, %weights, %data);
$families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];

my %category_mapping = (
    MOZART  =>  {
        TNO     => {    BIGALK  => [ '0.130 Butanes', '0.165 Chlorinated', '0.215 Esters', '0.078 Ethers', '0.411 Higher_alkanes' ],
                        TOLUENE => [ '0.050 Other_aromatics', '0.494 Toluene', '0.026 Trimethylbenzenes', '0.430 Xylenes' ],
                    },
        IPCC    => {    BIGALK  => [ '0.341 Butanes', '0.008 Chlorinated', '0.111 Esters', '0.033 Ethers', '0.271 Higher_alkanes', '0.237 Pentanes' ],
                        TOLUENE => [ '0.244 Benzene', '0.164 Other_aromatics', '0.339 Toluene', '0.037 Trimethylbenzenes', '0.216 Xylenes' ],
                    },
        EMEP  =>    {   BIGALK  => [ '1.0 Butanes' ],
                        TOLUENE => [ '1.0 Xylenes' ],
                    },
        DE94    =>  {   BIGALK  => [ '0.277 Butanes', '0.037 Chlorinated', '0.200 Esters', '0.023 Ethers', '0.462 Higher_alkanes' ],
                        TOLUENE => [ '0.259 Other_aromatics', '0.337 Toluene', '0.134 Trimethylbenzenes', '0.270 Xylenes' ],
                    },
        GR95    =>  {   BIGALK  => [ '0.104 Chlorinated', '0.205 Esters', '0.691 Higher_alkanes' ],
                        TOLUENE => [ '0.009 Benzene', '0.558 Other_aromatics', '0.298 Toluene', '0.086 Trimethylbenzenes', '0.050 Xylenes' ],
                    },
        GR05    =>  {   BIGALK  => [ '0.092 Chlorinated', '0.352 Esters', '0.556 Higher_alkanes' ],
                        TOLUENE => [ '0.691 Other_aromatics', '0.196 Toluene', '0.066 Trimethylbenzenes', '0.047 Xylenes' ],
                    },
        UK98    =>  {   BIGALK  => [ '0.116 Butanes', '0.211 Chlorinated', '0.161 Esters', '0.089 Ethers', '0.317 Higher_alkanes', '0.107 Pentanes' ],
                        TOLUENE => [ '0.281 Other_aromatics', '0.384 Toluene', '0.093 Trimethylbenzenes', '0.242 Xylenes' ],
                    },
        UK08    =>  {   BIGALK  => [ '0.158 Butanes', '0.072 Chlorinated', '0.121 Esters', '0.099 Ethers', '0.547 Higher_alkanes' ],
                        TOLUENE => [ '0.000 Benzene', '0.336 Other_aromatics', '0.192 Toluene', '0.157 Trimethylbenzenes', '0.314 Xylenes' ],
                    },
                },
    RADM2   =>  {
        TNO     => {    HC3     => [ '0.362 Alcohols', '0.175 Butanes', '0.222 Chlorinated', '0.190 Esters', '0.007 Ethers', '0.045 Propane' ],
                        HC5     => [ '0.108 Alcohols', '0.275 Esters', '0.461 Higher_alkanes', '0.156 Ketones' ],
                        HC8     => [ '0.086 Alcohols', '0.186 Ethers', '0.728 Higher_alkanes' ],
                        KET     => [ '0.011 Alcohols', '0.989 Ketones' ],
                        TOL     => [ '0.057 Other_aromatics', '0.943 Toluene' ],
                        XYL     => [ '0.042 Other_aromatics', '0.054 Trimethylbenzenes', '0.903 Xylenes' ],
                    },
        IPCC    =>  {   HC3     => [ '0.112 Alcohols', '0.363 Butanes', '0.008 Chlorinated', '0.108 Esters', '0.002 Ethers', '0.164 Other_alkenes,_alkynes,_dienes', '0.242 Propane'] ,
                        HC5     => [ '0.032 Alcohols', '0.029 Esters', '0.232 Higher_alkanes', '0.025 Ketones', '0.681 Pentanes' ],
                        HC8     => [ '0.057 Alcohols', '0.130 Ethers', '0.813 Higher_alkanes' ],
                        KET     => [ '0.021 Alcohols', '0.979 Ketones' ],
                        TOL     => [ '0.358 Benzene', '0.144 Other_aromatics', '0.498 Toluene' ],
                        XYL     => [ '0.208 Other_aromatics', '0.116 Trimethylbenzenes', '0.676 Xylenes' ],
                    },
        EMEP    =>  {   HC3     => [ '0.419 Alcohols', '0.581 Butanes' ],
                    }, 
        DE94    =>  {   HC3     => [ '0.661 Alcohols', '0.181 Butanes', '0.086 Esters', '0.001 Ethers', '0.047 Propane' ],
                        HC5     => [ '0.317 Alcohols', '0.201 Esters', '0.402 Higher_alkanes', '0.080 Ketones' ],
                        HC8     => [ '0.539 Alcohols', '0.028 Ethers', '0.433 Higher_alkanes' ],
                        KET     => [ '0.063 Alcohols', '0.937 Ketones' ],
                        TOL     => [ '0.314 Other_aromatics', '0.686 Toluene' ],
                        XYL     => [ '0.205 Other_aromatics', '0.264 Trimethylbenzenes', '0.531 Xylenes' ],
                    },
        GR95    =>  {   HC3     => [ '0.753 Alcohols', '0.093 Chlorinated', '0.155 Esters' ],
                        HC5     => [ '0.508 Esters', '0.492 Ketones' ],
                        HC8     => [ '0.086 Alcohols', '0.914 Higher_alkanes' ],
                        OLT     => [ '1.0 Other_alkenes,_alkynes,_dienes'],
                        TOL     => [ '0.014 Benzene', '0.526 Other_aromatics', '0.460 Toluene' ],
                        XYL     => [ '0.615 Other_aromatics', '0.244 Trimethylbenzenes', '0.141 Xylenes' ],
                    },
        GR05    =>  {   HC3     => [ '0.853 Alcohols', '0.039 Chlorinated', '0.107 Esters' ],
                        HC5     => [ '0.932 Esters', '0.068 Ketones' ],
                        HC8     => [ '0.402 Alcohols', '0.598 Higher_alkanes' ],
                        KET     => [ '0.007 Alcohols', '0.993 Ketones' ],
                        OLT     => [ '1.0 Other_alkenes,_alkynes,_dienes'],
                        TOL     => [ '0.728 Other_aromatics', '0.272 Toluene' ],
                        XYL     => [ '0.596 Other_aromatics', '0.236 Trimethylbenzenes', '0.168 Xylenes' ],
                    },
        UK98    =>  {   HC3     => [ '0.492 Alcohols', '0.103 Butanes', '0.188 Chlorinated', '0.088 Esters', '0.005 Ethers', '0.124 Propane' ],
                        HC5     => [ '0.354 Alcohols', '0.154 Esters', '0.134 Higher_alkanes', '0.091 Ketones', '0.267 Pentanes' ] ,
                        HC8     => [ '0.296 Alcohols', '0.168 Ethers', '0.536 Higher_alkanes' ],
                        KET     => [ '0.059 Alcohols', '0.941 Ketones' ],
                        TOL     => [ '0.292 Other_aromatics', '0.708 Toluene' ],
                        XYL     => [ '0.268 Other_aromatics', '0.203 Trimethylbenzenes', '0.529 Xylenes' ],
                    },
        UK08    =>  {   HC3     => [ '0.763 Alcohols', '0.118 Butanes', '0.054 Chlorinated', '0.027 Esters', '0.002 Ethers', '0.000 Other_alkenes,_alkynes,_dienes', '0.036 Propane' ],
                        HC5     => [ '0.413 Alcohols', '0.220 Esters', '0.209 Higher_alkanes', '0.147 Ketones', '0.010 Pentanes' ],
                        HC8     => [ '0.283 Alcohols', '0.123 Ethers', '0.594 Higher_alkanes' ],
                        KET     => [ '0.078 Alcohols', '0.922 Ketones' ],
                        TOL     => [ '0.000 Benzene', '0.479 Other_aromatics', '0.521 Toluene' ],
                        XYL     => [ '0.253 Other_aromatics', '0.249 Trimethylbenzenes', '0.498 Xylenes' ],
                    },
                }
);

my $mecca = MECCA->new("$base/RADM2/TNO_Solvents_Only/boxmodel");
my $ntime = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $n_per_day = 86400 / $dt;
my $n_days = int $ntime / $n_per_day;

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my (%production_rates, %consumption_rates);
        if ($mechanism =~ /MOZART|RADM2/) { 
            my $boxmodel = "$base/$mechanism/${speciation}_tagged_solvents_only/boxmodel";
            my $mecca = MECCA->new($boxmodel);
            my $eqn_file = "$base/$mechanism/${speciation}_tagged_solvents_only/gas.eqn";
            my $kpp = KPP->new($eqn_file); 
            my $RO2_file = "$base/$mechanism/${speciation}_tagged_solvents_only/RO2_species.txt";
            my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
            $families{"Ox"} = [ qw( O3 NO2 O O1D NO3 N2O5 HO2NO2 ), @no2_reservoirs ];
            $weights{"Ox"} = { NO3 => 2, N2O5 => 3 };
            $data{$mechanism}{$speciation} = get_data($mecca, $kpp, $mechanism, $speciation);
        } else {
            my $mcm_base = "$base/$mechanism";
            opendir DIR, $mcm_base or die "can't open $mcm_base : $!";
            my @runs = grep { $_ =~ /${speciation}_tagged_solvents_only/ } readdir DIR;
            closedir DIR;
            my $no_runs = scalar @runs;
            foreach my $run (@runs) {
                my $boxmodel = "$base/$mechanism/$run/boxmodel";
                my $mecca = MECCA->new($boxmodel);
                my $eqn = "$base/$mechanism/$run/gas.eqn";
                my $kpp = KPP->new($eqn);
                my $RO2_file = "$base/$mechanism/$run/RO2_species.txt";
                my @no2_reservoirs = get_no2_reservoirs($kpp, $RO2_file);
                $families{"Ox"} = [ qw( O3 NO2 O O1D NO3 N2O5 HO2NO2 ), @no2_reservoirs ];
                $weights{"Ox"} = { NO3 => 2, N2O5 => 3 };
                $data{$mechanism}{$speciation}{$run} = get_data($mecca, $kpp, $mechanism, $speciation, $run, $no_runs);
            }
        }
    }
}

foreach my $speciation (keys %{$data{"MCM"}}) {
    my %allocated;
    foreach my $run (keys %{$data{'MCM'}{$speciation}}) {
        foreach my $process (sort keys %{$data{'MCM'}{$speciation}{$run}}) {
            $allocated{$process} += $data{'MCM'}{$speciation}{$run}{$process};
        }
    }
    $data{'MCM'}{$speciation} = \%allocated;
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
        q` library(grid) `,
        q` library(scales) `,
);

$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    #print "$mechanism\n";
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        #print "\t$speciation\n";
        $R->set('speciation', $speciation);
        foreach my $process (sort keys %{$data{$mechanism}{$speciation}}) {
            #print "\t\t$process $data{$mechanism}{$speciation}{$process}\n";
            $R->set('process', $process);
            $R->set('first.day', $data{$mechanism}{$speciation}{$process}->at(0));
            $R->set('cumulative', $data{$mechanism}{$speciation}{$process}->sum);
            $R->run(q` pre = data.frame(Process = process, First.Day = first.day, Cumulative = cumulative, Mechanism = mechanism, Speciation = speciation)  `,
                    q` data = rbind(data, pre) `,
            );
        }
    }
}
$R->run(q` data = gather(data, Time, Ox.Production, -Process, -Mechanism, -Speciation) `);
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
#$R->run(q` write.table(data, file = "Ox_production_allocated.csv", sep = ",", col.names = FALSE, row.names = FALSE, quote = FALSE) `);

$R->run(q` my.colours = c("Acids" = "#0352cb", "Alcohols" = "#011e4b", "Benzene" = "#0e5c28", "Butanes" = "#6e254f", "Chlorinated" = "#ef6638", "CO" = "#898989", "Esters" = "#0357d8", "Ethane" = "#a15daf", "Ethene" = "#dfb100", "Ethers" = "#2b9eb3", "Higher alkanes" = "#38103d", "Ketones" = "#02388b", "Aldehydes" = "#0499fc", "Other alkenes, alkynes, dienes" = "#f9c500", "Other aromatics" = "#216105", "Others" = "#2b9eb3", "Pentanes" = "#8f2ac9", "Propane" = "#5c1b54", "Propene" = "#796000", "Terpenes" = "#b99300", "Toluene" = "#0b6956", "Trimethylbenzenes" = "#109c00", "Xylenes" = "#0c734b", "CO" = "#6d6537", "Methane" = "#77aecc", "Inorganic" = "#000000") `);
#$R->run(q` my.colours = c("Acids" = "#cc6329", "Alcohols" = "#6c254f", "Benzene" = "#8c6238", "Butanes" = "#86b650", "Chlorinated" = "#f9c500", "CO" = "#898989", "Esters" = "#f3aa7f", "Ethane" = "#77aecc", "Ethene" = "#1c3e3d", "Ethers" = "#ba8b01", "Higher alkanes" = "#0e5c28", "Ketones" = "#ef6638", "Aldehydes" = "#8ed6d2", "Other alkenes, alkynes, dienes" = "#b569b3", "Other aromatics" = "#e7e85e", "Others" = "#2b9eb3", "Pentanes" = "#8c1531", "Propane" = "#9bb18d", "Propene" = "#623812", "Terpenes" = "#c9a415", "Toluene" = "#0352cb", "Trimethylbenzenes" = "#ae4901", "Xylenes" = "#1b695b", "CO" = "#6d6537", "Methane" = "#0c3f78", "Inorganic" = "#000000") `);
#$R->run(q` my.colours = c("Alkanes" = "#6c254f", "Alkenes,Alkynes" = "#f9c500", "Aromatics" = "#0e5c28", "Chlorinated" = "#ef6638", "CO" = "#6d6537", "Inorganic" = "#000000", "Methane" = "#0c3f78", "Oxygenated" = "#2b9eb3") `);
$R->run(q` plot.lines = function () { list( theme_tufte() ,
                                            ylab("Ox Production (molecules cm-3)") ,
                                            scale_y_continuous(expand = c(0, 0)) ,
                                            scale_x_discrete(expand = c(0, 0)) ,
                                            theme(axis.ticks.x = element_blank()) ,
                                            theme(legend.title = element_blank()) ,
                                            theme(axis.line = element_line(colour = "black")) ,
                                            theme(plot.title = element_text(face = "bold")) ,
                                            theme(axis.title.y = element_text(face = "bold")) ,
                                            theme(axis.title.x = element_blank()) ,
                                            theme(strip.text = element_text(face = "bold")) ,
                                            theme(strip.text.y = element_text(angle = 0)) ,
                                            theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 0.7, vjust = 0.8)) ,
                                            scale_fill_manual(values = my.colours, limits = rev(levels(data$Process))) ,
                                            theme(panel.margin = unit(5, "mm")) ) } `);

$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
#$R->run(q` data$Process = factor(data$Process, levels = c("Acids", "Alcohols", "Aldehydes", "Benzene", "Butanes", "CO", "Chlorinated", "Esters", "Ethane", "Ethene", "Ethers", "Higher alkanes", "Inorganic", "Ketones", "Methane", "Other alkenes, alkynes, dienes", "Other aromatics", "Pentanes", "Propane", "Propene", "Terpenes", "Toluene", "Trimethylbenzenes", "Xylenes")) `);
$R->run(q` data$Process = factor(data$Process, levels = c("Methane", "CO", "Inorganic", "Ethane", "Propane", "Butanes", "Pentanes", "Higher alkanes", "Ethene", "Propene", "Terpenes", "Other alkenes, alkynes, dienes", "Benzene", "Toluene", "Trimethylbenzenes", "Xylenes", "Other aromatics", "Acids", "Alcohols", "Aldehydes", "Esters", "Ethers", "Ketones", "Chlorinated")) `);
$R->run(q` plot.data = mutate(data, time = factor(Time, labels = c("First\nDay", "Cumulative"))) `);
$R->run(q` plot.data = mutate(plot.data, mechanism = factor(Mechanism, labels = c("MCM v3.2", "MOZART-4", "RADM2"))) `);
#my $p = $R->run(q` print(levels(data$Process)) `);
#print $p, "\n";

$R->run(q` plot = ggplot(plot.data, aes(x = Speciation, y = Ox.Production, fill = Process, order = Process)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_grid(time ~ mechanism) `,
        q` plot = plot + plot.lines() `,
);

$R->run(q` CairoPDF(file = "Absolute_Ox_budget_allocated_facet_mechanism.pdf", width = 8.6, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->run(q` plot = ggplot(plot.data, aes(x = mechanism, y = Ox.Production, fill = Process)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_grid(time ~ Speciation) `,
        q` plot = plot + plot.lines() `,
);

$R->run(q` CairoPDF(file = "Absolute_Ox_budget_allocated_facet_speciation.pdf", width = 8.6, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

###percent contributions
$R->run(q` get.fraction = function (d) {    d = d %>%   group_by(Speciation) %>%
                                                        mutate(First.Day.Total = sum(First.Day), Cumulative.Total = sum(Cumulative)) %>% 
                                                        group_by(Process) %>%
                                                        mutate(First.Day.Fraction = First.Day / First.Day.Total, Cumulative.Fraction = Cumulative / Cumulative.Total) %>%
                                                        select(-First.Day, -Cumulative, -First.Day.Total, -Cumulative.Total)
                                                        return(d) } `);
$R->run(q` wide = data %>% spread(Time, Ox.Production) `,
        q` radm2 = wide %>% filter(Mechanism == "RADM2") %>% arrange(Speciation) `, ##get mechanism data
        q` mozart = wide %>% filter(Mechanism == "MOZART") %>% arrange(Speciation) `, ##get mechanism data
        q` radm2 = get.fraction(radm2) `,
        q` mozart = get.fraction(mozart) `,
        q` fractions = rbind(mozart, radm2) `,
        q` fractions = gather(fractions, Time, Ox.Fraction, -Mechanism, -Speciation, -Process) `,
        q` plot.fractions = mutate(fractions, time = factor(Time, labels = c("First\nDay", "Cumulative"))) `,
        q` plot.fractions = mutate(plot.fractions, mechanism = factor(Mechanism, labels = c("MCM v3.2", "MOZART-4", "RADM2"))) `,
);
#my $p = $R->run(q` print(fractions) `);
#print $p, "\n";

$R->run(q` p = ggplot(plot.fractions, aes(x = Speciation, y = Ox.Fraction, fill = Process)) `,
        q` p = p + geom_bar(stat = "identity") `,
        q` p = p + facet_grid(time ~ mechanism) `,
        q` p = p + plot.lines() `,
        q` p = p + scale_y_continuous(expand = c(0, 0), labels = percent) `,
);

$R->run(q` CairoPDF(file = "Fractional_Ox_budget_allocated_facet_mechanism.pdf", width = 8.6, height = 6) `,
        q` print(p) `,
        q` dev.off() `,
);

$R->run(q` p = ggplot(plot.data, aes(x = mechanism, y = Ox.Production, fill = Process)) `,
        q` p = p + geom_bar(stat = "identity") `,
        q` p = p + facet_grid(time ~ Speciation) `,
        q` p = p + plot.lines() `,
        q` p = p + scale_y_continuous(expand = c(0, 0), labels = percent) `,
);

$R->run(q` CairoPDF(file = "Fractional_Ox_budget_allocated_facet_speciation.pdf", width = 8.6, height = 6) `,
        q` print(p) `,
        q` dev.off() `,
);

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

sub get_category {
    my ($process, $mechanism, $speciation) = @_;
    my $category;
    if ($process =~ /CO \+ OH/) {
        $category = "CO";
    } elsif ($process =~ /\+/) {
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
    } elsif ($process eq "HC8" or $process =~ /NC6|NC7|NC8|NC9|NC1\d|M\dPE|M\dHEX|CHEX/) {
        $category = "Higher alkanes";
    } elsif ($process =~ /CH3OH|C2H5OH|IPROPOL|NBUTOL|NPROPOL|BUT2OL|IBUTOL|MIBKAOH|C6H5CH2OH|ETHGLY|PROPGLY/) {
        $category = "Alcohols";
    } elsif ($process =~ /CH3OCH3|BUOX2ETOH|PR2OHMOX|EOX2EOL|MO2EOL/) {
        $category = "Ethers";
    } elsif ($process =~ /ACET$/) {
        $category = "Esters";
    } elsif ($process =~ /TRICLETH|CH3CCL3|CH2CL2|TCE/) {
        $category = "Chlorinated";
    } elsif ($process =~ /HCOOH|CH3COOH|PROPACID|CH3CO2H/) {
        $category = "Acids";
    } elsif ($process eq "KET" or $process =~ /MEK|CH3COCH3|MIBK|CYHEXONE/) {
        $category = "Ketones";
    } elsif ($process eq "ALD" or $process =~ /MACR|CH3CHO|C2H5CHO|C3H7CHO|IPRCHO|C4H9CHO|ACR|C4ALDB|HCHO|^CH2O/) {
        $category = "Aldehydes";
    } elsif ($process =~ /C5H8|ISO|BIGENE|BUT1ENE|C2H2/) {
        $category = "Other alkenes, alkynes, dienes";
    } elsif ($process eq "HC5" or $process eq "BIGALK" or $process =~ /NC5|IC5|NEOP/) {
        $category = "Pentanes";
    } elsif ($process eq "HC3" or $process eq "C3H8") {
        $category = "Propane";
    } elsif ($process eq "OLT" or $process eq "C3H6") {
        $category = "Propene";
    } elsif ($process =~ /C10H16|LIMONENE|PINENE|OLI/) {
        $category = "Terpenes";
    } elsif ($process =~ /Inorganic|Emissions/) {
        $category = $process;
    } else {
        $category = "Inorganic";
        #print "No category found for $process in $mechanism, $speciation\n";
    }
    return $category;
}

sub get_data {
    my ($mecca, $kpp, $mechanism, $speciation, $run, $no_runs) = @_;
    $no_runs = 1 unless (defined $no_runs);
    my (%production_rates, %consumption_rates, $producers, $producer_yields, $consumers, $consumer_yields);
    foreach my $species (qw( Ox HO2x )) {
        $kpp->family({
                name    => $species,
                members => $families{$species},
                weights => $weights{$species},
        });
        $producers = $kpp->producing($species);
        $producer_yields = $kpp->effect_on($species, $producers);
        $consumers = $kpp->consuming($species);
        $consumer_yields = $kpp->effect_on($species, $consumers);

        for (0..$#$producers) {
            my $reaction = $producers->[$_];
            my $reaction_number = $kpp->reaction_number($reaction);
            my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
            next if ($rate->sum == 0);
            my ($number, $parent) = split /_/, $reaction;
            if (defined $parent) {
                $production_rates{$species}{$parent} += $rate(1:$ntime-2);
            } else {
                my $reaction_string = $kpp->reaction_string($reaction);
                if ($reaction_string =~ /_notag/) {
                    $production_rates{$species}{"notag"} += $rate(1:$ntime-2);
                } else {
                    my ($reactants, $products) = split / = /, $reaction_string;
                    $production_rates{$species}{$reactants} += $rate(1:$ntime-2);
                }
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
                my $reaction_string = $kpp->reaction_string($reaction);
                if ($reaction_string =~ /_notag/) {
                    $consumption_rates{$species}{"notag"} += $rate(1:$ntime-2);
                } else {
                    my ($reactants, $products) = split / = /, $reaction_string;
                    $consumption_rates{$species}{$reactants} += $rate(1:$ntime-2);
                }
            }
        }
    } 
    remove_common_processes($production_rates{"HO2x"}, $consumption_rates{"HO2x"});
    my $total_ho2x_production;
    $total_ho2x_production += $production_rates{"HO2x"}{$_} foreach (keys %{$production_rates{"HO2x"}}); 
    foreach my $reaction (keys %{$production_rates{'HO2x'}}) {
        $production_rates{"Ox"}{$reaction} += $production_rates{"Ox"}{'HO2 + NO'} * $production_rates{'HO2x'}{$reaction} / $total_ho2x_production;
    }
    delete $production_rates{"Ox"}{'HO2 + NO'};
    delete $production_rates{"Ox"}{'notag'};
    remove_common_processes($production_rates{"Ox"}, $consumption_rates{"Ox"});

    my %final_categories;
    foreach my $process (keys %{$production_rates{'Ox'}}) {
        if (exists $category_mapping{$mechanism}) {
            if (exists $category_mapping{$mechanism}{$speciation}) {
                if (exists $category_mapping{$mechanism}{$speciation}{$process}) {
                    foreach my $entry (@{$category_mapping{$mechanism}{$speciation}{$process}}) {
                        my ($factor, $category) = split / /, $entry;
                        $category =~ s/_/ /g;
                        $final_categories{$category} += $factor * $production_rates{"Ox"}{$process};
                    }
                } else {
                    my $category = get_category($process, $mechanism, $speciation);
                    $final_categories{$category} += $production_rates{'Ox'}{$process};
                }
            }
        } else {
            my $category = get_category($process, $mechanism, $speciation);
            $final_categories{$category} += $production_rates{"Ox"}{$process};
        }
    }
    foreach my $category (keys %final_categories) {
        my $reshape = $final_categories{$category}->reshape($n_per_day, $n_days);
        $final_categories{$category} = $reshape->sumover;
        if ($category eq "Methane" or $category eq "Inorganic" or $category eq "CO") {
            #print "$category and $no_dirs\n";
            $final_categories{$category} /= $no_runs;
        } 
    }
#    my %last;
#    foreach my $category (keys %final_categories) {
#        if ($category =~ /Chlo/) {
#            $last{"Chlorinated"} += $final_categories{$category};
#        } elsif ($category =~ /Toluene|Benzene|Other\sArom|Xylen|Trime/) {
#            $last{"Aromatics"} += $final_categories{$category};
#        } elsif ($category =~ /Acid|Alcohol|Aldeh|Keton|Este|Ether/) {
#            $last{"Oxygenated"} += $final_categories{$category};
#        } elsif ($category =~ /Ethen|Other\sAlk|Prope|Terpe/) {
#            $last{"Alkenes,Alkynes"} += $final_categories{$category};
#        } elsif ($category =~ /Etha|Propa|Buta|Pent|High/) {
#            $last{"Alkanes"} += $final_categories{$category};
#            #} elsif ($category =~ /Meth|Inor|CO/) {
#            #$last{$category} = $final_categories{$category}
#        }
#    }
    #print "$run\n";
    #print "$_ : $final_categories{$_}\n" foreach keys %final_categories;
    #return \%last;
    return \%final_categories;
}
