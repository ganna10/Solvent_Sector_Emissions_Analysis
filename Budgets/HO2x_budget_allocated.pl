#! /usr/bin/env perl
# Allocate Cumulative HO2x production budgets back to categories used to defined speciations
# Jane Coates 4/4/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions";
#my $base = "/work/users/jco/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
#my @speciations = qw( UK08 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my (%families, %weights, %data);
$families{"HO2x"} = [ qw( HO2 HO2NO2 ) ];

my $mecca = MECCA->new("$base/RADM2/TNO_Solvents_Only/boxmodel");
my $ntime = $mecca->time->nelem;

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


foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my (%production_rates, %consumption_rates);
        if ($mechanism =~ /MOZART|RADM2/) { 
            my $boxmodel = "$base/$mechanism/${speciation}_tagged_solvents_only/boxmodel";
            my $mecca = MECCA->new($boxmodel);
            my $eqn_file = "$base/$mechanism/${speciation}_tagged_solvents_only/gas.eqn";
            my $kpp = KPP->new($eqn_file); 
            $data{$mechanism}{$speciation} = get_data($mecca, $kpp, $mechanism, $speciation);
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
                $data{$mechanism}{$speciation}{$run} = get_data($mecca, $kpp, $mechanism, $speciation, $run, $run_number);
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
        q` library(Cairo) `,
        q` library(ggthemes) `,
        q` library(grid) `,
);
$R->run(q` my.colours = c("Acids" = "#cc6329", "Alcohols" = "#6c254f", "Benzene" = "#8c6238", "Butanes" = "#86b650", "Chlorinated" = "#f9c500", "CO" = "#898989", "Esters" = "#f3aa7f", "Ethane" = "#77aecc", "Ethene" = "#1c3e3d", "Ethers" = "#ba8b01", "Higher alkanes" = "#0e5c28", "Ketones" = "#ef6638", "Aldehydes" = "#8ed6d2", "Other alkenes, alkynes, dienes" = "#b569b3", "Other aromatics" = "#e7e85e", "Others" = "#2b9eb3", "Pentanes" = "#8c1531", "Propane" = "#9bb18d", "Propene" = "#623812", "Terpenes" = "#c9a415", "Toluene" = "#0352cb", "Trimethylbenzenes" = "#ae4901", "Xylenes" = "#1b695b", "CO" = "#6d6537", "Methane" = "#0c3f78", "Inorganic" = "#000000") `);

$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    #print "$mechanism\n";
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        #print "\t$speciation\n";
        $R->run(q` pre = data.frame(Mechanism = mechanism) `);
        $R->set('speciation', $speciation);
        $R->run(q` pre$Speciation = speciation `);
        foreach my $type (sort keys %{$data{$mechanism}{$speciation}}) {
            #print "\t\t$type\n";
            $R->set('type', $type);
            $R->set('HO2x.production', $data{$mechanism}{$speciation}{$type});
            $R->run(q` pre[type] = HO2x.production `);
        }
        $R->run(q` pre = gather(pre, Type, HO2x.Production, -Mechanism, -Speciation) `,
                q` data = rbind(data, pre) `,
        );
        #my $p = $R->run(q` print(pre) `);
        #print $p, "\n";
    }
}
$R->run(q` data$Type = factor(data$Type, levels = c("Acids", "Alcohols", "Aldehydes", "Benzene", "Butanes", "CO", "Chlorinated", "Esters", "Ethane", "Ethene", "Ethers", "Higher alkanes", "Inorganic", "Ketones", "Methane", "Other alkenes, alkynes, dienes", "Other aromatics", "Pentanes", "Propane", "Propene", "Terpenes", "Toluene", "Trimethylbenzenes", "Xylenes")) `);
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` plot.lines = function () { list( theme_tufte() ,
                                            ggtitle("Cumulative HO2x Production Budget") ,
                                            ylab("HO2x Production (molecules cm-3)") ,
                                            scale_y_continuous(expand = c(0, 0)) ,
                                            scale_x_discrete(expand = c(0, 0)) ,
                                            theme(axis.ticks.x = element_blank()) ,
                                            theme(legend.title = element_blank()) ,
                                            theme(axis.line = element_line(colour = "black")) ,
                                            theme(plot.title = element_text(face = "bold")) ,
                                            theme(axis.title.y = element_text(face = "bold")) ,
                                            theme(axis.title.x = element_blank()) ,
                                            theme(strip.text = element_text(face = "bold")) ,
                                            theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 0.7, vjust = 0.8)) ,
                                            scale_fill_manual(values = my.colours, limits = rev(levels(data$Type))) ,
                                            theme(panel.margin = unit(5, "mm")) ) } `);

$R->run(q` plot = ggplot(data, aes(x = Speciation, y = HO2x.Production, fill = Type)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism, nrow = 1 ) `,
        q` plot = plot + plot.lines() `,
);

$R->run(q` CairoPDF(file = "Cumulative_HO2x_budget_allocated_facet_mechanism.pdf", width = 8.6, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = HO2x.Production, fill = Type)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2 ) `,
        q` plot = plot + plot.lines() `,
);

$R->run(q` CairoPDF(file = "Cumulative_HO2x_budget_allocated_facet_speciation.pdf", width = 8.6, height = 6) `,
        q` print(plot) `,
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

sub get_data {
    my ($mecca, $kpp, $mechanism, $speciation, $dir, $no_runs) = @_;
    my (%production_rates, %consumption_rates);
    $kpp->family({
            name    => "HO2x",
            members => $families{"HO2x"},
            weights => $weights{"HO2x"},
    });
    my $producers = $kpp->producing("HO2x");
    my $producer_yields = $kpp->effect_on("HO2x", $producers);
    my $consumers = $kpp->consuming("HO2x");
    my $consumer_yields = $kpp->effect_on("HO2x", $consumers);
    print "No consumers in $mechanism and $speciation\n" if (@$consumers == 0);
    print "No producers in $mechanism and $speciation\n" if (@$producers == 0);

    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        next if ($rate->sum == 0);
        my ($number, $parent) = split /_/, $reaction;
        if (defined $parent) {
            $production_rates{$parent} += $rate(1:$ntime-2);
        } else {
            my ($reactants, $products) = split / = /, $kpp->reaction_string($reaction);
            if ($products =~ /notag/) {
                $production_rates{"notag"} += $rate(1:$ntime-2);
            } else {
                $production_rates{$reactants} += $rate(1:$ntime-2);
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
            $consumption_rates{$parent} += $rate(1:$ntime-2);
        } else {
            my ($reactants, $products) = split / = /, $kpp->reaction_string($reaction);
            if ($products =~ /notag/) {
                $consumption_rates{"notag"} += $rate(1:$ntime-2);
            } else {
                $consumption_rates{$reactants} += $rate(1:$ntime-2);
            }
        }
    }
    remove_common_processes(\%production_rates, \%consumption_rates);
    delete $production_rates{"notag"};
    my %final_categories;
    foreach my $process (keys %production_rates) {
        if (exists $category_mapping{$mechanism}) {
            if (exists $category_mapping{$mechanism}{$speciation}) {
                if (exists $category_mapping{$mechanism}{$speciation}{$process}) {
                    foreach my $entry (@{$category_mapping{$mechanism}{$speciation}{$process}}) {
                        my ($factor, $category) = split / /, $entry;
                        $category =~ s/_/ /g;
                        $final_categories{$category} += $factor * $production_rates{$process}->sum;
                    }
                } else {
                    my $category = get_category($process, $mechanism, $speciation);
                    $final_categories{$category} += $production_rates{$process}->sum;
                }
            }
        } else {
            my $category = get_category($process, $mechanism, $speciation);
            $final_categories{$category} += $production_rates{$process}->sum;
        }
    }
    if ($mechanism eq "MCM") {
        foreach my $category (keys %final_categories) {
            if ($category eq "Methane" or $category eq "CO" or $category eq "Inorganic") {
                $final_categories{$category} /= $no_runs;
            }
        }
    }
    return \%final_categories;
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
    } elsif ($process =~ /Deposition|Inorganic|Emissions/) {
        $category = $process;
    } else {
        print "No category found for $process in $mechanism, $speciation\n";
    }
    return $category;
}
