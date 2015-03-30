#! /usr/bin/env perl
# Compare allocated production budgets of HCHO in tagged and non-tagged runs.
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
my @VOC = qw(ACR ALD APINENE BENZENE BIGALK BIGENE BPINENE BUOX2ETOH BUT1ENE BUT2OL C10H16 C2H2 C2H4 C2H5CHO C2H5OH C2H6 C3H6 C3H7CHO C3H8 C4ALDB C4H9CHO C5H8 C6H5CH2OH CH2CL2 CH2O CH3CCL3 CH3CHO CH3CO2H CH3COCH3 CH3COOH CH3OCH3 CH3OH CHEX CYHEXONE DIME35EB EBENZ EOX2EOL ETH ETHACET ETHGLY HC3 HC5 HC8 HCHO HCOOH IBUTOL IC4H10 IC5H12 IPBENZ IPRCHO IPROACET IPROPOL ISO ISOP KET LIMONENE M2HEX M2PE M3HEX M3PE MACR MEK METHTOL MIBK MIBKAOH MO2EOL MXYL NBUTACET NBUTOL NC10H22 NC11H24 NC12H26 NC4H10 NC5H12 NC6H14 NC7H16 NC8H18 NC9H20 NEOP NPROACET NPROPOL OL2 OLI OLT ORA1 ORA2 OXYL PBENZ PETHTOL PR2OHMOX PROPACID PROPGLY PXYL STYRENE TCE TM123B TM124B TM135B TOL TOLUENE TRICLETH XYL );

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
                    $data{$mechanism}{$speciation}{$run}{$dir} = get_data($mecca, $kpp, "HCHO", $mechanism, $speciation, $run, $dir, $no_dirs);
                }
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
foreach my $speciation (sort keys %{$data{"MCM"}}) {
    my %allocated;
    foreach my $run (sort keys %{$data{"MCM"}{$speciation}}) {
        next unless ($run =~ /tagged/);
        foreach my $dir (sort keys %{$data{'MCM'}{$speciation}{$run}}) {
            foreach my $item (sort keys %{$data{'MCM'}{$speciation}{$run}{$dir}}) {
                $allocated{$item} += $data{'MCM'}{$speciation}{$run}{$dir}{$item};
            }
        }
        $data{'MCM'}{$speciation}{$run} = \%allocated;
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
    #print "$mechanism\n";
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        $R->set('speciation', $speciation);
        #print "\t$speciation\n";
        foreach my $run_type (sort keys %{$data{$mechanism}{$speciation}}) {
            #print "\t\t$run_type : $data{$mechanism}{$speciation}{$run_type}\n";
            $R->run(q` pre = data.frame(Mechanism = mechanism) `);
            $R->set('run', $run_type);
            $R->run(q` pre$Run = run `);
            foreach my $process (sort keys %{$data{$mechanism}{$speciation}{$run_type}}) {
                #print "$process: $data{$mechanism}{$speciation}{$run_type}{$process}\n";
                $R->set('process', $process);
                $R->set('total.prod', $data{$mechanism}{$speciation}{$run_type}{$process});
                $R->run(q` pre[process] = total.prod `);
            }
            $R->run(q` pre$Speciation = speciation `);
            $R->run(q` pre = gather(pre, Process, Total.Prod, -Run, -Mechanism, -Speciation) `,
                    q` data = rbind(data, pre) `,
            );
        }
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
#$R->run(q` my.colours = c("Non-Tagged" = "#2b9eb3", "Tagged" = "#ef6638") `);

$R->run(q` plot = ggplot(data, aes(x = Run, y = Total.Prod, fill = Process)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_grid(Mechanism ~ Speciation ) `,
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
        #q` plot = plot + scale_fill_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "HCHO_allocated_production_tagged_vs_non_tagged.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $species, $mechanism, $speciation, $run, $dir, $no_dirs) = @_;
    my $ntime = $mecca->time->nelem;

    my (%production_rates, %consumption_rates, $producers, $producer_yields, $consumers, $consumer_yields);
    if ($run =~ /tagged/) {
        my $spc;
        if ($mechanism eq "MCM") {
            $spc = "$dir/gas.spc";
        } else {
            $spc = "$base/$mechanism/${speciation}_$run/gas.spc";
            $no_dirs = 1;
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
        $no_dirs = 1;
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
        my ($number, $parent) = split /_/, $reaction;
        if (defined $parent) {
            $production_rates{$parent} += $rate(1:$ntime-2);
        } else {
            my $reactants = $kpp->reactants($reaction);
            my $VOC = undef;
            foreach (@$reactants) {
                $VOC = $_ if ($_ ~~ @VOC);
            }
            if (defined $VOC) {
                $production_rates{$VOC} += $rate(1:$ntime-2);
            } else {
                my ($reactants, $products) = split / = /, $kpp->reaction_string($reaction);
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
            $consumption_rates{$reactants} += $rate(1:$ntime-2);
        }
    }
    #remove_common_processes(\%production_rates, \%consumption_rates);
    my $others = 4e7;
    foreach my $process (keys %production_rates) {
        #print "$process: ", $production_rates{$process}->sum, "\n";
        if ($production_rates{$process}->sum < $others) {
            $production_rates{"Others"} += $production_rates{$process}->sum;
            delete $production_rates{$process};
        } else {
            $production_rates{$process} = $production_rates{$process}->sum;
        }
    }
    $production_rates{$_} = $production_rates{$_} / $no_dirs foreach (keys %production_rates);
    return \%production_rates;
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
