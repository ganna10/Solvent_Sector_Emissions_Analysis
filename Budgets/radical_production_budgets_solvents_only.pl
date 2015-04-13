#! /usr/bin/env perl
# Compare radical production budgets allocated to reactions in solvents only runs, facet by mechanism and speciation
# Version 0: Jane Coates 10/4/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/work/users/jco/Solvent_Emissions";
#my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
#my @speciations = qw( TNO IPCC );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my (%families, %weights, %data);

foreach my $mechanism (@mechanisms) { 
    foreach my $speciation (@speciations) {
        my $boxmodel = "$base/$mechanism/${speciation}_Solvents_Only/boxmodel";
        my $mecca = MECCA->new($boxmodel);
        my $eqn  = "$base/$mechanism/${speciation}_Solvents_Only/gas.eqn";
        my $kpp = KPP->new($eqn);
        my $radical_file = "$base/$mechanism/${speciation}_Solvents_Only/radicals.txt";
        my @radicals = get_radicals($radical_file);
        $families{"radicals"} = [ @radicals ];
        $data{$mechanism}{$speciation} = get_data($mecca, $kpp, $mechanism, $speciation);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        $R->run(q` pre = data.frame(Mechanism = mechanism) `);
        $R->set('speciation', $speciation);
        $R->run(q` pre$Speciation = speciation `); 
        foreach my $ref (@{$data{$mechanism}{$speciation}}) {
            foreach my $reaction (sort keys %$ref) {
                $R->set('reaction', $reaction);
                $R->set('radical.prod', $ref->{$reaction});
                $R->run(q` pre[reaction] = radical.prod `);
            }
        }
        $R->run(q` pre = gather(pre, Reaction, Radical.Prod, -Mechanism, -Speciation) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);

$R->run(q` plot = ggplot(data, aes(x = Speciation, y = Radical.Prod, fill = Reaction)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism, nrow = 1 ) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
);

$R->run(q` CairoPDF(file = "radical_production_facet_mechanism.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = Radical.Prod, fill = Reaction)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2 ) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
);

$R->run(q` CairoPDF(file = "radical_production_facet_speciation.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $mechanism, $speciation) = @_;
    my $ntime = $mecca->time->nelem;
    my (%production_rates, %consumption_rates);
    $kpp->family({
            name    => "radicals",
            members => $families{"radicals"},
            weights => $weights{"radicals"},
    });
    my $producers = $kpp->producing("radicals");
    my $producer_yields = $kpp->effect_on("radicals", $producers);
    my $consumers = $kpp->consuming("radicals");
    my $consumer_yields = $kpp->effect_on("radicals", $consumers);

    for (0..$#$producers) {
        my $reaction = $producers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $producer_yields->[$_];
        my ($reactants, $products) = split / = /, $kpp->reaction_string($reaction);
        $production_rates{$reactants} += $rate(1:$ntime-2);
    } 

    for (0..$#$consumers) {
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
        my ($reactants, $products) = split / = /, $kpp->reaction_string($reaction);
        $consumption_rates{$reactants} += $rate(1:$ntime-2);
    } 
    remove_common_processes(\%production_rates, \%consumption_rates);
    my $others = 8e7;
    foreach my $reaction (keys %production_rates) {
        if ($production_rates{$reaction}->sum < $others) {
            $production_rates{"Others"} += $production_rates{$reaction}->sum;
            delete $production_rates{$reaction};
        } else {
            $production_rates{$reaction} = $production_rates{$reaction}->sum;
        }
    }
    my @sorted_prod = reverse sort { $production_rates{$a} <=> $production_rates{$b} } keys %production_rates;
    my @final_sorted;
    foreach (@sorted_prod) {
        next if ($_ =~ /Others/);
        push @final_sorted, { $_ => $production_rates{$_} };
    }
    push @final_sorted, { "Others" => $production_rates{"Others"} } if (defined $production_rates{"Others"});
    return \@final_sorted;
}

sub get_radicals {
    my ($file) = @_;
    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file : $!";
    local $/ = undef; 
    my @radicals = split /\s+/, <$in>;
    close $in;
    return @radicals;
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
