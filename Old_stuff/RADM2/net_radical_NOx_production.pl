#! /usr/bin/env perl
# Plot net radical to NOx production for each speciation seeing which reactions contribute to NO source calculation
# Version 0: Jane Coates 17/11/2014

use strict; 
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions/RADM2";
my $mecca = MECCA->new("$base/DE94_Solvents_Only/boxmodel");
my $NTIME = $mecca->time->nelem;
my $DT = $mecca->dt->at(0);
my $N_PER_DAY = 43200 / $DT;
my $N_DAYS = int $NTIME / $N_PER_DAY;

opendir DIR, $base or die "Can't open $base\n";
my @runs = grep { $_ =~ /_Solvents_Only/ } readdir DIR;
closedir DIR;

my (%families, %weights, %plot_data);

foreach my $run (@runs) {
    (my $speciation = $run) =~ s/^(.*?)_Solvents_Only/$1/;
    my $boxmodel = "$base/$run/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn_file = "$base/$run/gas.eqn";
    my $kpp = KPP->new($eqn_file);
    my $radicals_file = "$base/$run/radicals.txt"; 
    my @radicals = get_species($radicals_file);
    $families{$speciation} = [ @radicals ];
    $plot_data{$speciation} = get_data($mecca, $kpp, $speciation);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(reshape2) `,
        q` library(Cairo) `,
        q` library(grid) `,
);

$R->set('Time', [ ("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7") ]);
$R->run(q` data = data.frame() `);

foreach my $speciation (sort keys %plot_data) {
    $R->run(q` pre = data.frame(Time) `);
    foreach my $ref (@{$plot_data{$speciation}}) {
        foreach my $reaction (sort keys %$ref) {
            $R->set('reaction', $reaction);
            $R->set('rate', [ map { $_ } $ref->{$reaction}->dog ]);
            $R->run(q` pre[reaction] = rate `);
        }
    }
    $R->set('speciation', $speciation);
    $R->run(q` pre$Speciation = rep(speciation, length(Time)) `,
            q` pre = melt(pre, id.vars = c("Time", "Speciation"), variable.name = "Reaction", value.name = "Rate") `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` my.colours = c(  "O1D" = "#6c254f",
                            "HCHO + hv" = "#f9c500",
                            "MGLY + hv" = "#0352cb",
                            "KET + hv" = "#ef6638",
                            "HC8 + OH" = "#0e5c28",
                            "TCO3 + NO" = "#77aecc",
                            "DCB + hv" = "#b569b3",
                            "Production Others" = "#8c6238" ) `);
$R->run(q` data$Reaction = factor(data$Reaction, levels = c("O1D", "HCHO + hv", "MGLY + hv", "KET + hv", "HC8 + OH", "TCO3 + NO", "DCB + hv", "Production Others")) `);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Rate, fill = Reaction)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Speciation ) `,
        q` plot = plot + theme_bw() `,
        q` plot = plot + theme(panel.grid = element_blank()) `,
        q` plot = plot + theme(strip.background = element_blank()) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + scale_fill_manual(values = my.colours, limits = rev(levels(data$Reaction))) `,
);

$R->run(q` CairoPDF(file = "net_radical_NOx_production.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp, $speciation) = @_;

    $families{'NOx'} = [ qw( NO NO2 NO3 N2O5 ) ];
    $weights{'NOx'} = { N2O5 => 2 };
    $kpp->family({ #NOx family
                    name    => 'NOx',
                    members =>$families{'NOx'},
                    weights => $weights{'NOx'},
    });

    my ($radical_producers, $radical_producer_yields, $NOx_yields, %production_rates);
    if (exists $families{$speciation}) { #radicals family
        $kpp->family({                                                                                                                           
                        name    => $speciation,
                        members => $families{$speciation},
                        weights => $weights{$speciation},
        }); 
        $radical_producers = $kpp->producing($speciation);
        $radical_producer_yields = $kpp->effect_on($speciation, $radical_producers);
        $NOx_yields = $kpp->effect_on('NOx', $radical_producers);
    } else {
        print "No radical family found for $speciation\n";
    }
    
    die "No producers found for $speciation\n" if (@$radical_producers == 0);

    for (0..$#$radical_producers) { #get rates for all radical producing reactions
        next if ($NOx_yields->[$_] < 0);
        my $reaction = $radical_producers->[$_];
        my $reaction_string = $kpp->reaction_string($reaction);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $net_radical_yield = $radical_producer_yields->[$_] - $NOx_yields->[$_];
        next if ($net_radical_yield == 0);
        my $rate = $net_radical_yield * $mecca->rate($reaction_number); 
        next if ($rate->sum == 0); # do not include reactions that do not occur 
        $reaction_string =~ s/_(.*?)\b//g;
        my ($reactants, $products) = split / = /, $reaction_string;
        $production_rates{$reactants} += $rate(1:$NTIME-2);
    }

    my $others = 3.0e7;
    foreach my $reaction (keys %production_rates) {
        if ($production_rates{$reaction}->sum < $others) {
            $production_rates{'Production Others'} += $production_rates{$reaction};
            delete $production_rates{$reaction};
        }
    }

    foreach my $reaction (keys %production_rates) {
        my $reshape = $production_rates{$reaction}->copy->reshape($N_PER_DAY, $N_DAYS);
        my $integrate = $reshape->sumover;
        $integrate = $integrate(0:13:2);
        $production_rates{$reaction} = $integrate;
    }

    my $sort_function = sub { $_[0]->sum };
    my @sorted_plot_data;
    my @sorted_prod = sort { &$sort_function($production_rates{$b}) <=> &$sort_function($production_rates{$a}) } keys %production_rates;
    foreach (@sorted_prod) { #sum up rates of reactions, starting with reaction with lowest sum, production others added separately 
        next if ($_ eq 'Production Others');
        push @sorted_plot_data, { $_ => $production_rates{$_} };
    }

    push @sorted_plot_data, { 'Production Others' => $production_rates{'Production Others'} } if (defined $production_rates{'Production Others'}); #add Production Others to the beginning 
    return \@sorted_plot_data; 
}

sub get_species {
    my ($file) = @_;
    open my $in, '<:encoding(utf-8)', $file or die $!;
    local $/ = undef;
    my $lines = <$in>;
    close $in;
    my @species = split /\s+/, $lines;
}
