#! /usr/bin/env perl
# Compare HOx production budgets from the higher alkanes between speciations and mechanisms, Solvent only runs, allocated to reactions
# Version 0: Jane Coates 4/5/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

#my $base = "/work/users/jco/Solvent_Emissions";
my $base = "/local/home/coates/Solvent_Emissions";
#my @mechanisms = qw( MOZART );
my @mechanisms = qw( MCM MOZART RADM2 );
#my @speciations = qw( TNO );
my @speciations = qw( TNO IPCC DE94 GR95 GR05 UK98 UK08 );
my (%families, %weights, %data);
$families{"HOx"} = [ qw( OH HO2 HO2NO2 ) ];

my $others_max = 6e7;
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
    foreach my $reactants (keys %allocated) {
        if ($allocated{$reactants} < $others_max) {
            $allocated{"Others"} += $allocated{$reactants};
            delete $allocated{$reactants};
        }
    }
    my @sorted = reverse sort { $allocated{$a} <=> $allocated{$b} } keys %allocated;
    my @final_sorted;
    foreach (@sorted) {
        next if ($_ =~ /Others/);
        push @final_sorted, { $_ => $allocated{$_} };
    }
    push @final_sorted, { "Others" => $allocated{"Others"} } if (defined $allocated{"Others"});
    $data{"MCM"}{$speciation} = \@final_sorted;
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
        foreach my $ref (@{$data{$mechanism}{$speciation}}) {
            foreach my $reactants (sort keys %$ref) {
                $R->set('reactants', $reactants);
                $R->set('hox', $ref->{$reactants});
                $R->run(q` pre[reactants] = hox `);
            }
        }
        $R->run(q` pre = gather(pre, Reactants, HOx, -Mechanism, -Speciation) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
#$R->run(q` data$Reactants = factor(data$Reactants, levels = c("NC6H14", "M2PE", "M3PE", "CHEX", "NC7H16", "M2HEX", "M3HEX", "NC8H18", "NC9H20", "NC10H22", "NC11H24", "NC12H26")) `);
#$R->run(q` my.colours = c("NC6H14" = "#6c254f", "M2PE" = "#f9c500", "M3PE" = "#0e5c28", "CHEX" = "#ef6638", "NC7H16" = "#2b9eb3", "M2HEX" = "#b569b3", "M3HEX" = "#f7c56c", "NC8H18" = "#0352cb", "NC9H20" = "#ae4901", "NC10H22" = "#4c9383", "NC11H24" = "#8c1531", "NC12H26" = "#77aecc") `); 
#$R->run(q` write.table(data, file = "HOx_higher_alkanes_in_MCM_species.csv", sep = ",", row.name = FALSE, quote = FALSE) `);

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = HOx, fill = Reactants)) `,
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
        #q` plot = plot + scale_fill_manual(values = my.colours, limits = rev(levels(data$Reactants))) `,
);

$R->run(q` CairoPDF(file = "Higher_alkanes_HOx_production_reactions.pdf", width = 10, height = 7) `,
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
        my ($reactants, $products) = split / = /, $kpp->reaction_string($reaction);
        $reactants =~ s/_(.*?)\b//g;
        $production_rates{$parent}{$reactants} += $rate(1:$ntime-2);
    }

    for (0..$#$consumers) {
        my $reaction = $consumers->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number) * $consumer_yields->[$_];
        next if ($rate->sum == 0); 
        my ($number, $parent) = split /_/, $reaction;
        next unless ($parent ~~ @higher_alkanes);
        my ($reactants, $products) = split / = /, $kpp->reaction_string($reaction);
        $reactants =~ s/_(.*?)\b//g;
        $consumption_rates{$parent}{$reactants} += $rate(1:$ntime-2);
    }

    print "before\n";
    foreach my $parent (keys %production_rates) {
        print "$_ : ", $production_rates{$parent}{$_}->sum, "\n" foreach keys %{$production_rates{$parent}};
        remove_common_processes($production_rates{$parent}, $consumption_rates{$parent});
    }

    my %contributions;
    if ($mechanism eq "RADM2" or $mechanism eq "MOZART") { #allocation of RADM2 and MOZART emitted species to Higher Alkanes
        open my $in, '<:encoding(utf-8)', "Higher_alkanes_RADM2_MOZ_fractional_contributions.csv" or die $!;
        my @lines = <$in>;
        close $in;
        foreach my $line (@lines) {
            chomp $line;
            my ($mechanism, $speciation, $species, $contribution) = split /,/, $line;
            $contributions{$mechanism}{$speciation}{$species} = $contribution;
        }

        print "after\n";
        foreach my $parent (keys %production_rates) {
            foreach my $reactants (keys %{$production_rates{$parent}}) {
                print "$reactants : ", $production_rates{$parent}{$reactants}->sum, "\n";
                if (defined $contributions{$mechanism}{$speciation}{$parent}) {
                    $production_rates{$parent}{$reactants} *= $contributions{$mechanism}{$speciation}{$parent};
                } else {
                    delete $production_rates{$parent}{$reactants};
                } 
            }
        }
    }
    my %final;
    foreach my $parent (keys %production_rates) {
        foreach my $reactants (keys %{$production_rates{$parent}}) {
            #print "$reactants : ", $production_rates{$parent}{$reactants}->sum / $no_dirs, "\n";
            $final{$reactants} += $production_rates{$parent}{$reactants}->sum / $no_dirs;
        }
    }

    if ($mechanism =~ /MOZ|RAD/) {
        foreach my $reaction (keys %final) {
            if ($final{$reaction} < $others_max) {
                $final{"Others"} += $final{$reaction};
                delete $final{$reaction};
            }
        }
        my @sorted = reverse sort { $final{$a} <=> $final{$b} } keys %final;
        my @final_sorted;
        foreach (@sorted) {
            next if ($_ =~ /Others/);
            push @final_sorted, { $_ => $final{$_} };
        }
        push @final_sorted, { "Others" => $final{"Others"} } if (defined $final{"Others"});
        return \@final_sorted;
    } else {
        return \%final;
    }
}
