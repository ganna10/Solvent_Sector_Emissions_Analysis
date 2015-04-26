#! /usr/bin/env perl
# look at emission fraction of VOC type multplied by cumulative TOPPs of that type
# Version 0: Jane Coates 26/4/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw( MCM );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my (%emissions, %TOPPs);

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $dir = "$base/$mechanism/${speciation}_Solvents_Only";
        my $mecca = MECCA->new("$dir/boxmodel");
        my $kpp = KPP->new("$dir/gas.eqn");
        my $TOPP_file = "${mechanism}_${speciation}_${speciation}_tagged_solvents_only_cumulative_TOPP_values.csv";
        $TOPPs{$mechanism}{$speciation} = get_TOPPs($TOPP_file);
        $emissions{$mechanism}{$speciation} = get_emission_fractions($kpp, $mecca);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(scales) `,
        q` library(grid) `,
        q` library(ggthemes) `,
);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %emissions) {
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$emissions{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->run(q` pre = data.frame(Mechanism = mechanism, Speciation = speciation) `);
        $R->set('total.emissions', $emissions{$mechanism}{$speciation}{"Total"});
        foreach my $type (sort keys %{$emissions{$mechanism}{$speciation}}) {
            next if ($type eq "Total");
            $R->set('type', $type);
            $R->set('emission.fraction', $emissions{$mechanism}{$speciation}{$type});
            $R->set('topp', $TOPPs{$mechanism}{$speciation}{$type});
            $R->run(q` pre$Type = type `,
                    q` pre$Emission.Fraction = emission.fraction `,
                    q` pre$TOPP = topp `,
            );
            $R->run(q` pre = mutate(pre, OF = TOPP * Emission.Fraction * total.emissions) `,
                    q` data = rbind(data, pre) `,
            );
        }
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` plot = ggplot(data, aes(x = Speciation, y = OF, fill = Type)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
);
$R->run(q` CairoPDF(file = "Ox_Formation_by_mechanism.pdf", width = 8.6, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_TOPPs {
    my ($file) = @_;
    my %topps;
    open my $in, '<:encoding(utf-8)', $file or die $!;
    my @lines = <$in>;
    close $in;
    foreach my $line (@lines) {
        next if ($line =~ /^VOC/);
        chomp $line;
        my ($VOC, $value) = split /,/, $line; 
        if ($VOC =~ /^NC|^IC|^NEOP|^M2|^M3/ or $VOC eq "C3H8" or $VOC eq "CHEX" or $VOC eq "C2H6") {
            #print "$VOC : alkane\n";
            $topps{"Alkanes"} += $value;
        } elsif ($VOC eq "C2H4" or $VOC eq "C3H6" or $VOC eq "C2H2" or $VOC eq "LIMONENE") {
            #print "$VOC : alkene or alkyne\n";
            $topps{"Alkenes, Alkynes"} += $value;
        } elsif ($VOC =~ /^TOL|XYL|^TM|BENZ|^STY|ETHTOL$|^DIM/) {
            #print "$VOC : aromatic\n";
            $topps{"Aromatics"} += $value;
        } elsif ($VOC eq "CH2CL2" or $VOC eq "CH3CCL3" or $VOC eq "TCE" or $VOC eq "TRICLETH") {
            #print "$VOC : chlorinated\n";
            $topps{"Chlorinated"} += $value;
        } else {
            #print "$VOC : oxygenated\n";
            $topps{"Oxygenated"} += $value;
        }
    }
    return \%topps;
}

sub get_emission_fractions {
    my ($kpp, $mecca) = @_;
    my %emissions;

    my $emitting_reactions = $kpp->consuming("UNITY");
    for (0..$#$emitting_reactions) {
        my $reaction = $emitting_reactions->[$_];
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number);
        next if ($rate->sum == 0);
        my $VOC = $kpp->products($reaction)->[0];
        next if ($VOC eq "NO"); 
        if ($VOC =~ /^NC|^IC|^NEOP|^M2|^M3/ or $VOC eq "C3H8" or $VOC eq "CHEX" or $VOC eq "C2H6") {
            #print "$VOC : alkane\n";
            $emissions{"Alkanes"} += $rate->sum;
        } elsif ($VOC eq "C2H4" or $VOC eq "C3H6" or $VOC eq "C2H2" or $VOC eq "LIMONENE") {
            #print "$VOC : alkene or alkyne\n";
            $emissions{"Alkenes, Alkynes"} += $rate->sum;
        } elsif ($VOC =~ /^TOL|XYL|^TM|BENZ|^STY|ETHTOL$|^DIM/) {
            #print "$VOC : aromatic\n";
            $emissions{"Aromatics"} += $rate->sum;
        } elsif ($VOC eq "CH2CL2" or $VOC eq "CH3CCL3" or $VOC eq "TCE" or $VOC eq "TRICLETH") {
            #print "$VOC : chlorinated\n";
            $emissions{"Chlorinated"} += $rate->sum;
        } else {
            #print "$VOC : oxygenated\n";
            $emissions{"Oxygenated"} += $rate->sum;
        }
    }
    my $total_emissions = 0;
    $total_emissions += $emissions{$_} foreach (keys %emissions);
    $emissions{$_} /= $total_emissions foreach (keys %emissions);
    $emissions{"Total"} = $total_emissions;
    return \%emissions;
}
