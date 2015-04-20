#! /usr/bin/env perl
# Get total emission rate of all VOC in speciations and allocate by molecular weight bins
# Version 0: Jane Coates 20/4/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions/MCM";
#my @speciations = qw( TNO UK08 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my %data;

my $weights_file = "MCM_molecular_weights.csv";
open my $in, '<:encoding(utf-8)', $weights_file or die "Can't open $weights_file : $!";
my @lines = <$in>;
close $in;
my %weights;
foreach my $line (@lines) {
    chomp $line;
    my ($species, $weight) = split /,/, $line;
    $weights{$species} = $weight;
}

foreach my $speciation (@speciations) {
    my $dir = "$base/${speciation}_Solvents_Only";
    my $boxmodel = "$dir/boxmodel";
    my $mecca = MECCA->new($boxmodel);
    my $eqn = "$dir/gas.eqn";
    my $kpp = KPP->new($eqn);
    $data{$speciation} = get_voc_emissions($mecca, $kpp, $speciation);
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
);

$R->run(q` data = data.frame() `);
foreach my $speciation (keys %data) {
    $R->set('speciation', $speciation);
    $R->run(q` pre = data.frame(Speciation = speciation) `);
    foreach my $VOC (sort keys %{$data{$speciation}}) {
        next if ($VOC eq "NO");
        print "No weight for $VOC\n" unless (defined $weights{$VOC});
        $R->set('voc', $VOC);
        $R->set('emissions', $data{$speciation}{$VOC});
        $R->set('weight', $weights{$VOC});
        $R->run(q` pre$VOC = voc `,
                q` pre$Emissions = emissions `,
                q` pre$Weight = weight `,
                q` data = rbind(data, pre) `,
        );
    }
}
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` data = data %>% mutate(Weight.Bins = ifelse(Weight < 40, "<40", ifelse(Weight >= 40 & Weight < 80, "40-80", ifelse(Weight >= 80 & Weight < 100, "80-100", ifelse(Weight >= 100 & Weight < 120, "100-120", ifelse(Weight >= 120 & Weight < 140, "120-140", ifelse(Weight >= 140 & Weight <180, "140-180", "other"))))))) `);
$R->run(q` data = data %>% group_by(Weight.Bins, Speciation) %>% summarize(Total.Emissions = sum(Emissions)) `);
$R->run(q` data = data %>% arrange(desc(Weight.Bins)) `);
$R->run(q` data$Weight.Bins = factor(data$Weight.Bins, levels = rev(c("<40", "40-80", "80-100", "100-120", "120-140", "140-180"))) `);
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` my.colours = c("<40" = "#62254f", "40-80" = "#f9c500", "80-100" = "#0e5628", "100-120" = "#ef6638", "120-140" = "#2b9eb3", "140-180" = "#b569b3") `);

$R->run(q` plot = ggplot(data, aes(x = Speciation, y = Total.Emissions, order = -as.numeric(Weight.Bins), fill = Weight.Bins)) `,
        q` plot = plot + geom_bar(stat = "identity", position = "stack") `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + ylab("Total VOC emissions (molecules cm-3)") `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.ticks.x = element_blank()) `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
        q` plot = plot + theme(axis.title = element_text(face = "bold")) `,
        q` plot = plot + theme(axis.text.x = element_text(face = "bold")) `,
        q` plot = plot + scale_fill_manual(values = my.colours) `,
);

$R->run(q` CairoPDF(file = "MCM_VOC_emissions_vs_Weights.pdf") `,
        q` print(plot) `,
        q` dev.off() `,
);

###Correlate total VOC emissions with total NO emissions
foreach my $speciation (keys %data) {
    my $total_VOC = 0;
    foreach my $VOC (keys %{$data{$speciation}}) {
        next if ($VOC eq "NO");
        $total_VOC += $data{$speciation}{$VOC};
        delete $data{$speciation}{$VOC};
    }
    $data{$speciation}{"Total.VOC"} = $total_VOC;
}

$R->run(q` no.data = data.frame() `);
foreach my $speciation (keys %data) {
    $R->set('speciation', $speciation);
    $R->run(q` pre = data.frame(Speciation = speciation) `);
    foreach my $item (sort keys %{$data{$speciation}}) {
        $R->set('emissions', $data{$speciation}{$item});
        if ($item eq "NO") {
            $R->run(q` pre$NO.emissions = emissions `);
        } elsif ($item =~ "Total") {
            $R->run(q` pre$VOC.emissions = emissions `);
        }
    }
    $R->run(q` no.data = rbind(no.data, pre) `);
}
#my $p = $R->run(q` print(no.data) `);
#print $p, "\n";
$R->run(q` speciation.colours = c("TNO" = "#000000", "IPCC" = "#2b9eb3", "EMEP" = "#ef6638", "DE94" = "#f9c500", "GR95" = "#6c254f", "GR05" = "#0352cb", "UK98" = "#0e5c28", "UK08" = "#b569b3") `);
$R->run(q` lm_eqn = function(m) { l <- list(a = format(coef(m)[1], digits = 2),
                                      b = format(abs(coef(m)[2]), digits = 2),
                                      r2 = format(summary(m)$r.squared, digits = 3));
                                  if (coef(m)[2] >= 0)  {
                                    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
                                  } else {
                                    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
                                  }
                                  as.character(as.expression(eq));                 
                                } `);

$R->run(q` p = ggplot(no.data, aes(x = VOC.emissions, y = NO.emissions, colour = Speciation)) `,
        q` p = p + geom_point(size = 4) `,
        q` p = p + theme_tufte() `,
        q` p = p + theme(axis.line = element_line(colour = "black")) `,
        q` p = p + geom_text(aes(x = 8e8, y = 2e9, label = lm_eqn(lm(VOC.emissions ~ NO.emissions, no.data))), colour = "black", parse = TRUE) + geom_smooth(method = "lm", se = FALSE, colour = "black") `,
        q` p = p + scale_colour_manual(values = speciation.colours) `,
);

$R->run(q` CairoPDF(file = "NO_vs_total_VOC_emissions.pdf") `,
        q` print(p) `,
        q` dev.off() `,
);

$R->stop();

sub get_voc_emissions {
    my ($mecca, $kpp, $speciation) = @_;

    my @no_match = qw( NO2 CO CH4 );
    my $emissions = $kpp->consuming("UNITY");
    my %emissions;
    for (0..$#$emissions) {
        my $reaction = $emissions->[$_];
        my $products = $kpp->products($reaction);
        my $VOC = $products->[0];
        next if ($VOC ~~ @no_match);
        my $reaction_number = $kpp->reaction_number($reaction);
        my $rate = $mecca->rate($reaction_number);
        next if ($rate->sum == 0);
        $emissions{$VOC} += $rate->sum;
    }
    return \%emissions;
}
