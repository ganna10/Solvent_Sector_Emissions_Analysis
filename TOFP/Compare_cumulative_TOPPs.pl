#! /usr/bin/env perl
# Compare TOPPs of individual VOC between speciations and mechanisms
# Version 0: Jane Coates 27/4/2015

use strict;
use diagnostics;
use Statistics::R;

my @mechanisms = qw( MCM MOZART RADM2 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my %data;

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $TOPP_file = "${mechanism}_${speciation}_${speciation}_tagged_solvents_only_cumulative_TOPP_values.csv";
        open my $in, '<:encoding(utf-8)', $TOPP_file or die $!;
        my @lines = <$in>;
        close $in;
        foreach my $line (@lines) {
            next if ($line =~ /^VOC/);
            chomp $line;
            my ($voc, $topp) = split /,/, $line;
            $data{$mechanism}{$speciation}{$voc} = $topp;
        }
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(tidyr) `,
        q` library(grid) `,
        q` library(scales) `,
        q` library(ggthemes) `,
        q` library(Cairo) `,
);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->run(q` pre = data.frame(Mechanism = mechanism, Speciation = speciation) `);
        foreach my $VOC (sort keys %{$data{$mechanism}{$speciation}}) {
            $R->set('voc', $VOC);
            $R->set('topp', $data{$mechanism}{$speciation}{$VOC});
            $R->run(q` pre[voc] = topp `);
        }
        $R->run(q` pre = gather(pre, VOC, TOPP, -Mechanism, -Speciation) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` mechanism.colours = c("MCM" = "#6c254f", "MOZART" = "#ef6638", "RADM2" = "#0e5c28") `);
$R->run(q` speciation.colours = c("TNO" = "#000000", "IPCC" = "#2b9eb3", "EMEP" = "#ef6638", "DE94" = "#f9c500", "GR95" = "#6c254f", "GR05" = "#0352cb", "UK98" = "#0e5c28", "UK08" = "#b569b3") `);
$R->run(q` plot.lines = function () { list( theme_bw() ,
                                            ylab("Ox Production (molecules cm-3)") ,
                                            theme(legend.title = element_blank()) ,
                                            theme(axis.line = element_line(colour = "black")) ,
                                            theme(axis.title = element_text(face = "bold")) ,
                                            theme(strip.text = element_text(face = "bold")) ,
                                            theme(panel.margin = unit(5, "mm")) ) } `);

$R->run(q` plot = ggplot(data, aes(x = Speciation, y = TOPP, colour = Mechanism)) `,
        q` plot = plot + geom_point(size = 4) `,
        q` plot = plot + facet_wrap( ~ VOC) `,
        q` plot = plot + plot.lines() `,
        q` plot = plot + scale_colour_manual(values = mechanism.colours) `,
);

$R->run(q` CairoPDF(file = "TOPP_comparison_betwee_mechanisms.pdf", width = 15, height = 11) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = TOPP, colour = Speciation)) `,
        q` plot = plot + geom_point(size = 4) `,
        q` plot = plot + facet_wrap( ~ VOC) `,
        q` plot = plot + plot.lines() `,
        q` plot = plot + scale_colour_manual(values = speciation.colours) `,
);

$R->run(q` CairoPDF(file = "TOPP_comparison_betwee_speciation.pdf", width = 15, height = 11) `,
        q` print(plot) `,
        q` dev.off() `,
); 

$R->stop();
