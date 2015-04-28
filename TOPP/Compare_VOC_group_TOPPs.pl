#! /usr/bin/env perl
#Plot TOPPs of VOC groups, TOPP is mean TOPP from the same mechanism runs with each speciation
#Version 0: Jane Coates 28/4/2015

use strict;
use diagnostics;
use Statistics::R;

my $file = "mean_TOPPs_of_VOC_mechanism.csv";
my %TOPPs;
open my $in, '<:encoding(utf-8)', $file or die $!;
my @lines = <$in>;
close $in;

my @higher_alkanes = qw( NC6H14 M2PE M3PE NC7H16 M2HEX M3HEX NC8H18 NC9H20 NC10H22 NC11H24 NC12H26 CHEX );
my @butanes = qw( NC4H10 IC4H10 );
my @pentanes = qw( NC5H12 IC5H12 NEOP );
my @alkenes = qw( BUT1ENE BIGENE OL2 OLT OLI C10H16 C2H2 ISOP ISO C5H8 );
my @terpenes = qw( LIMONENE APINENE BPINENE );
my @aromatics = qw( EBENZ PBENZ IPBENZ STYRENE METHTOL OETHTOL PETHTOL DIME35EB );
my @xylenes = qw( MXYL OXYL PXYL );
my @tmbs = qw( TM123B TM124B TM135B );
my @chlorinated = qw( CH2CL2 CH3CCL3 TCE TRICLETH );
my @acids = qw(HCOOH CH3CO2H PROPACID CH3COOH );
my @alcohols = qw( CH3OH C2H5OH IPROPOL NBUTOL NPROPOL BUT2OL IBUTOL MIBKAOH C6H5CH2OH ETHGLY PROPGLY );
my @aldehydes = qw( HCHO CH3CHO C2H5CHO C3H7CHO IPRCHO C4H9CHO ACR MACR C4ALDB );
my @esters = qw( ETHACET NBUTACET IPROACET NPROACET );
my @ethers = qw( BUOX2ETOH PR2OHMOX EOX2EOL CH3OCH3 MO2EOL );
my @ketones = qw( CH3COCH3 MEK MIBK CYHEXONE ); 

foreach my $line (@lines) {
    next if ($line =~ /^Mechani/);
    chomp $line;
    my ($mechanism, $VOC, $TOPP) = split /\s/,$line;
    if ($VOC eq "C2H6") {
        $TOPPs{$mechanism}{"Ethane"} += $TOPP;
    } elsif ($VOC eq "C3H8") {
        $TOPPs{$mechanism}{"Propane"} += $TOPP;
    } elsif ($VOC ~~ @butanes) {
        $TOPPs{$mechanism}{"Butanes"} += $TOPP;
    } elsif ($VOC ~~ @pentanes) {
        $TOPPs{$mechanism}{"Pentanes"} += $TOPP;
    } elsif ($VOC ~~ @higher_alkanes) {
        $TOPPs{$mechanism}{"Higher Alkanes"} += $TOPP;
    } elsif ($VOC eq "C2H4") {
        $TOPPs{$mechanism}{"Ethene"} += $TOPP;
    } elsif ($VOC eq "C3H6") {
        $TOPPs{$mechanism}{"Propene"} += $TOPP;
    } elsif ($VOC ~~ @alkenes) {
        $TOPPs{$mechanism}{"Other Alkenes, Alkynes"} += $TOPP;
    } elsif ($VOC ~~ @terpenes) {
        $TOPPs{$mechanism}{"Terpenes"} += $TOPP;
    } elsif ($VOC eq "BENZENE") {
        $TOPPs{$mechanism}{"Benzene"} += $TOPP;
    } elsif ($VOC eq "TOLUENE") {
        $TOPPs{$mechanism}{"Toluene"} += $TOPP;
    } elsif ($VOC ~~ @xylenes) {
        $TOPPs{$mechanism}{"Xylenes"} += $TOPP;
    } elsif ($VOC ~~ @tmbs) {
        $TOPPs{$mechanism}{"Trimethylbenzenes"} += $TOPP;
    } elsif ($VOC ~~ @aromatics) {
        $TOPPs{$mechanism}{"Other Aromatics"} += $TOPP;
    } elsif ($VOC ~~ @chlorinated) {
        $TOPPs{$mechanism}{"Chlorinated"} += $TOPP;
    } elsif ($VOC ~~ @esters) {
        $TOPPs{$mechanism}{"Esters"} += $TOPP;
    } elsif ($VOC ~~ @ketones) {
        $TOPPs{$mechanism}{"Ketones"} += $TOPP;
    } elsif ($VOC ~~ @ethers) {
        $TOPPs{$mechanism}{"Ethers"} += $TOPP;
    } elsif ($VOC ~~ @alcohols) {
        $TOPPs{$mechanism}{"Alcohols"} += $TOPP;
    } elsif ($VOC ~~ @aldehydes) {
        $TOPPs{$mechanism}{"Aldehydes"} += $TOPP;
    } elsif ($VOC ~~ @acids) {
        $TOPPs{$mechanism}{"Acids"} += $TOPP;
    } else {
        print "No match for $VOC\n";
    } 
}

my $R = Statistics::R->new();
$R->run(q` library(tidyr) `,
        q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
        q` library(grid) `,
);
$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %TOPPs) {
    $R->set('mechanism', $mechanism);
    $R->run(q` pre = data.frame(Mechanism = mechanism) `);
    foreach my $VOC (sort keys %{$TOPPs{$mechanism}}) {
        $R->set('voc', $VOC);
        $R->set('topp', $TOPPs{$mechanism}{$VOC});
        $R->run(q` pre[voc] = topp `);
    }
    $R->run(q` pre = gather(pre, VOC, TOPP, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` plot = ggplot(data, aes(x = VOC, y = TOPP)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism, nrow = 4) `,
);

$R->run(q` CairoPDF(file = "TOPPs_each_type_by_mechanism.pdf", width = 10, height = 7) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = TOPP)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ VOC, nrow = 3) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + theme(axis.line = element_line(colour = "black")) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 0.7, vjust = 0.8)) `,
        q` plot = plot + theme(strip.text = element_text(face = "bold")) `,
        q` plot = plot + theme(panel.margin = unit(5, "mm")) `,
);

$R->run(q` CairoPDF(file = "TOPPs_each_type_by_VOC.pdf", width = 10, height = 7) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();
