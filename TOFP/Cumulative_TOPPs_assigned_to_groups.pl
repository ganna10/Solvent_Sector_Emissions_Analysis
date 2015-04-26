#! /usr/bin/env perl
# Extract cumulative TOPPs, assign to individual groups, facet by mechanism and speciation
# Version 0: Jane Coates 26/4/2015

use strict;
use diagnostics;
use Statistics::R;

#my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my %data;

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $TOPP_file = "${mechanism}_${speciation}_${speciation}_tagged_solvents_only_cumulative_TOPP_values.csv";
        open my $in, '<:encoding(utf-8)', $TOPP_file or die "Can't open $TOPP_file : $!";
        my @lines = <$in>;
        close $in;
        foreach my $line (@lines) {
            next if ($line =~ /^VOC/);
            chomp $line;
            my ($VOC, $TOPP) = split /,/, $line;
            my $category = get_category($VOC);
            $data{$mechanism}{$speciation}{$category} += $TOPP;
        }
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
        q` library(grid) `,
        q` library(ggthemes) `,
);

$R->run(q` data = data.frame() `);
foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->run(q` pre = data.frame(Mechanism = mechanism, Speciation = speciation) `);
        foreach my $category (sort keys %{$data{$mechanism}{$speciation}}) {
            $R->set('category', $category);
            $R->set('topp', $data{$mechanism}{$speciation}{$category});
            $R->run(q` pre[category] = topp `);
        }
        $R->run(q` pre = gather(pre, Category, TOPP, -Mechanism, -Speciation) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` my.colours = c("Acids" = "#cc6329", "Alcohols" = "#6c254f", "Benzene" = "#8c6238", "Butanes" = "#86b650", "Chlorinated" = "#f9c500", "CO" = "#898989", "Esters" = "#f3aa7f", "Ethane" = "#77aecc", "Ethene" = "#1c3e3d", "Ethers" = "#ba8b01", "Higher alkanes" = "#0e5c28", "Ketones" = "#ef6638", "Aldehydes" = "#8ed6d2", "Other alkenes, alkynes, dienes" = "#b569b3", "Other aromatics" = "#e7e85e", "Others" = "#2b9eb3", "Pentanes" = "#8c1531", "Propane" = "#9bb18d", "Propene" = "#623812", "Terpenes" = "#c9a415", "Toluene" = "#0352cb", "Trimethylbenzenes" = "#ae4901", "Xylenes" = "#1b695b", "CO" = "#6d6537", "Methane" = "#0c3f78", "Inorganic" = "#000000") `);
$R->run(q` plot.lines = function () { list( theme_tufte() ,
                                            ylab("Cumulative TOPP (molecules(Ox) molecules(VOC)-1)") ,
                                            scale_y_continuous(expand = c(0, 0)) ,
                                            scale_x_discrete(expand = c(0, 0)) ,
                                            theme(axis.ticks.x = element_blank()) ,
                                            theme(legend.title = element_blank()) ,
                                            theme(axis.line = element_line(colour = "black")) ,
                                            theme(axis.title.y = element_text(face = "bold")) ,
                                            theme(axis.title.x = element_blank()) ,
                                            theme(strip.text = element_text(face = "bold")) ,
                                            theme(axis.text.x = element_text(face = "bold", angle = 45, hjust = 0.7, vjust = 0.8)) ,
                                            scale_fill_manual(values = my.colours, limits = rev(levels(data$Category))) ,
                                            theme(panel.margin = unit(5, "mm")) ) } `);
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` data$Category = factor(data$Category, levels = c("Acids", "Alcohols", "Aldehydes", "Benzene", "Butanes", "Chlorinated", "Esters", "Ethane", "Ethene", "Ethers", "Higher alkanes", "Ketones", "Other alkenes, alkynes, dienes", "Other aromatics", "Pentanes", "Propane", "Propene", "Terpenes", "Toluene", "Trimethylbenzenes", "Xylenes")) `);

$R->run(q` plot = ggplot(data, aes(x = Speciation, y = TOPP, fill = Category)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism, nrow = 1) `,
        q` plot = plot + plot.lines() `,
);

$R->run(q` CairoPDF(file = "Cumulative_TOPP_categories_facet_mechanism.pdf", width = 8.7, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = TOPP, fill = Category)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2) `,
        q` plot = plot + plot.lines() `,
);

$R->run(q` CairoPDF(file = "Cumulative_TOPP_categories_facet_speciation.pdf", width = 8.7, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_category {
    my ($process) = @_;
    my $category;
    if ($process =~ /^TOL/) {
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
    } else {
        print "No category found for $process\n";
    }
    return $category;
}
