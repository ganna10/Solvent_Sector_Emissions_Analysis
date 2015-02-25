#! /usr/bin/env perl
# Plot different attributions of the solvent speciations in each mechanism
# Version 0: Jane Coates 25/2/2015

use strict;
use diagnostics;
use Statistics::R;

my (%emissions, %groups);
my $base = "/work/users/jco/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
my @speciations = qw( DE94 EMEP GR05 GR95 IPCC TNO UK08 UK98 );
my @functional_gps = qw( Alkanes Alkenes Aromatics Aldehydes Alcohols Acids Ketones Terpenes Alkynes Ethers Esters Chlorinated );

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $boxmodel = "$base/$mechanism/${speciation}_Solvents_Only/boxmodel";
        my $emis_name;
        if ($mechanism eq "RADM2") {
            ($emis_name = $mechanism) =~ s/2//;
        } else {
            $emis_name = $mechanism;
        }
        my $emis_file = "$boxmodel/${emis_name}_EMIS.nml";
        my $emissions = get_emissions($emis_file);
        $emissions{$mechanism}{$speciation} = $emissions;
    }
}

foreach my $mechanism (sort keys %emissions) {
    foreach my $speciation (sort keys %{$emissions{$mechanism}}) {
        my $total_emissions = 0;
        $total_emissions += $emissions{$mechanism}{$speciation}{$_} foreach (keys %{$emissions{$mechanism}{$speciation}});
        foreach my $VOC (sort keys %{$emissions{$mechanism}{$speciation}}) {
            my $group = return_group($VOC);
            $groups{$mechanism}{$speciation}{$group} += $emissions{$mechanism}{$speciation}{$VOC} / $total_emissions;
        }
    }
}

foreach my $mechanism (sort keys %groups) {
    foreach my $speciation (sort keys %{$groups{$mechanism}}) {
        foreach (@functional_gps) {
            next if (exists $groups{$mechanism}{$speciation}{$_});
            $groups{$mechanism}{$speciation}{$_} = 0.00;
        }
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
        q` library(tidyr) `,
        q` library(scales) `,
);

$R->run(q` my.colours = c("Alkanes" = "#6c254f", "Alkenes" = "#f9c500", "Aromatics" = "#2b9eb3", "Aldehydes" = "#86b650", "Alcohols" = "#0352cb", "Acids" = "#cc6329", "Chlorinated" = "#8c1531", "Ketones" = "#c9a415", "Terpenes" = "#0e5c28", "Alkynes" = "#1c3e3d", "Ethers" = "#b569b3", "Esters" = "#ef6638") `);

foreach my $mechanism (sort keys %groups) {
    $R->set('filename', "${mechanism}_emissions_by_functional_gp.pdf");
    $R->set('title', "$mechanism: Percentage of Emissions by Functional Group");
    $R->run(q` data = data.frame() `);
    foreach my $speciation (sort keys %{$groups{$mechanism}}) {
        $R->run(q` pre = data.frame(dummy1 = c(1)) `);
        $R->set('speciation', $speciation);
        foreach my $gp (sort keys %{$groups{$mechanism}{$speciation}}) {
            $R->set('group', $gp);
            $R->set('ratio', $groups{$mechanism}{$speciation}{$gp});
            $R->run(q` pre[group] = ratio `);
        }
        $R->run(q` pre = gather(pre, Group, Ratio) `,
                q` pre = pre[-1,] `,
                q` pre$Speciation = rep(speciation, length(pre$Group)) `,
                q` data = rbind(data, pre) `,
        );
    }
    #my $p = $R->run(q` print(levels(data$Group)) `);
    #print $p, "\n";
    $R->run(q` data$Group = factor(data$Group, levels = c("Acids", "Alcohols", "Aldehydes", "Alkanes", "Alkenes", "Alkynes", "Aromatics", "Chlorinated", "Esters", "Ethers", "Ketones", "Terpenes")) `);
    $R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);

    $R->run(q` plot = ggplot(data, aes(x = Speciation, y = Ratio, fill = Group)) `,
            q` plot = plot + geom_bar(stat = "identity", position = "stack") `,
            q` plot = plot + ggtitle(title) `,
            q` plot = plot + theme_tufte() `,
            q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
            q` plot = plot + scale_y_continuous(labels = percent, expand = c(0, 0)) `,
            q` plot = plot + theme(axis.title = element_blank()) `,
            q` plot = plot + theme(legend.title = element_blank()) `,
            q` plot = plot + theme(axis.ticks.x = element_blank()) `,
            q` plot = plot + theme(plot.title = element_text(face = "bold")) `,
            q` plot = plot + theme(axis.text.x = element_text(face = "bold")) `,
            q` plot = plot + scale_fill_manual(values = my.colours, limits = rev(levels(data$Group))) `,
    );

    $R->run(q` CairoPDF(file = filename, width = 8.6, height = 6) `,
            q` print(plot) `,
            q` dev.off() `,
    );
}

$R->stop();

sub return_group {
    my ($VOC) = @_;
    my $group;
    if ($VOC =~ /^HC\d|^ETH$|BIGALK|^NC|^IC|C2H6|C3H8|NEOP|CHEX|M2PE|M3PE|M3HEX|M2HEX/) {
        $group = "Alkanes";
    } elsif ($VOC =~ /^OL|^ISO|BIGENE|C2H4|C3H6|C5H8|1ENE$/) {
        $group = "Alkenes";
    } elsif ($VOC =~ /^TOL|ETHTOL$|XYL|ENZ|^TM|DIME35EB|STYRENE/) {
        $group = "Aromatics";
    } elsif ($VOC =~ /HCHO|ALD|^CH2O|CHO$|ACR$/) {
        $group = "Aldehydes";
    } elsif ($VOC =~ /C2H2/) {
        $group = "Alkynes";
    } elsif ($VOC =~ /KET|MEK|CH3COCH3|MIBK|CYHEXONE/) {
        $group = "Ketones";
    } elsif ($VOC =~ /ORA|HCOOH|CH3COOH|CH3CO2H|ACID/) {
        $group = "Acids";
    } elsif ($VOC =~ /C2H5OH|CH3OH|C6H5CH2OH|IBUTOL|NBUTOL|2OL$|POL$|MIBKAOH|GLY$/) {
        $group = "Alcohols";
    } elsif ($VOC =~ /C10H16|LIMONENE|PINENE/) {
        $group = "Terpenes";
    } elsif ($VOC =~ /TRICLETH|CH3CCL3|CH2CL2|TCE/) {
        $group = "Chlorinated";
    } elsif ($VOC =~ /ACET$/) {
        $group = "Esters";
    } elsif ($VOC =~ /CH3OCH3|EOL$|PR2OHMOX|BUOX2ETOH/) {
        $group = "Ethers";
    } else {
        print "No group for $VOC\n";
    }
    #print "$VOC = $group\n";
    return $group;
}

sub get_emissions {
    my ($file) = @_;

    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file : $!";
    local $/ = undef;
    my $data = <$in>;
    my %data;
    my @lines = split /\s+/, $data;
    foreach my $line (@lines) {
        next unless ($line =~ /^EMIS_/);
        $line =~ s/EMIS_|,//g;
        my ($species, $emission) = split /=/, $line;
        $data{$species} = $emission;
    } 
    close $in;
    return \%data;
}
