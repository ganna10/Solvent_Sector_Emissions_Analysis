#! /usr/bin/env perl
# Plot concentrations (molecule cm-3) of initial VOC in each speciation and mechanism
# Version 0: Jane Coates 14/4/2015

use strict;
use diagnostics;
use MECCA;
use Statistics::R;
use PDL;
use PDL::NiceSlice;

my $base = "/work/users/jco/Solvent_Emissions";
my @mechanisms = qw( MCM );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my %concentrations;

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $boxmodel = "$base/$mechanism/${speciation}_Solvents_Only/boxmodel";
        my $mecca = MECCA->new($boxmodel);
        my $label;
        if ($mechanism eq "RADM2") {
            $label = "RADM";
        } else {
            $label = $mechanism;
        }
        my $emis_file = "$boxmodel/${label}_EMIS.nml";
        my @VOC = get_VOC($emis_file);
        my $cair = $mecca->cair->at(0);
        foreach my $VOC (@VOC) {
            my $tracer = $mecca->tracer($VOC);
            my $category = get_category($VOC);
            $concentrations{$mechanism}{$speciation}{$category} += $tracer * $cair;
        }
    }
}
my $mecca = MECCA->new("$base/MCM/EMEP_Solvents_Only/boxmodel");
my $times = $mecca->time;
$times -= $times->at(0);
$times /= 86400;

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
        q` library(ggthemes) `,
);
$R->run(q` data = data.frame() `);
$R->set('Time', [map { $_ } $times->dog]);

foreach my $mechanism (sort keys %concentrations) {
    $R->set('mechanism', $mechanism);
    foreach my $speciation (sort keys %{$concentrations{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->run(q` pre = data.frame(Time) `,
                q` pre$Speciation = rep(speciation, length(Time)) `,
                q` pre$Mechanism = rep(mechanism, length(Time)) `,
        );
        foreach my $process (sort keys %{$concentrations{$mechanism}{$speciation}}) {
            $R->set('process', $process);
            $R->set('concentration', [ map { $_ } $concentrations{$mechanism}{$speciation}{$process}->dog ]);
            $R->run(q` pre[process] = concentration `);
        }
        $R->run(q` pre = gather(pre, Process, Concentration, -Mechanism, -Speciation, -Time) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` data$Process = factor(data$Process, levels = c("Acids", "Alcohols", "Aldehydes", "Benzene", "Butanes", "Chlorinated", "Emissions", "Esters", "Ethane", "Ethene", "Ethers", "Higher alkanes", "Ketones", "Methane", "Other alkenes, alkynes, dienes", "Other aromatics", "Pentanes", "Propane", "Propene", "Terpenes", "Toluene", "Trimethylbenzenes", "Xylenes")) `);

$R->run(q` plot = ggplot(data, aes(x = Time, y = Concentration)) `,
        q` plot = plot + geom_line() `,
        q` plot = plot + facet_grid( Speciation ~ Process, scales = "free" ) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + scale_x_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
);

$R->run(q` CairoPDF(file = "Concentrations_initial_VOCs.pdf", width = 38.7, height = 26) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_VOC {
    my ($file) = @_;
    my @VOC;
    open my $in, '<:encoding(utf-8)', $file or die "Can't open $file : $!";
    my @lines = <$in>;
    close $in;
    foreach my $line (@lines) {
        next unless ($line =~ /^EMIS_/);
        chomp $line;
        $line =~ s/^EMIS_(.*?)=(.*?)$/$1/;
        push @VOC, $line;
    }
    return @VOC;
}

sub get_category {
    my ($process, $mechanism, $speciation) = @_;
    my $category;
    if ($process =~ /CO \+ OH/) {
        $category = "CO";
    } elsif ($process =~ /=/) {
        $category = "Inorganic";
    } elsif ($process eq "CH4") {
        $category = "Methane";
    } elsif ($process =~ /^TOL/) {
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
    } elsif ($process =~ /Deposition|Inorganic|Emissions/) {
        $category = $process;
    } else {
        print "No category found for $process in $mechanism, $speciation\n";
    }
    return $category;
}
