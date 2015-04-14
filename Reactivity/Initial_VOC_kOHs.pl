#! /usr/bin/env perl
# Sum total kOH of initial VOC reactions
# Version 0: Jane Coates 14/4/2015
# Version 1: Jane Coates 14/4/2015 using sum of emitted VOC concentrations between 06:00 and 12:00 of first day, i.e. when emissions are constant

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

#my $base = "/local/home/coates/Solvent_Emissions";
my $base = "/work/users/jco/Solvent_Emissions";
my @mechanisms = qw( MCM );
#my @speciations = qw( EMEP );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my %data;

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $eqn = "$base/$mechanism/${speciation}_Solvents_Only/gas.eqn";
        my $kpp = KPP->new($eqn);
        my $boxmodel = "$base/$mechanism/${speciation}_Solvents_Only/boxmodel";
        my $mecca = MECCA->new($boxmodel);
        my $label;
        if ($mechanism eq "RADM2") {
            $label = "RADM";
        } else {
            $label = $mechanism;
        }
        my $emis_file = "$base/$mechanism/${speciation}_Solvents_Only/boxmodel/${label}_EMIS.nml";
        $data{$mechanism}{$speciation} = get_data($kpp, $mecca, $emis_file);
    }
}

my $R = Statistics::R->new();
$R->run(q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(tidyr) `,
        q` library(ggthemes) `,
);
$R->run(q` data = data.frame() `);

foreach my $mechanism (sort keys %data) {
    $R->set('mechanism', $mechanism);
    #print "$mechanism\n";
    foreach my $speciation (sort keys %{$data{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->run(q` pre = data.frame(Mechanism = mechanism) `);
        $R->run(q` pre$Speciation = speciation `);
        #print "\t$speciation\n";
        foreach my $process (sort keys %{$data{$mechanism}{$speciation}}) {
            $R->set('process', $process);
            $R->set('koh', $data{$mechanism}{$speciation}{$process});
            $R->run(q` pre[process] = koh `);
            #print "\t\t$VOC : $data{$mechanism}{$speciation}{$VOC}\n";
        }
        $R->run(q` pre = gather(pre, Process, kOH, -Mechanism, -Speciation) `,
                q` data = rbind(data, pre) `,
        );
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";

$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);
$R->run(q` data$Process = factor(data$Process, levels = c("Acids", "Alcohols", "Aldehydes", "Benzene", "Butanes", "Chlorinated", "Emissions", "Esters", "Ethane", "Ethene", "Ethers", "Higher alkanes", "Ketones", "Methane", "Other alkenes, alkynes, dienes", "Other aromatics", "Pentanes", "Propane", "Propene", "Terpenes", "Toluene", "Trimethylbenzenes", "Xylenes")) `);
$R->run(q` my.colours = c("Acids" = "#cc6329", "Alcohols" = "#6c254f", "Benzene" = "#8c6238", "Butanes" = "#86b650", "Chlorinated" = "#f9c500", "Esters" = "#f3aa7f", "Ethane" = "#77aecc", "Ethene" = "#1c3e3d", "Ethers" = "#ba8b01", "Higher alkanes" = "#0e5c28", "Ketones" = "#ef6638", "Aldehydes" = "#8ed6d2", "Other alkenes, alkynes, dienes" = "#b569b3", "Other aromatics" = "#e7e85e", "Others" = "#2b9eb3", "Pentanes" = "#8c1531", "Propane" = "#9bb18d", "Propene" = "#623812", "Terpenes" = "#c9a415", "Toluene" = "#0352cb", "Trimethylbenzenes" = "#ae4901", "Xylenes" = "#1b695b", "CO" = "#6d6537", "Methane" = "#0c3f78", "Emissions" = "#000000") `);

$R->run(q` plot = ggplot(data, aes(x = Speciation, y = kOH, fill = Process)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism, nrow = 1) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
        q` plot = plot + theme(axis.ticks.x = element_blank()) `,
        q` plot = plot + theme(axis.title.x = element_blank()) `,
        q` plot = plot + theme(legend.title = element_blank()) `,
        q` plot = plot + ylab("OH Reactivity (s-1)") `,
        q` plot = plot + scale_fill_manual(values = my.colours, limits = rev(levels(data$Process))) `,
);

$R->run(q` CairoPDF(file = "initial_VOC_kOH.pdf", width = 8.6, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_rate_constant {
    my ($reaction, $kpp, $cair) = @_;
    my $rate_constant = $kpp->rate_string($reaction);
    #print $rate_constant, "\n";
    if ($rate_constant =~ /^ARR2/) {
        $rate_constant =~ s/ARR2\(|_dp|, TEMP\)//g;
        my ($factor, $exponent) = split /, /, $rate_constant;
        if ($exponent =~ /^-/) {
            $exponent =~ s/^-|\s$//g;
            $rate_constant = "$factor*exp($exponent/TEMP)";
        } else {
            $exponent =~ s/\s$//g;
            $rate_constant = "$factor*exp(-$exponent/TEMP)";
        }
    } elsif ($rate_constant =~ /^THER/) {
        $rate_constant =~ s/THERMAL_T2\(|_dp, TEMP\)//g;
        my ($factor, $exponent) = split /, /, $rate_constant;
        if ($exponent =~ /^-/) {
            $exponent =~ s/^-|\s$//g;
            $rate_constant = "(TEMP**2)*$factor*exp($exponent/TEMP)";
        } else {
            $exponent =~ s/\s$//g;
            $rate_constant = "(TEMP**2)*$factor*exp(-$exponent/TEMP)";
        }
    } elsif ($rate_constant =~ /RC2H2OH/) {
        my $K0C2H2 = 5.5e-30;
        my $KIC2H2 = 8.3e-13*(293/300)**2;
        my $AC2H2 = $K0C2H2*$cair/$KIC2H2;
        my $BC2H2 = log10 $AC2H2;
        my $MC2H2 = $K0C2H2*$cair/(1+$AC2H2);
        my $EC2H2 = (1+$BC2H2**2)**(-1);
        $rate_constant = $MC2H2*0.6**$EC2H2;
    } elsif ($rate_constant =~ /RC2H4OH/) {
        my $K0C2H4 = 1.0e-28*(293/300)**(-0.8);
        my $KIC2H4 = 8.8e-12;
        my $AC2H4 = $K0C2H4*$cair/$KIC2H4;
        my $BC2H4 = log10($AC2H4);
        my $MC2H4 = $K0C2H4*$cair/(1+$AC2H4);
        my $EC2H4 = (1+$BC2H4**2)**(-1);
        $rate_constant = $MC2H4*0.6**$EC2H4;
    } elsif ($rate_constant =~ /RC3H6OH/) {
        my $K0C3H6 = 8.0e-27*(293/300)**(-3.5);
        my $KIC3H6 = 3.0e-11;
        my $AC3H6 = $K0C3H6*$cair/$KIC3H6;
        my $BC3H6 = log10($AC3H6);
        my $MC3H6 = $K0C3H6*$cair/(1+$AC3H6);
        my $EC3H6 = (1+$BC3H6**2)**(-1);
        $rate_constant = $MC3H6*0.5**$EC3H6;
    } elsif ($rate_constant =~ /RMPANOH/) {
        my $K0MPAN = 8.0e-27*(293/300)**(-3.5);
        my $KIMPAN = 3.0e-11;
        my $AMPAN = $K0MPAN*$cair/$KIMPAN;
        my $BMPAN = log10($AMPAN);
        my $MMPAN = $K0MPAN*$cair/(1+$AMPAN);
        my $EMPAN = (1+$BMPAN**2)**(-1);
        $rate_constant = $MMPAN*0.5**$EMPAN;
    } elsif ($rate_constant =~ /RHCNOH/) {
        my $K0HCN = 4.28e-33;
        my $KIHCN = 9.3e-15*(293/300)**4.42;
        my $AHCN = $K0HCN*$cair/$KIHCN;
        my $BHCN = log10($AHCN);
        my $MHCN = $K0HCN*$cair/(1+$AHCN);
        my $EHCN = (1+$BHCN**2)**(-1);
        $rate_constant = $MHCN*0.8**$EHCN;
    } elsif ($rate_constant =~ /KMT05/) {
        my ($KMT, $factor) = split /KMT05/, $rate_constant;
        $KMT = 1.44e-13*(1+($cair/4.2e+19));
        if (defined $factor) {
            $rate_constant = $KMT . $factor;
        } else {
            $rate_constant = $KMT;
        } 
    } elsif ($rate_constant =~ /KMT07/) {
        my ($KMT, $factor) = split /KMT07/, $rate_constant;
        my $K70 = 7.4e-31*$cair*(293/300)**(-2.4);
        my $K7I = 3.3e-11*(293/300)**(-0.3);
        my $FC7 = exp(-293/1420);
        my $KR7 = $K70/$K7I;
        my $NC7 = 0.75-1.27*(log10($FC7));
        my $F7 = 10**(log10($FC7)/(1+(log10($KR7)/$NC7)**(2)));
        $KMT = ($K70*$K7I)*$F7/($K70+$K7I);
        if (defined $factor) {
            $rate_constant = $KMT . $factor;
        } else {
            $rate_constant = $KMT;
        } 
    } elsif ($rate_constant =~ /KMT08/) {
        my ($KMT, $factor) = split /KMT08/, $rate_constant;
        my $K80 = 3.3e-30*$cair*(293/300)**(-3.0);
        my $K8I = 4.1e-11;
        my $FC8 = 0.4;
        my $KR8 = $K80/$K8I;
        my $NC8 = 0.75-1.27*(log10($FC8));
        my $F8 = 10**(log10($FC8)/(1+(log10($KR8)/$NC8)**(2)));
        $KMT = ($K80*$K8I)*$F8/($K80+$K8I);
        if (defined $factor) {
            $rate_constant = $KMT . $factor;
        } else {
            $rate_constant = $KMT;
        } 
    } elsif ($rate_constant =~ /KMT11/) {
        my ($KMT, $factor) = split /KMT11/, $rate_constant;
        my $K1 = 2.40e-14*exp(460/293);
        my $K3 = 6.50e-34*exp(1335/293);
        my $K4 = 2.70e-17*exp(2199/293);
        my $K2 = ($K3*$cair)/(1+($K3*$cair/$K4));
        $KMT = $K1 + $K2;
        if (defined $factor) {
            $rate_constant = $KMT . $factor;
        } else {
            $rate_constant = $KMT;
        } 
    } elsif ($rate_constant =~ /KMT12/) {
        my ($KMT, $factor) = split /KMT12/, $rate_constant;
        my $K120 = 4.5e-31*$cair*(293/300)**(-3.9);
        my $K12I = 1.3e-12*(293/300)**(-0.7);
        my $FC12 = 0.525;
        my $NC12 = 0.75-1.27*(log10($FC12));
        my $KR12 = $K120/$K12I;
        my $F12 = 10**(log10($FC12)/(1.0+(log10($KR12)/$NC12)**(2)));
        $KMT = ($K120*$K12I*$F12)/($K120+$K12I);
        if (defined $factor) {
            $rate_constant = $KMT . $factor;
        } else {
            $rate_constant = $KMT;
        } 
    } elsif ($rate_constant =~ /KMT15/) {
        my ($KMT, $factor) = split /KMT15/, $rate_constant;
        my $K150 = 8.6e-29*$cair*(293/300)**(-3.1);
        my $K15I = 9.0e-12*(293/300)**(-0.85);
        my $KR15 = $K150/$K15I;
        my $FC15 = 0.48;
        my $NC15 = 0.75-1.27*(log10($FC15));
        my $F15 = 10**(log10($FC15)/(1+(log10($KR15)/$NC15)**(2)));
        $KMT = ($K150*$K15I)*$F15/($K150+$K15I);
        if (defined $factor) {
            $rate_constant = $KMT . $factor;
        } else {
            $rate_constant = $KMT;
        } 
    } elsif ($rate_constant =~ /KMT16/) {
        my ($KMT, $factor) = split /KMT16/, $rate_constant;
        my $K160 = 8e-27*$cair*(293/300)**(-3.5);
        my $K16I = 3.0e-11*(293/300)**(-1);
        my $KR16 = $K160/$K16I;
        my $FC16 = 0.5;
        my $NC16 = 0.75-1.27*(log10($FC16));
        my $F16 = 10**(log10($FC16)/(1+(log10($KR16)/$NC16)**(2)));
        $KMT = ($K160*$K16I)*$F16/($K160+$K16I);
        if (defined $factor) {
            $rate_constant = $KMT . $factor;
        } else {
            $rate_constant = $KMT;
        } 
    } elsif ($rate_constant =~ /KMT17/) {
        my ($KMT, $factor) = split /KMT17/, $rate_constant;
        my $K170 = 5.0e-30*$cair*(293/300)**(-1.5);
        my $K17I = 1.0e-12;
        my $KR17 = $K170/$K17I;
        my $FC17 = 0.17*exp(-51/293)+exp(-293/204);
        my $NC17 = 0.75-1.27*(log10($FC17));
        my $F17 = 10**(log10($FC17)/(1.0+(log10($KR17)/$NC17)**(2)));
        $KMT = ($K170*$K17I*$F17)/($K170+$K17I);
        if (defined $factor) {
            $rate_constant = $KMT . $factor;
        } else {
            $rate_constant = $KMT;
        } 
    } elsif ($rate_constant =~ /KMT18/) {
        my ($KMT, $factor) = split /KMT18/, $rate_constant;
        $KMT = 9.5e-39*210e-03*$cair*exp(5270/293)/(1+7.5e-29*210e-03*$cair*exp(5610/293));
        if (defined $factor) {
            $rate_constant = $KMT . $factor;
        } else {
            $rate_constant = $KMT;
        } 
    } elsif ($rate_constant =~ /_VD\(KPP_O3\)/) {
        $rate_constant = 0.5/(1000*100);
    }
    $rate_constant =~ s/EXP/exp/g;
    $rate_constant =~ s/TEMP/293/g;
    $rate_constant =~ s/D/e/g;
    #print $rate_constant, "\n";
    return $rate_constant;
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

sub get_data {
    my ($kpp, $mecca, $emis_file) = @_;
    my $cair = $mecca->cair->at(0);

    open my $in, '<:encoding(utf-8)', $emis_file or die "Can't open $emis_file : $!";
    my @lines = <$in>;
    close $in;
    my @VOC;
    foreach my $line (@lines) {
        next unless ($line =~ /^EMIS/);
        chomp $line;
        $line =~ s/^EMIS_(.*?)=(.*?)$/$1/;
        push @VOC, $line;
    }

    my %kOH;
    foreach my $VOC (@VOC) {
        my $concentration = $mecca->tracer($VOC) * $cair;
        #$concentration = $concentration(1:19);
        my $oh_reactions = $kpp->reacting_with($VOC, 'OH');
        for (0..$#$oh_reactions) {
            my $reaction = $oh_reactions->[$_];
            my $rate_constant = get_rate_constant($reaction, $kpp, $cair);
            my $category = get_category($VOC);
            $kOH{$category} += eval($rate_constant) * $concentration->sum;
        }
    }
    return \%kOH;
}
