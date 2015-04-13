#! /usr/bin/env perl
# Calculate initial total reactivity of initial VOC with OH i.e. k(VOC + OH) x emission rate
# Version 0: Jane Coates 13/4/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/work/users/jco/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
#my @speciations = qw( IPCC );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my %data;
my @inorganic = qw( NO NO2 O3 );

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        my $dir = "$base/$mechanism/${speciation}_Solvents_Only";
        my $boxmodel = "$dir/boxmodel";
        my $mecca = MECCA->new($boxmodel);
        my $eqn = "$dir/gas.eqn";
        my $kpp = KPP->new($eqn);
        $data{$mechanism}{$speciation} = get_data($mecca, $kpp);
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
    $R->run(q` pre = data.frame(Mechanism = mechanism) `);
    foreach my $speciation (keys %{$data{$mechanism}}) {
        $R->set('speciation', $speciation);
        $R->set('oh.reactivity', $data{$mechanism}{$speciation});
        $R->run(q` pre[speciation] = oh.reactivity `);
    }
    $R->run(q` pre = gather(pre, Speciation, OH.Reactivity, -Mechanism) `,
            q` data = rbind(data, pre) `,
    );
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);

$R->run(q` plot = ggplot(data, aes(x = Speciation, y = OH.Reactivity)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Mechanism, nrow = 1) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
);

$R->run(q` CairoPDF(file = "Total_OH_reactivity_facet_Mechanism.pdf", width = 8.7, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->run(q` plot = ggplot(data, aes(x = Mechanism, y = OH.Reactivity)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_wrap( ~ Speciation, nrow = 2) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
);

$R->run(q` CairoPDF(file = "Total_OH_reactivity_facet_Speciation.pdf", width = 8.7, height = 6) `,
        q` print(plot) `,
        q` dev.off() `,
);

$R->stop();

sub get_data {
    my ($mecca, $kpp) = @_;
    my $total_reactivity;
    my $cair = $mecca->cair->at(0);
    my @oxidants = qw( OH O3 );
    foreach my $oxidant (@oxidants) {
        my $consumers = $kpp->consuming($oxidant);
        my $reactivity;
        foreach my $reaction (@$consumers) {
            my $reactants = $kpp->reactants($reaction);
            my ($other_reactant) = grep { $_ ne $oxidant } @$reactants;
            next unless (defined $other_reactant);
            next if ($other_reactant eq 'hv');
            next if ($other_reactant ~~ @inorganic);
            my $other_reactant_conc = $mecca->tracer($other_reactant) * $cair;
            my $rate_constant = get_rate_constant($reaction, $kpp, $cair);
            $rate_constant =  eval $rate_constant;
            $reactivity += $rate_constant * $other_reactant_conc;
            #print $reactivity->sum, "\n";
        }
        $total_reactivity += $reactivity->sum;
    }
    #print $total_reactivity, "\n";
    return $total_reactivity;
}

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
