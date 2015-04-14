#! /usr/bin/env perl
# Version 0: Jane Coates 13/4/2015

use strict;
use diagnostics;
use KPP;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $base = "/local/home/coates/Solvent_Emissions";
#my $base = "/work/users/jco/Solvent_Emissions";
my @mechanisms = qw( MCM MOZART RADM2 );
#my @speciations = qw( GR05 );
my @speciations = qw( TNO IPCC EMEP DE94 GR95 GR05 UK98 UK08 );
my @runs = qw( Solvents_Only tagged_solvents_only );
my %data;
my @inorganic = qw( NO NO2 O3 NO3 HONO H2 HO2 HNO3 HO2NO2 H2O2 SO2 OHOP );

my %category_mapping = (
    MOZART  =>  {
        TNO     => {    BIGALK  => [ '0.130 Butanes', '0.165 Chlorinated', '0.215 Esters', '0.078 Ethers', '0.411 Higher_alkanes' ],
                        TOLUENE => [ '0.050 Other_aromatics', '0.494 Toluene', '0.026 Trimethylbenzenes', '0.430 Xylenes' ],
                    },
        IPCC    => {    BIGALK  => [ '0.341 Butanes', '0.008 Chlorinated', '0.111 Esters', '0.033 Ethers', '0.271 Higher_alkanes', '0.237 Pentanes' ],
                        TOLUENE => [ '0.244 Benzene', '0.164 Other_aromatics', '0.339 Toluene', '0.037 Trimethylbenzenes', '0.216 Xylenes' ],
                    },
        EMEP  =>    {   BIGALK  => [ '1.0 Butanes' ],
                        TOLUENE => [ '1.0 Xylenes' ],
                    },
        DE94    =>  {   BIGALK  => [ '0.277 Butanes', '0.037 Chlorinated', '0.200 Esters', '0.023 Ethers', '0.462 Higher_alkanes' ],
                        TOLUENE => [ '0.259 Other_aromatics', '0.337 Toluene', '0.134 Trimethylbenzenes', '0.270 Xylenes' ],
                    },
        GR95    =>  {   BIGALK  => [ '0.104 Chlorinated', '0.205 Esters', '0.691 Higher_alkanes' ],
                        TOLUENE => [ '0.009 Benzene', '0.558 Other_aromatics', '0.298 Toluene', '0.086 Trimethylbenzenes', '0.050 Xylenes' ],
                    },
        GR05    =>  {   BIGALK  => [ '0.092 Chlorinated', '0.352 Esters', '0.556 Higher_alkanes' ],
                        TOLUENE => [ '0.691 Other_aromatics', '0.196 Toluene', '0.066 Trimethylbenzenes', '0.047 Xylenes' ],
                    },
        UK98    =>  {   BIGALK  => [ '0.116 Butanes', '0.211 Chlorinated', '0.161 Esters', '0.089 Ethers', '0.317 Higher_alkanes', '0.107 Pentanes' ],
                        TOLUENE => [ '0.281 Other_aromatics', '0.384 Toluene', '0.093 Trimethylbenzenes', '0.242 Xylenes' ],
                    },
        UK08    =>  {   BIGALK  => [ '0.158 Butanes', '0.072 Chlorinated', '0.121 Esters', '0.099 Ethers', '0.547 Higher_alkanes' ],
                        TOLUENE => [ '0.000 Benzene', '0.336 Other_aromatics', '0.192 Toluene', '0.157 Trimethylbenzenes', '0.314 Xylenes' ],
                    },
                },
    RADM2   =>  {
        TNO     => {    HC3     => [ '0.362 Alcohols', '0.175 Butanes', '0.222 Chlorinated', '0.190 Esters', '0.007 Ethers', '0.045 Propane' ],
                        HC5     => [ '0.108 Alcohols', '0.275 Esters', '0.461 Higher_alkanes', '0.156 Ketones' ],
                        HC8     => [ '0.086 Alcohols', '0.186 Ethers', '0.728 Higher_alkanes' ],
                        KET     => [ '0.011 Alcohols', '0.989 Ketones' ],
                        TOL     => [ '0.057 Other_aromatics', '0.943 Toluene' ],
                        XYL     => [ '0.042 Other_aromatics', '0.054 Trimethylbenzenes', '0.903 Xylenes' ],
                    },
        IPCC    =>  {   HC3     => [ '0.112 Alcohols', '0.363 Butanes', '0.008 Chlorinated', '0.108 Esters', '0.002 Ethers', '0.164 Other_alkenes,_alkynes,_dienes', '0.242 Propane'] ,
                        HC5     => [ '0.032 Alcohols', '0.029 Esters', '0.232 Higher_alkanes', '0.025 Ketones', '0.681 Pentanes' ],
                        HC8     => [ '0.057 Alcohols', '0.130 Ethers', '0.813 Higher_alkanes' ],
                        KET     => [ '0.021 Alcohols', '0.979 Ketones' ],
                        TOL     => [ '0.358 Benzene', '0.144 Other_aromatics', '0.498 Toluene' ],
                        XYL     => [ '0.208 Other_aromatics', '0.116 Trimethylbenzenes', '0.676 Xylenes' ],
                    },
        EMEP    =>  {   HC3     => [ '0.419 Alcohols', '0.581 Butanes' ],
                    }, 
        DE94    =>  {   HC3     => [ '0.661 Alcohols', '0.181 Butanes', '0.086 Esters', '0.001 Ethers', '0.047 Propane' ],
                        HC5     => [ '0.317 Alcohols', '0.201 Esters', '0.402 Higher_alkanes', '0.080 Ketones' ],
                        HC8     => [ '0.539 Alcohols', '0.028 Ethers', '0.433 Higher_alkanes' ],
                        KET     => [ '0.063 Alcohols', '0.937 Ketones' ],
                        TOL     => [ '0.314 Other_aromatics', '0.686 Toluene' ],
                        XYL     => [ '0.205 Other_aromatics', '0.264 Trimethylbenzenes', '0.531 Xylenes' ],
                    },
        GR95    =>  {   HC3     => [ '0.753 Alcohols', '0.093 Chlorinated', '0.155 Esters' ],
                        HC5     => [ '0.508 Esters', '0.492 Ketones' ],
                        HC8     => [ '0.086 Alcohols', '0.914 Higher_alkanes' ],
                        OLT     => [ '1.0 Other_alkenes,_alkynes,_dienes'],
                        TOL     => [ '0.014 Benzene', '0.526 Other_aromatics', '0.460 Toluene' ],
                        XYL     => [ '0.615 Other_aromatics', '0.244 Trimethylbenzenes', '0.141 Xylenes' ],
                    },
        GR05    =>  {   HC3     => [ '0.853 Alcohols', '0.039 Chlorinated', '0.107 Esters' ],
                        HC5     => [ '0.932 Esters', '0.068 Ketones' ],
                        HC8     => [ '0.402 Alcohols', '0.598 Higher_alkanes' ],
                        KET     => [ '0.007 Alcohols', '0.993 Ketones' ],
                        OLT     => [ '1.0 Other_alkenes,_alkynes,_dienes'],
                        TOL     => [ '0.728 Other_aromatics', '0.272 Toluene' ],
                        XYL     => [ '0.596 Other_aromatics', '0.236 Trimethylbenzenes', '0.168 Xylenes' ],
                    },
        UK98    =>  {   HC3     => [ '0.492 Alcohols', '0.103 Butanes', '0.188 Chlorinated', '0.088 Esters', '0.005 Ethers', '0.124 Propane' ],
                        HC5     => [ '0.354 Alcohols', '0.154 Esters', '0.134 Higher_alkanes', '0.091 Ketones', '0.267 Pentanes' ] ,
                        HC8     => [ '0.296 Alcohols', '0.168 Ethers', '0.536 Higher_alkanes' ],
                        KET     => [ '0.059 Alcohols', '0.941 Ketones' ],
                        TOL     => [ '0.292 Other_aromatics', '0.708 Toluene' ],
                        XYL     => [ '0.268 Other_aromatics', '0.203 Trimethylbenzenes', '0.529 Xylenes' ],
                    },
        UK08    =>  {   HC3     => [ '0.763 Alcohols', '0.118 Butanes', '0.054 Chlorinated', '0.027 Esters', '0.002 Ethers', '0.000 Other_alkenes,_alkynes,_dienes', '0.036 Propane' ],
                        HC5     => [ '0.413 Alcohols', '0.220 Esters', '0.209 Higher_alkanes', '0.147 Ketones', '0.010 Pentanes' ],
                        HC8     => [ '0.283 Alcohols', '0.123 Ethers', '0.594 Higher_alkanes' ],
                        KET     => [ '0.078 Alcohols', '0.922 Ketones' ],
                        TOL     => [ '0.000 Benzene', '0.479 Other_aromatics', '0.521 Toluene' ],
                        XYL     => [ '0.253 Other_aromatics', '0.249 Trimethylbenzenes', '0.498 Xylenes' ],
                    },
                }
);

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        foreach my $run (@runs) {
            if ($run =~ "tagged" and $mechanism eq "MCM") {
                my $mcm_base = "$base/$mechanism";
                opendir DIR, $mcm_base or die "can't open $mcm_base : $!";
                my @dirs = grep { $_ =~ /${speciation}_tagged_solvents_only/ } readdir DIR;
                closedir DIR;
                my $no_runs = scalar @dirs;
                foreach my $dir (@dirs) {
                    #print $dir, "\n";
                    my $boxmodel = "$base/$mechanism/$dir/boxmodel";
                    my $mecca = MECCA->new($boxmodel);
                    my $eqn = "$base/$mechanism/$dir/gas.eqn";
                    my $kpp = KPP->new($eqn);
                    (my $label = $dir) =~ s/${speciation}_//;
                    $data{$mechanism}{$speciation}{$label}{$dir} = get_data($mecca, $kpp, $mechanism, $speciation, $run, $no_runs);
                }
            } else {
                #print $run, "\n";
                my $dir = "$base/$mechanism/${speciation}_$run";
                my $boxmodel = "$dir/boxmodel";
                my $mecca = MECCA->new($boxmodel);
                my $eqn = "$dir/gas.eqn";
                my $kpp = KPP->new($eqn);
                $data{$mechanism}{$speciation}{$run} = get_data($mecca, $kpp, $mechanism, $speciation, $run);
            }
        }
    }
}
foreach my $speciation (sort keys %{$data{"MCM"}}) {
    my %allocated;
    foreach my $run (sort keys %{$data{"MCM"}{$speciation}}) {
        next unless ($run =~ /tagged/);
        foreach my $dir (sort keys %{$data{'MCM'}{$speciation}{$run}}) {
            foreach my $item (sort keys %{$data{'MCM'}{$speciation}{$run}{$dir}}) {
                $allocated{$item} += $data{'MCM'}{$speciation}{$run}{$dir}{$item};
            }
        }
        $data{'MCM'}{$speciation}{"tagged_solvents_only"} = \%allocated;
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
    foreach my $speciation (keys %{$data{$mechanism}}) {
        $R->set('speciation', $speciation);
        #print "\t$speciation\n";
        foreach my $run (keys %{$data{$mechanism}{$speciation}}) {
            next unless ($run =~ /nly$/);
            #print "\t\t$run\n";
            $R->set('run', $run);
            $R->run(q` pre = data.frame(Mechanism = mechanism) `);
            $R->run(q` pre$Speciation = speciation `);
            $R->run(q` pre$Run = run `);
            foreach my $VOC (keys %{$data{$mechanism}{$speciation}{$run}}) {
                #print "$VOC\n";
                #print "\t\t\t$VOC : $data{$mechanism}{$speciation}{$run}{$VOC}\n";
                $R->set('voc', $VOC);
                $R->set('reactivity', $data{$mechanism}{$speciation}{$run}{$VOC});
                $R->run(q` pre[voc] = reactivity `);
            }
        $R->run(q` pre = gather(pre, VOC, Reactivity, -Run, -Speciation, -Mechanism) `,
                q` data = rbind(data, pre) `,
        );
        }
    }
}
#my $p = $R->run(q` print(data) `);
#print $p, "\n";
$R->run(q` data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08")) `);

$R->run(q` plot = ggplot(data, aes(x = Run, y = Reactivity, fill = VOC)) `,
        q` plot = plot + geom_bar(stat = "identity") `,
        q` plot = plot + facet_grid( Speciation ~ Mechanism) `,
        q` plot = plot + theme_tufte() `,
        q` plot = plot + scale_y_continuous(expand = c(0, 0)) `,
        q` plot = plot + scale_x_discrete(expand = c(0, 0)) `,
);

$R->run(q` CairoPDF(file = "Total_OH_reactivity_facet_Mechanism.pdf", width = 8.7, height = 6) `,
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
    } elsif ($process =~ /^O$|^NO3$|^H2$|^H2O2$|^SO2$|^HONO$/) {
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
    } elsif ($process =~ /Inorganic|Emissions/) {
        $category = $process;
    } else {
        $category = "Inorganic";
        #print "No category found for $process in $mechanism, $speciation\n";
    }
    return $category;
}

sub get_data {
    my ($mecca, $kpp, $mechanism, $speciation, $run, $no_runs) = @_;
    $no_runs = 1 unless (defined $no_runs);
    my %total_reactivity;
    my $cair = $mecca->cair->at(0);
    my @oxidants = qw( OH );
    foreach my $oxidant (@oxidants) {
        my $consumers = $kpp->consuming($oxidant);
        my $reactivity = 0;
        foreach my $reaction (@$consumers) {
            my $reactants = $kpp->reactants($reaction);
            my ($other_reactant) = grep { $_ ne $oxidant } @$reactants;
            next unless (defined $other_reactant);
            next if ($other_reactant eq 'hv');
            #next if ($other_reactant ~~ @inorganic);
            my $other_reactant_conc = $mecca->tracer($other_reactant) * $cair;
            my $rate_constant = get_rate_constant($reaction, $kpp, $cair);
            $rate_constant =  eval $rate_constant;
            $reactivity = $rate_constant * $other_reactant_conc;
            next if ($reactivity->sum == 0);
            $other_reactant =~ s/(.*?)_(.*?)$/$2/g;
            $other_reactant = "Inorganic" if ($other_reactant ~~ @inorganic);
            $other_reactant = "Inorganic" if ($run =~ /Solvents_Only/);
            $other_reactant = "notag" if ($kpp->reaction_string($reaction) =~ /notag/);
            $total_reactivity{$other_reactant} += $reactivity->sum;
        } 
    }
    delete $total_reactivity{"notag"};

    my %final_categories;
    foreach my $process (keys %total_reactivity) {
        if (exists $category_mapping{$mechanism}) {
            if (exists $category_mapping{$mechanism}{$speciation}) {
                if (exists $category_mapping{$mechanism}{$speciation}{$process}) {
                    foreach my $entry (@{$category_mapping{$mechanism}{$speciation}{$process}}) {
                        my ($factor, $category) = split / /, $entry;
                        $category =~ s/_/ /g;
                        $final_categories{$category} += $factor * $total_reactivity{$process};
                    }
                } else {
                    my $category = get_category($process, $mechanism, $speciation);
                    $final_categories{$category} += $total_reactivity{$process};
                }
            }
        } else {
            my $category = get_category($process, $mechanism, $speciation);
            $final_categories{$category} += $total_reactivity{$process};
        }
    }
    foreach my $process (keys %final_categories) {
        if ($process eq "Methane" or $process eq "CO" or $process eq "Inorganic") {
            $final_categories{$process} /= $no_runs;
        }
    }
    #print "$_ : $final_categories{$_}\n" for keys %final_categories;
    return \%final_categories;
}
