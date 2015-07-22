#! /usr/bin/env perl
# statistics of ozone mixing ratio differences in each mechanism between the speciations
# Version 0: Jane Coates 22/7/2015

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $species = "O3";

my $base = "/work/users/jco/Solvent_Emissions";
my @mechanisms = qw(MCM MOZART RADM2);
my @speciations = qw( DE94 EMEP GR05 GR95 IPCC TNO UK08 UK98 );
my @runs = qw( all_sectors );
my %data;

my $mecca = MECCA->new("$base/RADM2/DE94_all_sectors/boxmodel");
my $ntime = $mecca->time->nelem;
my $dt = $mecca->dt->at(0);
my $n_per_day = 86400 / $dt;
my $n_days = int $ntime / $n_per_day;
my $times = $mecca->time;
$times = $times(1:$ntime-2);

foreach my $mechanism (@mechanisms) {
    foreach my $speciation (@speciations) {
        foreach my $run (@runs) {
            my $dir = "$base/$mechanism/${speciation}_$run";
            my $mecca = MECCA->new("$dir/boxmodel");
            my $mr = $mecca->tracer($species); 
            $data{$mechanism}{$speciation} = $mr(1:$ntime-2);
        }
    }
}

my $R = Statistics::R->new();
$R->run(q` library(tidyr) `,
        q` library(dplyr) `,
        q` library(ggplot2) `,
        q` library(Cairo) `,
        q` library(ggthemes) `,
);

$R->set('Time', [ map { $_ } $times->dog ]);
$R->run(q` d = data.frame() `);
    foreach my $mechanism (sort keys %data) {
        $R->set('mechanism', $mechanism);
        $R->run(q` pre = data.frame(Mechanism = rep(mechanism, length(Time))) `);
        foreach my $speciation (sort keys %{$data{$mechanism}}) {
            $R->set('speciation', $speciation);
            $R->set('mixing.ratio', [  map { $_ } $data{$mechanism}{$speciation}->dog ]);
            $R->run(q` pre[speciation] = mixing.ratio `);
        }
#        $R->run(q` pre = gather(pre, Speciation, Mixing.Ratio, -Time, -Run, -Mechanism) `,
        $R->run(q` d = rbind(d, pre) `);
#        );
}
#$R->run(q` write.table(d, file = "ozone_mixing_ratios.csv", quote = FALSE, row.names = FALSE, sep = ",") `);
#By Mechanism analysis
$R->run(q` radm2 = d %>% filter(Mechanism == "RADM2") %>% select(-Mechanism) `,
        q` mozart = d %>% filter(Mechanism == "MOZART") %>% select(-Mechanism) `,
        q` mcm = d %>% filter(Mechanism == "MCM") %>% select(-Mechanism) `,
);

$R->run(q` radm2 = radm2 %>% rowwise() %>% mutate(Min = min(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) %>% mutate(Max = max(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) `,
        q` radm2 = radm2 %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `, 
        q` radm2.sd = sd(radm2$Largest.Diff) `,
        q` radm2.mean = mean(radm2$Largest.Diff) `,
);

$R->run(q` mozart = mozart %>% rowwise() %>% mutate(Min = min(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) %>% mutate(Max = max(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) `,
        q` mozart = mozart %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `, 
        q` mozart.sd = sd(mozart$Largest.Diff) `,
        q` mozart.mean = mean(mozart$Largest.Diff) `,
);

$R->run(q` mcm = mcm %>% rowwise() %>% mutate(Min = min(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) %>% mutate(Max = max(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) `,
        q` mcm = mcm %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `, 
        q` mcm.sd = sd(mcm$Largest.Diff) `,
        q` mcm.mean = mean(mcm$Largest.Diff) `,
);
#my $p = $R->run(q` print(radm2.day2.mean) `);
#print $p, "\n";
#$p = $R->run(q` print(radm2.day3.mean) `);
#print $p, "\n";

$R->run(q` entire.data = data.frame(Mechanism = "MCM", St.Dev = mcm.sd, Average = mcm.mean) `,
        q` entire.data = rbind(entire.data, data.frame(Mechanism = 'MOZART', St.Dev = mozart.sd, Average = mozart.mean)) `,
        q` entire.data = rbind(entire.data, data.frame(Mechanism = 'RADM2', St.Dev = radm2.sd, Average = radm2.mean)) `,
        q` entire.data$St.Dev = sprintf("%.3e", entire.data$St.Dev) `,
        q` entire.data$Average = sprintf("%.3e", entire.data$Average) `,
);
$R->run(q` write.table(entire.data, file = "all_sectors_Ozone_time_series_entire_difference_statistics_by_mechanism.txt", row.names = FALSE, quote = FALSE, sep = "\t") `);

#By Speciation analysis
$R->run(q` d$Time = Time `,
        q` d = d %>% gather(Speciation, Mixing.Ratio, -Mechanism, -Time) `,
        q` d = spread(d, Mechanism, Mixing.Ratio) `,

        q` de94 = d %>% filter(Speciation == "DE94") %>% select(-Speciation, -Time) `,
        q` de94 = de94 %>% rowwise() %>% mutate(Min = min(MCM, MOZART, RADM2)) %>% mutate(Max = max(MCM, MOZART, RADM2)) `,
        q` de94 = de94 %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `,
        q` de94.entire.sd = sd(de94$Largest.Diff) `,
        q` de94.entire.mean = mean(de94$Largest.Diff) `,

        q` emep = d %>% filter(Speciation == "EMEP") %>% select(-Speciation, -Time) `,
        q` emep = emep %>% rowwise() %>% mutate(Min = min(MCM, MOZART, RADM2)) %>% mutate(Max = max(MCM, MOZART, RADM2)) `,
        q` emep = emep %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `,
        q` emep.entire.sd = sd(emep$Largest.Diff) `,
        q` emep.entire.mean = mean(emep$Largest.Diff) `,

        q` gr05 = d %>% filter(Speciation == "GR05") %>% select(-Speciation, -Time) `,
        q` gr05 = gr05 %>% rowwise() %>% mutate(Min = min(MCM, MOZART, RADM2)) %>% mutate(Max = max(MCM, MOZART, RADM2)) `,
        q` gr05 = gr05 %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `,
        q` gr05.entire.sd = sd(gr05$Largest.Diff) `,
        q` gr05.entire.mean = mean(gr05$Largest.Diff) `,

        q` gr95 = d %>% filter(Speciation == "GR95") %>% select(-Speciation, -Time) `,
        q` gr95 = gr95 %>% rowwise() %>% mutate(Min = min(MCM, MOZART, RADM2)) %>% mutate(Max = max(MCM, MOZART, RADM2)) `,
        q` gr95 = gr95 %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `,
        q` gr95.entire.sd = sd(gr95$Largest.Diff) `,
        q` gr95.entire.mean = mean(gr95$Largest.Diff) `,

        q` ipcc = d %>% filter(Speciation == "IPCC") %>% select(-Speciation, -Time) `,
        q` ipcc = ipcc %>% rowwise() %>% mutate(Min = min(MCM, MOZART, RADM2)) %>% mutate(Max = max(MCM, MOZART, RADM2)) `,
        q` ipcc = ipcc %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `,
        q` ipcc.entire.sd = sd(ipcc$Largest.Diff) `,
        q` ipcc.entire.mean = mean(ipcc$Largest.Diff) `,

        q` tno = d %>% filter(Speciation == "TNO") %>% select(-Speciation, -Time) `,
        q` tno = tno %>% rowwise() %>% mutate(Min = min(MCM, MOZART, RADM2)) %>% mutate(Max = max(MCM, MOZART, RADM2)) `,
        q` tno = tno %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `,
        q` tno.entire.sd = sd(tno$Largest.Diff) `,
        q` tno.entire.mean = mean(tno$Largest.Diff) `,

        q` uk98 = d %>% filter(Speciation == "UK98") %>% select(-Speciation, -Time) `,
        q` uk98 = uk98 %>% rowwise() %>% mutate(Min = min(MCM, MOZART, RADM2)) %>% mutate(Max = max(MCM, MOZART, RADM2)) `,
        q` uk98 = uk98 %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `,
        q` uk98.entire.sd = sd(uk98$Largest.Diff) `,
        q` uk98.entire.mean = mean(uk98$Largest.Diff) `,

        q` uk08 = d %>% filter(Speciation == "UK08") %>% select(-Speciation, -Time) `,
        q` uk08 = uk08 %>% rowwise() %>% mutate(Min = min(MCM, MOZART, RADM2)) %>% mutate(Max = max(MCM, MOZART, RADM2)) `,
        q` uk08 = uk08 %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `,
        q` uk08.entire.sd = sd(uk08$Largest.Diff) `,
        q` uk08.entire.mean = mean(uk08$Largest.Diff) `,
);
#my $p = $R->run(q` print(length(de94$Max)) `);
#print $p, "\n";

$R->run(q` entire.data = data.frame(Speciation = "DE94", St.Dev = de94.entire.sd, Average = de94.entire.mean) `,
        q` entire.data = rbind(entire.data, data.frame(Speciation = 'EMEP', St.Dev = emep.entire.sd, Average = emep.entire.mean)) `,
        q` entire.data = rbind(entire.data, data.frame(Speciation = 'GR05', St.Dev = gr05.entire.sd, Average = gr05.entire.mean)) `,
        q` entire.data = rbind(entire.data, data.frame(Speciation = 'GR95', St.Dev = gr95.entire.sd, Average = gr95.entire.mean)) `,
        q` entire.data = rbind(entire.data, data.frame(Speciation = 'IPCC', St.Dev = ipcc.entire.sd, Average = ipcc.entire.mean)) `,
        q` entire.data = rbind(entire.data, data.frame(Speciation = 'TNO', St.Dev = tno.entire.sd, Average = tno.entire.mean)) `,
        q` entire.data = rbind(entire.data, data.frame(Speciation = 'UK08', St.Dev = uk08.entire.sd, Average = uk08.entire.mean)) `,
        q` entire.data = rbind(entire.data, data.frame(Speciation = 'UK98', St.Dev = uk98.entire.sd, Average = uk98.entire.mean)) `,
        q` entire.data$St.Dev = sprintf("%.3e", entire.data$St.Dev) `,
        q` entire.data$Average = sprintf("%.3e", entire.data$Average) `,
);
$R->run(q` write.table(entire.data, file = "all_sectors_Ozone_time_series_entire_difference_statistics_by_speciation.txt", row.names = FALSE, quote = FALSE, sep = "\t") `);

$R->stop();
