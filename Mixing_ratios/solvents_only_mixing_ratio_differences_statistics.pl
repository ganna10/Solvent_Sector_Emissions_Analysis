#! /usr/bin/env perl
# statistics of ozone mixing ratio differences in each mechanism between the speciations
# Version 0: Jane Coates 25/6/2015

use strict;
use diagnostics;
use MECCA;
use PDL;
use PDL::NiceSlice;
use Statistics::R;

my $species = "O3";

my $base = "/local/home/coates/Solvent_Emissions";
my @mechanisms = qw(MCM MOZART RADM2);
my @speciations = qw( DE94 EMEP GR05 GR95 IPCC TNO UK08 UK98 );
my @runs = qw( Solvents_Only );
my %data;

my $mecca = MECCA->new("$base/RADM2/DE94_Solvents_Only/boxmodel");
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
$R->run(q` write.table(d, file = "solvents_only_ozone_mixing_ratios.csv", quote = FALSE, row.names = FALSE, sep = ",") `);
#By Mechanism analysis
#$R->run(q` radm2 = d %>% filter(Mechanism == "RADM2") %>% select(-Mechanism) `,
#        q` mozart = d %>% filter(Mechanism == "MOZART") %>% select(-Mechanism) `,
#        q` mcm = d %>% filter(Mechanism == "MCM") %>% select(-Mechanism) `,
#);
#
#$R->run(q` radm2 = radm2 %>% rowwise() %>% mutate(Min = min(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) %>% mutate(Max = max(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) `,
#        q` radm2 = radm2 %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `, 
#        q` radm2.sd = sd(radm2$Largest.Diff) `,
#        q` radm2.mean = mean(radm2$Largest.Diff) `,
#
#        q` radm2.day1 = radm2[1:72,] `,
#        q` radm2.day1.sd = sd(radm2.day1$Largest.Diff) `,
#        q` radm2.day1.mean = mean(radm2.day1$Largest.Diff) `,
#
#        q` radm2.day2 = radm2[73:144,] `,
#        q` radm2.day2.sd = sd(radm2.day2$Largest.Diff) `,
#        q` radm2.day2.mean = mean(radm2.day2$Largest.Diff) `,
#
#        q` radm2.day3 = radm2[145:216,] `,
#        q` radm2.day3.sd = sd(radm2.day3$Largest.Diff) `,
#        q` radm2.day3.mean = mean(radm2.day3$Largest.Diff) `,
#
#        q` radm2.day4 = radm2[217:288,] `,
#        q` radm2.day4.sd = sd(radm2.day4$Largest.Diff) `,
#        q` radm2.day4.mean = mean(radm2.day4$Largest.Diff) `,
#
#        q` radm2.day5 = radm2[289:360,] `,
#        q` radm2.day5.sd = sd(radm2.day5$Largest.Diff) `,
#        q` radm2.day5.mean = mean(radm2.day5$Largest.Diff) `,
#
#        q` radm2.day6 = radm2[361:432,] `,
#        q` radm2.day6.sd = sd(radm2.day6$Largest.Diff) `,
#        q` radm2.day6.mean = mean(radm2.day6$Largest.Diff) `,
#
#        q` radm2.day7 = radm2[433:504,] `,
#        q` radm2.day7.sd = sd(radm2.day7$Largest.Diff) `,
#        q` radm2.day7.mean = mean(radm2.day7$Largest.Diff) `,
#);
#
#$R->run(q` mozart = mozart %>% rowwise() %>% mutate(Min = min(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) %>% mutate(Max = max(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) `,
#        q` mozart = mozart %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `, 
#        q` mozart.sd = sd(mozart$Largest.Diff) `,
#        q` mozart.mean = mean(mozart$Largest.Diff) `,
#
#        q` mozart.day1 = mozart[1:72,] `,
#        q` mozart.day1.sd = sd(mozart.day1$Largest.Diff) `,
#        q` mozart.day1.mean = mean(mozart.day1$Largest.Diff) `,
#
#        q` mozart.day2 = mozart[73:144,] `,
#        q` mozart.day2.sd = sd(mozart.day2$Largest.Diff) `,
#        q` mozart.day2.mean = mean(mozart.day2$Largest.Diff) `,
#
#        q` mozart.day3 = mozart[145:216,] `,
#        q` mozart.day3.sd = sd(mozart.day3$Largest.Diff) `,
#        q` mozart.day3.mean = mean(mozart.day3$Largest.Diff) `,
#
#        q` mozart.day4 = mozart[217:288,] `,
#        q` mozart.day4.sd = sd(mozart.day4$Largest.Diff) `,
#        q` mozart.day4.mean = mean(mozart.day4$Largest.Diff) `,
#
#        q` mozart.day5 = mozart[289:360,] `,
#        q` mozart.day5.sd = sd(mozart.day5$Largest.Diff) `,
#        q` mozart.day5.mean = mean(mozart.day5$Largest.Diff) `,
#
#        q` mozart.day6 = mozart[361:432,] `,
#        q` mozart.day6.sd = sd(mozart.day6$Largest.Diff) `,
#        q` mozart.day6.mean = mean(mozart.day6$Largest.Diff) `,
#
#        q` mozart.day7 = mozart[433:504,] `,
#        q` mozart.day7.sd = sd(mozart.day7$Largest.Diff) `,
#        q` mozart.day7.mean = mean(mozart.day7$Largest.Diff) `,
#);
#
#$R->run(q` mcm = mcm %>% rowwise() %>% mutate(Min = min(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) %>% mutate(Max = max(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) `,
#        q` mcm = mcm %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `, 
#        q` mcm.sd = sd(mcm$Largest.Diff) `,
#        q` mcm.mean = mean(mcm$Largest.Diff) `,
#
#        q` mcm.day1 = mcm[1:72,] `,
#        q` mcm.day1.sd = sd(mcm.day1$Largest.Diff) `,
#        q` mcm.day1.mean = mean(mcm.day1$Largest.Diff) `,
#
#        q` mcm.day2 = mcm[73:144,] `,
#        q` mcm.day2.sd = sd(mcm.day2$Largest.Diff) `,
#        q` mcm.day2.mean = mean(mcm.day2$Largest.Diff) `,
#
#        q` mcm.day3 = mcm[145:216,] `,
#        q` mcm.day3.sd = sd(mcm.day3$Largest.Diff) `,
#        q` mcm.day3.mean = mean(mcm.day3$Largest.Diff) `,
#
#        q` mcm.day4 = mcm[217:288,] `,
#        q` mcm.day4.sd = sd(mcm.day4$Largest.Diff) `,
#        q` mcm.day4.mean = mean(mcm.day4$Largest.Diff) `,
#
#        q` mcm.day5 = mcm[289:360,] `,
#        q` mcm.day5.sd = sd(mcm.day5$Largest.Diff) `,
#        q` mcm.day5.mean = mean(mcm.day5$Largest.Diff) `,
#
#        q` mcm.day6 = mcm[361:432,] `,
#        q` mcm.day6.sd = sd(mcm.day6$Largest.Diff) `,
#        q` mcm.day6.mean = mean(mcm.day6$Largest.Diff) `,
#
#        q` mcm.day7 = mcm[433:504,] `,
#        q` mcm.day7.sd = sd(mcm.day7$Largest.Diff) `,
#        q` mcm.day7.mean = mean(mcm.day7$Largest.Diff) `,
#);
##my $p = $R->run(q` print(radm2.day2.mean) `);
##print $p, "\n";
##$p = $R->run(q` print(radm2.day3.mean) `);
##print $p, "\n";
#
#$R->run(q` entire.data = data.frame(Mechanism = "MCM", St.Dev = mcm.sd, Average = mcm.mean) `,
#        q` entire.data = rbind(entire.data, data.frame(Mechanism = 'MOZART', St.Dev = mozart.sd, Average = mozart.mean)) `,
#        q` entire.data = rbind(entire.data, data.frame(Mechanism = 'RADM2', St.Dev = radm2.sd, Average = radm2.mean)) `,
#        q` entire.data$St.Dev = sprintf("%.3e", entire.data$St.Dev) `,
#        q` entire.data$Average = sprintf("%.3e", entire.data$Average) `,
#);
#$R->run(q` write.table(entire.data, file = "Solvents_Only_Ozone_time_series_entire_difference_statistics_by_mechanism.txt", row.names = FALSE, quote = FALSE, sep = "\t") `);
#
#$R->set('Days', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
#$R->run(q` radm2.data = data.frame(Mechanism = rep("RADM2", 7), Days) `,
#        q` radm2.data$Std.Dev = c(radm2.day1.sd, radm2.day2.sd, radm2.day3.sd, radm2.day4.sd, radm2.day5.sd, radm2.day6.sd, radm2.day7.sd) `, 
#        q` radm2.data$Average = c(radm2.day1.mean, radm2.day2.mean, radm2.day3.mean, radm2.day4.mean, radm2.day5.mean, radm2.day6.mean, radm2.day7.mean) `, 
#);
#
#$R->run(q` mozart.data = data.frame(Mechanism = rep("MOZART", 7), Days) `,
#        q` mozart.data$Std.Dev = c(mozart.day1.sd, mozart.day2.sd, mozart.day3.sd, mozart.day4.sd, mozart.day5.sd, mozart.day6.sd, mozart.day7.sd) `, 
#        q` mozart.data$Average = c(mozart.day1.mean, mozart.day2.mean, mozart.day3.mean, mozart.day4.mean, mozart.day5.mean, mozart.day6.mean, mozart.day7.mean) `, 
#);
#
#$R->run(q` mcm.data = data.frame(Mechanism = rep("MCM", 7), Days) `,
#        q` mcm.data$Std.Dev = c(mcm.day1.sd, mcm.day2.sd, mcm.day3.sd, mcm.day4.sd, mcm.day5.sd, mcm.day6.sd, mcm.day7.sd) `, 
#        q` mcm.data$Average = c(mcm.day1.mean, mcm.day2.mean, mcm.day3.mean, mcm.day4.mean, mcm.day5.mean, mcm.day6.mean, mcm.day7.mean) `, 
#);
#
#$R->run(q` print.data = rbind(mcm.data, mozart.data, radm2.data) `,
#        q` print.data$Std.Dev = sprintf("%.3e", print.data$Std.Dev) `,
#        q` print.data$Average = sprintf("%.3e", print.data$Average) `,
#);
#$R->run(q` write.table(print.data, file = "Solvents_Only_Ozone_time_series_daily_difference_statistics_by_mechanism.txt", row.names = FALSE, quote = FALSE, sep = "\t") `);
#
#$R->run(q` plot.data = print.data %>% gather(Variable, Value, -Mechanism, -Days) `,
#        q` plot.data$Value = as.numeric(as.character(plot.data$Value)) `,
#);
##my $p = $R->run(q` print(plot.data) `);
##print $p, "\n";
#$R->run(q` p = ggplot(plot.data, aes( x = Days, y = Value, colour = Mechanism)) `,
#        q` p = p + geom_point() `,
#        q` p = p + facet_wrap( ~ Variable, scales = "free_y") `,
#        q` p = p + theme_hc() `,
#        q` p = p + ggtitle("Solvents Only - by Mechanism") `,
#        q` p = p + theme(axis.line = element_line(colour = "black")) `,
#        q` p = p + theme(strip.text = element_text(face = "bold")) `,
#        q` p = p + theme(strip.background = element_blank()) `,
#);
#
#$R->run(q` CairoPDF(file = "Solvents_only_Ozone_statistics_plots_by_mechanism.pdf", width = 11, height = 7) `,
#        q` print(p) `,
#        q` dev.off() `,
#);

#By Speciation analysis
#$R->run(q` d$Time = Time `,
#        q` d = d %>% gather(Speciation, Mixing.Ratio, -Mechanism, -Time) `,
#        q` d = spread(d, Mechanism, Mixing.Ratio) `,
#
#        q` de94 = d %>% filter(Speciation == "DE94") %>% select(-Speciation, -Time) `,
#        q` de94 = de94 %>% rowwise() %>% mutate(Min = min(MCM, MOZART, RADM2)) %>% mutate(Max = max(MCM, MOZART, RADM2)) `,
#        q` de94 = de94 %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `,
#        q` de94.entire.sd = sd(de94$Largest.Diff) `,
#        q` de94.entire.mean = mean(de94$Largest.Diff) `,
#
#        q` de94.day1 = de94[1:72,] `,
#        q` de94.day1.sd = sd(de94.day1$Largest.Diff) `,
#        q` de94.day1.mean = mean(de94.day1$Largest.Diff) `,
#
#        q` de94.day2 = de94[73:144,] `,
#        q` de94.day2.sd = sd(de94.day2$Largest.Diff) `,
#        q` de94.day2.mean = mean(de94.day2$Largest.Diff) `,
#
#        q` de94.day3 = de94[145:216,] `,
#        q` de94.day3.sd = sd(de94.day3$Largest.Diff) `,
#        q` de94.day3.mean = mean(de94.day3$Largest.Diff) `,
#
#        q` de94.day4 = de94[217:288,] `,
#        q` de94.day4.sd = sd(de94.day4$Largest.Diff) `,
#        q` de94.day4.mean = mean(de94.day4$Largest.Diff) `,
#
#        q` de94.day5 = de94[289:360,] `,
#        q` de94.day5.sd = sd(de94.day5$Largest.Diff) `,
#        q` de94.day5.mean = mean(de94.day5$Largest.Diff) `,
#
#        q` de94.day6 = de94[361:432,] `,
#        q` de94.day6.sd = sd(de94.day6$Largest.Diff) `,
#        q` de94.day6.mean = mean(de94.day6$Largest.Diff) `,
#
#        q` de94.day7 = de94[433:504,] `,
#        q` de94.day7.sd = sd(de94.day7$Largest.Diff) `,
#        q` de94.day7.mean = mean(de94.day7$Largest.Diff) `, 
#
#        q` emep = d %>% filter(Speciation == "EMEP") %>% select(-Speciation, -Time) `,
#        q` emep = emep %>% rowwise() %>% mutate(Min = min(MCM, MOZART, RADM2)) %>% mutate(Max = max(MCM, MOZART, RADM2)) `,
#        q` emep = emep %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `,
#        q` emep.entire.sd = sd(emep$Largest.Diff) `,
#        q` emep.entire.mean = mean(emep$Largest.Diff) `,
#
#        q` emep.day1 = emep[1:72,] `,
#        q` emep.day1.sd = sd(emep.day1$Largest.Diff) `,
#        q` emep.day1.mean = mean(emep.day1$Largest.Diff) `,
#
#        q` emep.day2 = emep[73:144,] `,
#        q` emep.day2.sd = sd(emep.day2$Largest.Diff) `,
#        q` emep.day2.mean = mean(emep.day2$Largest.Diff) `,
#
#        q` emep.day3 = emep[145:216,] `,
#        q` emep.day3.sd = sd(emep.day3$Largest.Diff) `,
#        q` emep.day3.mean = mean(emep.day3$Largest.Diff) `,
#
#        q` emep.day4 = emep[217:288,] `,
#        q` emep.day4.sd = sd(emep.day4$Largest.Diff) `,
#        q` emep.day4.mean = mean(emep.day4$Largest.Diff) `,
#
#        q` emep.day5 = emep[289:360,] `,
#        q` emep.day5.sd = sd(emep.day5$Largest.Diff) `,
#        q` emep.day5.mean = mean(emep.day5$Largest.Diff) `,
#
#        q` emep.day6 = emep[361:432,] `,
#        q` emep.day6.sd = sd(emep.day6$Largest.Diff) `,
#        q` emep.day6.mean = mean(emep.day6$Largest.Diff) `,
#
#        q` emep.day7 = emep[433:504,] `,
#        q` emep.day7.sd = sd(emep.day7$Largest.Diff) `,
#        q` emep.day7.mean = mean(emep.day7$Largest.Diff) `, 
#
#        q` gr05 = d %>% filter(Speciation == "GR05") %>% select(-Speciation, -Time) `,
#        q` gr05 = gr05 %>% rowwise() %>% mutate(Min = min(MCM, MOZART, RADM2)) %>% mutate(Max = max(MCM, MOZART, RADM2)) `,
#        q` gr05 = gr05 %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `,
#        q` gr05.entire.sd = sd(gr05$Largest.Diff) `,
#        q` gr05.entire.mean = mean(gr05$Largest.Diff) `,
#
#        q` gr05.day1 = gr05[1:72,] `,
#        q` gr05.day1.sd = sd(gr05.day1$Largest.Diff) `,
#        q` gr05.day1.mean = mean(gr05.day1$Largest.Diff) `,
#
#        q` gr05.day2 = gr05[73:144,] `,
#        q` gr05.day2.sd = sd(gr05.day2$Largest.Diff) `,
#        q` gr05.day2.mean = mean(gr05.day2$Largest.Diff) `,
#
#        q` gr05.day3 = gr05[145:216,] `,
#        q` gr05.day3.sd = sd(gr05.day3$Largest.Diff) `,
#        q` gr05.day3.mean = mean(gr05.day3$Largest.Diff) `,
#
#        q` gr05.day4 = gr05[217:288,] `,
#        q` gr05.day4.sd = sd(gr05.day4$Largest.Diff) `,
#        q` gr05.day4.mean = mean(gr05.day4$Largest.Diff) `,
#
#        q` gr05.day5 = gr05[289:360,] `,
#        q` gr05.day5.sd = sd(gr05.day5$Largest.Diff) `,
#        q` gr05.day5.mean = mean(gr05.day5$Largest.Diff) `,
#
#        q` gr05.day6 = gr05[361:432,] `,
#        q` gr05.day6.sd = sd(gr05.day6$Largest.Diff) `,
#        q` gr05.day6.mean = mean(gr05.day6$Largest.Diff) `,
#
#        q` gr05.day7 = gr05[433:504,] `,
#        q` gr05.day7.sd = sd(gr05.day7$Largest.Diff) `,
#        q` gr05.day7.mean = mean(gr05.day7$Largest.Diff) `, 
#
#        q` gr95 = d %>% filter(Speciation == "GR95") %>% select(-Speciation, -Time) `,
#        q` gr95 = gr95 %>% rowwise() %>% mutate(Min = min(MCM, MOZART, RADM2)) %>% mutate(Max = max(MCM, MOZART, RADM2)) `,
#        q` gr95 = gr95 %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `,
#        q` gr95.entire.sd = sd(gr95$Largest.Diff) `,
#        q` gr95.entire.mean = mean(gr95$Largest.Diff) `,
#
#        q` gr95.day1 = gr95[1:72,] `,
#        q` gr95.day1.sd = sd(gr95.day1$Largest.Diff) `,
#        q` gr95.day1.mean = mean(gr95.day1$Largest.Diff) `,
#
#        q` gr95.day2 = gr95[73:144,] `,
#        q` gr95.day2.sd = sd(gr95.day2$Largest.Diff) `,
#        q` gr95.day2.mean = mean(gr95.day2$Largest.Diff) `,
#
#        q` gr95.day3 = gr95[145:216,] `,
#        q` gr95.day3.sd = sd(gr95.day3$Largest.Diff) `,
#        q` gr95.day3.mean = mean(gr95.day3$Largest.Diff) `,
#
#        q` gr95.day4 = gr95[217:288,] `,
#        q` gr95.day4.sd = sd(gr95.day4$Largest.Diff) `,
#        q` gr95.day4.mean = mean(gr95.day4$Largest.Diff) `,
#
#        q` gr95.day5 = gr95[289:360,] `,
#        q` gr95.day5.sd = sd(gr95.day5$Largest.Diff) `,
#        q` gr95.day5.mean = mean(gr95.day5$Largest.Diff) `,
#
#        q` gr95.day6 = gr95[361:432,] `,
#        q` gr95.day6.sd = sd(gr95.day6$Largest.Diff) `,
#        q` gr95.day6.mean = mean(gr95.day6$Largest.Diff) `,
#
#        q` gr95.day7 = gr95[433:504,] `,
#        q` gr95.day7.sd = sd(gr95.day7$Largest.Diff) `,
#        q` gr95.day7.mean = mean(gr95.day7$Largest.Diff) `, 
#
#        q` ipcc = d %>% filter(Speciation == "IPCC") %>% select(-Speciation, -Time) `,
#        q` ipcc = ipcc %>% rowwise() %>% mutate(Min = min(MCM, MOZART, RADM2)) %>% mutate(Max = max(MCM, MOZART, RADM2)) `,
#        q` ipcc = ipcc %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `,
#        q` ipcc.entire.sd = sd(ipcc$Largest.Diff) `,
#        q` ipcc.entire.mean = mean(ipcc$Largest.Diff) `,
#
#        q` ipcc.day1 = ipcc[1:72,] `,
#        q` ipcc.day1.sd = sd(ipcc.day1$Largest.Diff) `,
#        q` ipcc.day1.mean = mean(ipcc.day1$Largest.Diff) `,
#
#        q` ipcc.day2 = ipcc[73:144,] `,
#        q` ipcc.day2.sd = sd(ipcc.day2$Largest.Diff) `,
#        q` ipcc.day2.mean = mean(ipcc.day2$Largest.Diff) `,
#
#        q` ipcc.day3 = ipcc[145:216,] `,
#        q` ipcc.day3.sd = sd(ipcc.day3$Largest.Diff) `,
#        q` ipcc.day3.mean = mean(ipcc.day3$Largest.Diff) `,
#
#        q` ipcc.day4 = ipcc[217:288,] `,
#        q` ipcc.day4.sd = sd(ipcc.day4$Largest.Diff) `,
#        q` ipcc.day4.mean = mean(ipcc.day4$Largest.Diff) `,
#
#        q` ipcc.day5 = ipcc[289:360,] `,
#        q` ipcc.day5.sd = sd(ipcc.day5$Largest.Diff) `,
#        q` ipcc.day5.mean = mean(ipcc.day5$Largest.Diff) `,
#
#        q` ipcc.day6 = ipcc[361:432,] `,
#        q` ipcc.day6.sd = sd(ipcc.day6$Largest.Diff) `,
#        q` ipcc.day6.mean = mean(ipcc.day6$Largest.Diff) `,
#
#        q` ipcc.day7 = ipcc[433:504,] `,
#        q` ipcc.day7.sd = sd(ipcc.day7$Largest.Diff) `,
#        q` ipcc.day7.mean = mean(ipcc.day7$Largest.Diff) `, 
#
#        q` tno = d %>% filter(Speciation == "TNO") %>% select(-Speciation, -Time) `,
#        q` tno = tno %>% rowwise() %>% mutate(Min = min(MCM, MOZART, RADM2)) %>% mutate(Max = max(MCM, MOZART, RADM2)) `,
#        q` tno = tno %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `,
#        q` tno.entire.sd = sd(tno$Largest.Diff) `,
#        q` tno.entire.mean = mean(tno$Largest.Diff) `,
#
#        q` tno.day1 = tno[1:72,] `,
#        q` tno.day1.sd = sd(tno.day1$Largest.Diff) `,
#        q` tno.day1.mean = mean(tno.day1$Largest.Diff) `,
#
#        q` tno.day2 = tno[73:144,] `,
#        q` tno.day2.sd = sd(tno.day2$Largest.Diff) `,
#        q` tno.day2.mean = mean(tno.day2$Largest.Diff) `,
#
#        q` tno.day3 = tno[145:216,] `,
#        q` tno.day3.sd = sd(tno.day3$Largest.Diff) `,
#        q` tno.day3.mean = mean(tno.day3$Largest.Diff) `,
#
#        q` tno.day4 = tno[217:288,] `,
#        q` tno.day4.sd = sd(tno.day4$Largest.Diff) `,
#        q` tno.day4.mean = mean(tno.day4$Largest.Diff) `,
#
#        q` tno.day5 = tno[289:360,] `,
#        q` tno.day5.sd = sd(tno.day5$Largest.Diff) `,
#        q` tno.day5.mean = mean(tno.day5$Largest.Diff) `,
#
#        q` tno.day6 = tno[361:432,] `,
#        q` tno.day6.sd = sd(tno.day6$Largest.Diff) `,
#        q` tno.day6.mean = mean(tno.day6$Largest.Diff) `,
#
#        q` tno.day7 = tno[433:504,] `,
#        q` tno.day7.sd = sd(tno.day7$Largest.Diff) `,
#        q` tno.day7.mean = mean(tno.day7$Largest.Diff) `, 
#
#        q` uk98 = d %>% filter(Speciation == "UK98") %>% select(-Speciation, -Time) `,
#        q` uk98 = uk98 %>% rowwise() %>% mutate(Min = min(MCM, MOZART, RADM2)) %>% mutate(Max = max(MCM, MOZART, RADM2)) `,
#        q` uk98 = uk98 %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `,
#        q` uk98.entire.sd = sd(uk98$Largest.Diff) `,
#        q` uk98.entire.mean = mean(uk98$Largest.Diff) `,
#
#        q` uk98.day1 = uk98[1:72,] `,
#        q` uk98.day1.sd = sd(uk98.day1$Largest.Diff) `,
#        q` uk98.day1.mean = mean(uk98.day1$Largest.Diff) `,
#
#        q` uk98.day2 = uk98[73:144,] `,
#        q` uk98.day2.sd = sd(uk98.day2$Largest.Diff) `,
#        q` uk98.day2.mean = mean(uk98.day2$Largest.Diff) `,
#
#        q` uk98.day3 = uk98[145:216,] `,
#        q` uk98.day3.sd = sd(uk98.day3$Largest.Diff) `,
#        q` uk98.day3.mean = mean(uk98.day3$Largest.Diff) `,
#
#        q` uk98.day4 = uk98[217:288,] `,
#        q` uk98.day4.sd = sd(uk98.day4$Largest.Diff) `,
#        q` uk98.day4.mean = mean(uk98.day4$Largest.Diff) `,
#
#        q` uk98.day5 = uk98[289:360,] `,
#        q` uk98.day5.sd = sd(uk98.day5$Largest.Diff) `,
#        q` uk98.day5.mean = mean(uk98.day5$Largest.Diff) `,
#
#        q` uk98.day6 = uk98[361:432,] `,
#        q` uk98.day6.sd = sd(uk98.day6$Largest.Diff) `,
#        q` uk98.day6.mean = mean(uk98.day6$Largest.Diff) `,
#
#        q` uk98.day7 = uk98[433:504,] `,
#        q` uk98.day7.sd = sd(uk98.day7$Largest.Diff) `,
#        q` uk98.day7.mean = mean(uk98.day7$Largest.Diff) `, 
#
#        q` uk08 = d %>% filter(Speciation == "UK08") %>% select(-Speciation, -Time) `,
#        q` uk08 = uk08 %>% rowwise() %>% mutate(Min = min(MCM, MOZART, RADM2)) %>% mutate(Max = max(MCM, MOZART, RADM2)) `,
#        q` uk08 = uk08 %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `,
#        q` uk08.entire.sd = sd(uk08$Largest.Diff) `,
#        q` uk08.entire.mean = mean(uk08$Largest.Diff) `,
#
#        q` uk08.day1 = uk08[1:72,] `,
#        q` uk08.day1.sd = sd(uk08.day1$Largest.Diff) `,
#        q` uk08.day1.mean = mean(uk08.day1$Largest.Diff) `,
#
#        q` uk08.day2 = uk08[73:144,] `,
#        q` uk08.day2.sd = sd(uk08.day2$Largest.Diff) `,
#        q` uk08.day2.mean = mean(uk08.day2$Largest.Diff) `,
#
#        q` uk08.day3 = uk08[145:216,] `,
#        q` uk08.day3.sd = sd(uk08.day3$Largest.Diff) `,
#        q` uk08.day3.mean = mean(uk08.day3$Largest.Diff) `,
#
#        q` uk08.day4 = uk08[217:288,] `,
#        q` uk08.day4.sd = sd(uk08.day4$Largest.Diff) `,
#        q` uk08.day4.mean = mean(uk08.day4$Largest.Diff) `,
#
#        q` uk08.day5 = uk08[289:360,] `,
#        q` uk08.day5.sd = sd(uk08.day5$Largest.Diff) `,
#        q` uk08.day5.mean = mean(uk08.day5$Largest.Diff) `,
#
#        q` uk08.day6 = uk08[361:432,] `,
#        q` uk08.day6.sd = sd(uk08.day6$Largest.Diff) `,
#        q` uk08.day6.mean = mean(uk08.day6$Largest.Diff) `,
#
#        q` uk08.day7 = uk08[433:504,] `,
#        q` uk08.day7.sd = sd(uk08.day7$Largest.Diff) `,
#        q` uk08.day7.mean = mean(uk08.day7$Largest.Diff) `, 
#);
#my $p = $R->run(q` print(length(de94$Max)) `);
#print $p, "\n";

#$R->run(q` entire.data = data.frame(Speciation = "DE94", St.Dev = de94.entire.sd, Average = de94.entire.mean) `,
#        q` entire.data = rbind(entire.data, data.frame(Speciation = 'EMEP', St.Dev = emep.entire.sd, Average = emep.entire.mean)) `,
#        q` entire.data = rbind(entire.data, data.frame(Speciation = 'GR05', St.Dev = gr05.entire.sd, Average = gr05.entire.mean)) `,
#        q` entire.data = rbind(entire.data, data.frame(Speciation = 'GR95', St.Dev = gr95.entire.sd, Average = gr95.entire.mean)) `,
#        q` entire.data = rbind(entire.data, data.frame(Speciation = 'IPCC', St.Dev = ipcc.entire.sd, Average = ipcc.entire.mean)) `,
#        q` entire.data = rbind(entire.data, data.frame(Speciation = 'TNO', St.Dev = tno.entire.sd, Average = tno.entire.mean)) `,
#        q` entire.data = rbind(entire.data, data.frame(Speciation = 'UK08', St.Dev = uk08.entire.sd, Average = uk08.entire.mean)) `,
#        q` entire.data = rbind(entire.data, data.frame(Speciation = 'UK98', St.Dev = uk98.entire.sd, Average = uk98.entire.mean)) `,
#        q` entire.data$St.Dev = sprintf("%.3e", entire.data$St.Dev) `,
#        q` entire.data$Average = sprintf("%.3e", entire.data$Average) `,
#);
#$R->run(q` write.table(entire.data, file = "Solvents_Only_Ozone_time_series_entire_difference_statistics_by_speciation.txt", row.names = FALSE, quote = FALSE, sep = "\t") `);
#
#$R->run(q` de94.data = data.frame(Speciation = rep("DE94", 7), Days) `,
#        q` de94.data$Std.Dev = c(de94.day1.sd, de94.day2.sd, de94.day3.sd, de94.day4.sd, de94.day5.sd, de94.day6.sd, de94.day7.sd) `, 
#        q` de94.data$Average = c(de94.day1.mean, de94.day2.mean, de94.day3.mean, de94.day4.mean, de94.day5.mean, de94.day6.mean, de94.day7.mean) `, 
#);
#
#$R->run(q` emep.data = data.frame(Speciation = rep("EMEP", 7), Days) `,
#        q` emep.data$Std.Dev = c(emep.day1.sd, emep.day2.sd, emep.day3.sd, emep.day4.sd, emep.day5.sd, emep.day6.sd, emep.day7.sd) `, 
#        q` emep.data$Average = c(emep.day1.mean, emep.day2.mean, emep.day3.mean, emep.day4.mean, emep.day5.mean, emep.day6.mean, emep.day7.mean) `, 
#);
#
#$R->run(q` gr95.data = data.frame(Speciation = rep("GR95", 7), Days) `,
#        q` gr95.data$Std.Dev = c(gr95.day1.sd, gr95.day2.sd, gr95.day3.sd, gr95.day4.sd, gr95.day5.sd, gr95.day6.sd, gr95.day7.sd) `, 
#        q` gr95.data$Average = c(gr95.day1.mean, gr95.day2.mean, gr95.day3.mean, gr95.day4.mean, gr95.day5.mean, gr95.day6.mean, gr95.day7.mean) `, 
#);
#
#$R->run(q` gr05.data = data.frame(Speciation = rep("GR05", 7), Days) `,
#        q` gr05.data$Std.Dev = c(gr05.day1.sd, gr05.day2.sd, gr05.day3.sd, gr05.day4.sd, gr05.day5.sd, gr05.day6.sd, gr05.day7.sd) `, 
#        q` gr05.data$Average = c(gr05.day1.mean, gr05.day2.mean, gr05.day3.mean, gr05.day4.mean, gr05.day5.mean, gr05.day6.mean, gr05.day7.mean) `, 
#);
#
#$R->run(q` ipcc.data = data.frame(Speciation = rep("IPCC", 7), Days) `,
#        q` ipcc.data$Std.Dev = c(ipcc.day1.sd, ipcc.day2.sd, ipcc.day3.sd, ipcc.day4.sd, ipcc.day5.sd, ipcc.day6.sd, ipcc.day7.sd) `, 
#        q` ipcc.data$Average = c(ipcc.day1.mean, ipcc.day2.mean, ipcc.day3.mean, ipcc.day4.mean, ipcc.day5.mean, ipcc.day6.mean, ipcc.day7.mean) `, 
#);
#
#$R->run(q` tno.data = data.frame(Speciation = rep("TNO", 7), Days) `,
#        q` tno.data$Std.Dev = c(tno.day1.sd, tno.day2.sd, tno.day3.sd, tno.day4.sd, tno.day5.sd, tno.day6.sd, tno.day7.sd) `, 
#        q` tno.data$Average = c(tno.day1.mean, tno.day2.mean, tno.day3.mean, tno.day4.mean, tno.day5.mean, tno.day6.mean, tno.day7.mean) `, 
#);
#
#$R->run(q` uk08.data = data.frame(Speciation = rep("UK08", 7), Days) `,
#        q` uk08.data$Std.Dev = c(uk08.day1.sd, uk08.day2.sd, uk08.day3.sd, uk08.day4.sd, uk08.day5.sd, uk08.day6.sd, uk08.day7.sd) `, 
#        q` uk08.data$Average = c(uk08.day1.mean, uk08.day2.mean, uk08.day3.mean, uk08.day4.mean, uk08.day5.mean, uk08.day6.mean, uk08.day7.mean) `, 
#);
#
#$R->run(q` uk98.data = data.frame(Speciation = rep("UK98", 7), Days) `,
#        q` uk98.data$Std.Dev = c(uk98.day1.sd, uk98.day2.sd, uk98.day3.sd, uk98.day4.sd, uk98.day5.sd, uk98.day6.sd, uk98.day7.sd) `, 
#        q` uk98.data$Average = c(uk98.day1.mean, uk98.day2.mean, uk98.day3.mean, uk98.day4.mean, uk98.day5.mean, uk98.day6.mean, uk98.day7.mean) `, 
#);
#
#$R->run(q` print.data = rbind(de94.data, emep.data, gr05.data, gr95.data, ipcc.data, tno.data, uk08.data, uk98.data) `,
#        q` print.data$Std.Dev = sprintf("%.3e", print.data$Std.Dev) `,
#        q` print.data$Average = sprintf("%.3e", print.data$Average) `,
#);
#$R->run(q` write.table(print.data, file = "Solvents_Only_Ozone_time_series_daily_difference_statistics_by_speciation.txt", row.names = FALSE, quote = FALSE, sep = "\t") `);
#
#$R->run(q` plot.data = print.data %>% gather(Variable, Value, -Speciation, -Days) `,
#        q` plot.data$Value = as.numeric(as.character(plot.data$Value)) `,
#);
#my $p = $R->run(q` print(plot.data) `);
#print $p, "\n";
#$R->run(q` p = ggplot(plot.data, aes( x = Days, y = Value, colour = Speciation)) `,
#        q` p = p + geom_point() `,
#        q` p = p + facet_wrap( ~ Variable, scales = "free_y") `,
#        q` p = p + theme_hc() `,
#        q` p = p + ggtitle("Solvents Only - by Speciation") `,
#        q` p = p + theme(axis.line = element_line(colour = "black")) `,
#        q` p = p + theme(strip.text = element_text(face = "bold")) `,
#        q` p = p + theme(strip.background = element_blank()) `,
#);
#
#$R->run(q` CairoPDF(file = "Solvents_only_Ozone_statistics_plots_by_speciation.pdf", width = 11, height = 7) `,
#        q` print(p) `,
#        q` dev.off() `,
#);

$R->stop();
