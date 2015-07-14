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
$R->run(q` radm2 = d %>% filter(Mechanism == "RADM2") %>% select(-Mechanism) `,
        q` mozart = d %>% filter(Mechanism == "MOZART") %>% select(-Mechanism) `,
        q` mcm = d %>% filter(Mechanism == "MCM") %>% select(-Mechanism) `,
);

$R->run(q` radm2 = radm2 %>% rowwise() %>% mutate(Min = min(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) %>% mutate(Max = max(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) `,
        q` radm2 = radm2 %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `, 
        q` radm2.sd = sd(radm2$Largest.Diff) `,
        q` radm2.mean = mean(radm2$Largest.Diff) `,

        q` radm2.day1 = radm2[1:72,] `,
        q` radm2.day1.sd = sd(radm2.day1$Largest.Diff) `,
        q` radm2.day1.mean = mean(radm2.day1$Largest.Diff) `,

        q` radm2.day2 = radm2[73:144,] `,
        q` radm2.day2.sd = sd(radm2.day2$Largest.Diff) `,
        q` radm2.day2.mean = mean(radm2.day2$Largest.Diff) `,

        q` radm2.day3 = radm2[145:216,] `,
        q` radm2.day3.sd = sd(radm2.day3$Largest.Diff) `,
        q` radm2.day3.mean = mean(radm2.day3$Largest.Diff) `,

        q` radm2.day4 = radm2[217:288,] `,
        q` radm2.day4.sd = sd(radm2.day4$Largest.Diff) `,
        q` radm2.day4.mean = mean(radm2.day4$Largest.Diff) `,

        q` radm2.day5 = radm2[289:360,] `,
        q` radm2.day5.sd = sd(radm2.day5$Largest.Diff) `,
        q` radm2.day5.mean = mean(radm2.day5$Largest.Diff) `,

        q` radm2.day6 = radm2[361:432,] `,
        q` radm2.day6.sd = sd(radm2.day6$Largest.Diff) `,
        q` radm2.day6.mean = mean(radm2.day6$Largest.Diff) `,

        q` radm2.day7 = radm2[433:504,] `,
        q` radm2.day7.sd = sd(radm2.day7$Largest.Diff) `,
        q` radm2.day7.mean = mean(radm2.day7$Largest.Diff) `,
);

$R->run(q` mozart = mozart %>% rowwise() %>% mutate(Min = min(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) %>% mutate(Max = max(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) `,
        q` mozart = mozart %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `, 
        q` mozart.sd = sd(mozart$Largest.Diff) `,
        q` mozart.mean = mean(mozart$Largest.Diff) `,

        q` mozart.day1 = mozart[1:72,] `,
        q` mozart.day1.sd = sd(mozart.day1$Largest.Diff) `,
        q` mozart.day1.mean = mean(mozart.day1$Largest.Diff) `,

        q` mozart.day2 = mozart[73:144,] `,
        q` mozart.day2.sd = sd(mozart.day2$Largest.Diff) `,
        q` mozart.day2.mean = mean(mozart.day2$Largest.Diff) `,

        q` mozart.day3 = mozart[145:216,] `,
        q` mozart.day3.sd = sd(mozart.day3$Largest.Diff) `,
        q` mozart.day3.mean = mean(mozart.day3$Largest.Diff) `,

        q` mozart.day4 = mozart[217:288,] `,
        q` mozart.day4.sd = sd(mozart.day4$Largest.Diff) `,
        q` mozart.day4.mean = mean(mozart.day4$Largest.Diff) `,

        q` mozart.day5 = mozart[289:360,] `,
        q` mozart.day5.sd = sd(mozart.day5$Largest.Diff) `,
        q` mozart.day5.mean = mean(mozart.day5$Largest.Diff) `,

        q` mozart.day6 = mozart[361:432,] `,
        q` mozart.day6.sd = sd(mozart.day6$Largest.Diff) `,
        q` mozart.day6.mean = mean(mozart.day6$Largest.Diff) `,

        q` mozart.day7 = mozart[433:504,] `,
        q` mozart.day7.sd = sd(mozart.day7$Largest.Diff) `,
        q` mozart.day7.mean = mean(mozart.day7$Largest.Diff) `,
);

$R->run(q` mcm = mcm %>% rowwise() %>% mutate(Min = min(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) %>% mutate(Max = max(DE94, EMEP, GR05, GR95, IPCC, TNO, UK08, UK98)) `,
        q` mcm = mcm %>% rowwise() %>% mutate(Largest.Diff = Max - Min) `, 
        q` mcm.sd = sd(mcm$Largest.Diff) `,
        q` mcm.mean = mean(mcm$Largest.Diff) `,

        q` mcm.day1 = mcm[1:72,] `,
        q` mcm.day1.sd = sd(mcm.day1$Largest.Diff) `,
        q` mcm.day1.mean = mean(mcm.day1$Largest.Diff) `,

        q` mcm.day2 = mcm[73:144,] `,
        q` mcm.day2.sd = sd(mcm.day2$Largest.Diff) `,
        q` mcm.day2.mean = mean(mcm.day2$Largest.Diff) `,

        q` mcm.day3 = mcm[145:216,] `,
        q` mcm.day3.sd = sd(mcm.day3$Largest.Diff) `,
        q` mcm.day3.mean = mean(mcm.day3$Largest.Diff) `,

        q` mcm.day4 = mcm[217:288,] `,
        q` mcm.day4.sd = sd(mcm.day4$Largest.Diff) `,
        q` mcm.day4.mean = mean(mcm.day4$Largest.Diff) `,

        q` mcm.day5 = mcm[289:360,] `,
        q` mcm.day5.sd = sd(mcm.day5$Largest.Diff) `,
        q` mcm.day5.mean = mean(mcm.day5$Largest.Diff) `,

        q` mcm.day6 = mcm[361:432,] `,
        q` mcm.day6.sd = sd(mcm.day6$Largest.Diff) `,
        q` mcm.day6.mean = mean(mcm.day6$Largest.Diff) `,

        q` mcm.day7 = mcm[433:504,] `,
        q` mcm.day7.sd = sd(mcm.day7$Largest.Diff) `,
        q` mcm.day7.mean = mean(mcm.day7$Largest.Diff) `,
);

$R->run(q` entire.data = data.frame(Mechanism = "MCM", St.Dev = mcm.sd, Average = mcm.mean) `,
        q` entire.data = rbind(entire.data, data.frame(Mechanism = 'MOZART', St.Dev = mozart.sd, Average = mozart.mean)) `,
        q` entire.data = rbind(entire.data, data.frame(Mechanism = 'RADM2', St.Dev = radm2.sd, Average = radm2.mean)) `,
        q` entire.data$St.Dev = sprintf("%.3e", entire.data$St.Dev) `,
        q` entire.data$Average = sprintf("%.3e", entire.data$Average) `,
);
$R->run(q` write.table(entire.data, file = "Ozone_time_series_entire_difference_statistics.txt", row.names = FALSE, quote = FALSE, sep = "\t") `);

$R->set('Time', [("Day 1", "Day 2", "Day 3", "Day 4", "Day 5", "Day 6", "Day 7")]);
$R->run(q` radm2.data = data.frame(Mechanism = rep("RADM2", 7), Time) `,
        q` radm2.data$Std.Dev = c(radm2.day1.sd, radm2.day2.sd, radm2.day2.sd, radm2.day3.sd, radm2.day5.sd, radm2.day6.sd, radm2.day7.sd) `, 
        q` radm2.data$Average = c(radm2.day1.mean, radm2.day2.mean, radm2.day2.mean, radm2.day3.mean, radm2.day5.mean, radm2.day6.mean, radm2.day7.mean) `, 
);

$R->run(q` mozart.data = data.frame(Mechanism = rep("MOZART", 7), Time) `,
        q` mozart.data$Std.Dev = c(mozart.day1.sd, mozart.day2.sd, mozart.day2.sd, mozart.day3.sd, mozart.day5.sd, mozart.day6.sd, mozart.day7.sd) `, 
        q` mozart.data$Average = c(mozart.day1.mean, mozart.day2.mean, mozart.day2.mean, mozart.day3.mean, mozart.day5.mean, mozart.day6.mean, mozart.day7.mean) `, 
);

$R->run(q` mcm.data = data.frame(Mechanism = rep("MCM", 7), Time) `,
        q` mcm.data$Std.Dev = c(mcm.day1.sd, mcm.day2.sd, mcm.day2.sd, mcm.day3.sd, mcm.day5.sd, mcm.day6.sd, mcm.day7.sd) `, 
        q` mcm.data$Average = c(mcm.day1.mean, mcm.day2.mean, mcm.day2.mean, mcm.day3.mean, mcm.day5.mean, mcm.day6.mean, mcm.day7.mean) `, 
);

$R->run(q` print.data = rbind(mcm.data, mozart.data, radm2.data) `,
        q` print.data$Std.Dev = sprintf("%.3e", print.data$Std.Dev) `,
        q` print.data$Average = sprintf("%.3e", print.data$Average) `,
);
#$R->run(q` write.table(print.data, file = "Ozone_time_series_daily_difference_statistics.txt", row.names = FALSE, quote = FALSE, sep = "\t") `);

$R->run(q` plot.data = print.data %>% gather(Variable, Value, -Mechanism, -Time) `);
#my $p = $R->run(q` print(plot.data) `);
#print $p, "\n";
$R->run(q` p = ggplot(plot.data, aes( x = Time, y = Value, colour = Mechanism)) `,
        q` p = p + geom_point() `,
        q` p = p + facet_wrap( ~ Variable, scales = "free_y") `,
        q` p = p + theme_hc() `,
        q` p = p + theme(axis.line = element_line(colour = "black")) `,
        q` p = p + theme(strip.text = element_text(face = "bold")) `,
        q` p = p + theme(strip.background = element_blank()) `,
);

$R->run(q` CairoPDF(file = "Ozone_statistics_plots.pdf", width = 11, height = 7) `,
        q` print(p) `,
        q` dev.off() `,
);

$R->stop();
