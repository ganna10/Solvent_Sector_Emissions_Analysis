library(tidyr)
library(dplyr)

d = read.csv("Solvents_Only_mixing_ratios.csv", sep = ";", header = TRUE)
d = tbl_df(d)
d = d %>% select(-mechanism)

mcm = d %>% filter(Mechanism == "MCM") %>% select(-Mechanism)
mcm.sd = mcm %>% group_by(Time) %>% summarise(St.Dev = sd(Mixing.Ratio))
mean(mcm.sd$St.Dev)

mozart = d %>% filter(Mechanism == "MOZART") %>% select(-Mechanism)
mozart.sd = mozart %>% group_by(Time) %>% summarise(St.Dev = sd(Mixing.Ratio))
mean(mozart.sd$St.Dev)

radm2 = d %>% filter(Mechanism == "RADM2") %>% select(-Mechanism)
radm2.sd = radm2 %>% group_by(Time) %>% summarise(St.Dev = sd(Mixing.Ratio))
mean(radm2.sd$St.Dev)
