# calculate average O3 mixing ratios in each mechanism
# JC 09/11/2015

setwd("~/Documents//Analysis//2014_Solvent_Sector_Emissions/Mixing_ratios/")
so.data = read.csv(file = "solvents_only_ozone_mixing_ratios.csv")
so.data = so.data %>% gather(Speciation, O3, -Mechanism)
so.day.night.data = so.data %>% group_by(Mechanism) %>% 
  summarise(Mean.O3 = mean(O3)*1e9)

mean.NO.data = read.csv(file = "mean_NO_source_ozone_mixing_ratios.csv")
mean.NO.data = mean.NO.data %>% gather(Speciation, O3, -Mechanism)
mean.NO.day.night.data = mean.NO.data %>% group_by(Mechanism) %>% 
  summarise(Mean.O3 = mean(O3)*1e9)