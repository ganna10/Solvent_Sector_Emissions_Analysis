library(tidyr)
library(dplyr)

data = read.csv(file = "HO2x_cumulative_production.csv", sep = ";", header = TRUE)
data = tbl_df(data)
#data = data %>% group_by(Speciation) %>% spread(Mechanism, HO2x.Production)

inorganic = data %>% filter(Type == "Inorganic") %>% spread(Mechanism, HO2x.Production) %>% mutate(MOZART.rel.diff = (MCM - MOZART)/MCM, RADM2.rel.diff = (MCM - RADM2)/MCM)
print.data.frame(inorganic)

higher.alkanes = data %>% filter(Type == "Higher alkanes") 
higher.alkanes = higher.alkanes %>% spread(Mechanism, HO2x.Production) %>% mutate(MOZART.rel.diff = (MCM - MOZART)/MCM, RADM2.rel.diff = (MCM - RADM2)/MCM)
#print.data.frame(higher.alkanes)

methane = data %>% filter(Type == "Methane") 
methane = methane %>% spread(Mechanism, HO2x.Production) %>% mutate(MOZART.rel.diff = (MCM - MOZART)/MCM, RADM2.rel.diff = (MCM - RADM2)/MCM)
print.data.frame(methane)

ketones = data %>% filter(Type == "Ketones") 
ketones = ketones %>% spread(Mechanism, HO2x.Production) %>% mutate(MOZART.rel.diff = (MCM - MOZART)/MCM, RADM2.rel.diff = (MCM - RADM2)/MCM)
#print.data.frame(ketones)
