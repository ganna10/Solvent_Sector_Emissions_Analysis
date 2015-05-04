library(tidyr)
library(dplyr)

data = read.csv(file = "Ox_production_allocated.csv", header = FALSE)
data = tbl_df(data)
colnames(data) = c("Mechanism", "Speciation", "Type", "Ox")

ch4 = data %>% filter(Type == "Methane") %>% select(-Type)
ch4 = ch4 %>% spread(Mechanism, Ox)
ch4 = ch4 %>% mutate(MOZART.rel.diff = (MCM - MOZART)/MCM, RADM2.rel.diff = (MCM - RADM2)/MCM)

co = data %>% filter(Type == "CO") %>% select(-Type)
co = co %>% spread(Mechanism, Ox)
write.table(co, file = "Background_CO.csv", sep = ",", row.name = FALSE, quote = FALSE)
co = co %>% mutate(MOZART.rel.diff = (MCM - MOZART)/MCM, RADM2.rel.diff = (MCM - RADM2)/MCM)
