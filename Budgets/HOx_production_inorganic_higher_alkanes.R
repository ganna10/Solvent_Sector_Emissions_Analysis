library(tidyr)
library(dplyr)

data = read.csv(file = "HOx_allocated_production_data.csv", header = TRUE, sep = ";")
data = tbl_df(data)

inorganic = data %>% filter(Category == "Inorganic") %>% select(-Category)
inorganic = inorganic %>% spread(Mechanism, HOx.Prod) %>% mutate(MOZART.rel.diff = (MCM - MOZART)/MCM, RADM2.rel.diff = (MCM - RADM2)/MCM)

higher.alkanes = data %>% filter(Category == "Higher alkanes") %>% select(-Category)
higher.alkanes = higher.alkanes %>% spread(Mechanism, HOx.Prod) %>% mutate(MOZART.rel.diff = (MCM - MOZART)/MCM, RADM2.rel.diff = (MCM - RADM2)/MCM)
