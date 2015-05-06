library(tidyr)
library(dplyr)

data = read.csv(file = "HOx_allocated_production_data.csv", header = TRUE)
data = tbl_df(data)

data %>% group_by(Speciation) %>% spread(Mechanism, HOx.Prod)

inorganic = data %>% filter(Category == "Inorganic") %>% select(-Category)
inorganic = inorganic %>% spread(Mechanism, HOx.Prod) %>% mutate(MOZART.rel.diff = (MCM - MOZART)/MCM, RADM2.rel.diff = (MCM - RADM2)/MCM)
#inorganic

higher.alkanes = data %>% filter(Category == "Higher alkanes") %>% select(-Category)
higher.alkanes = higher.alkanes %>% spread(Mechanism, HOx.Prod) %>% mutate(MOZART.rel.diff = (MCM - MOZART)/MCM, RADM2.rel.diff = (MCM - RADM2)/MCM)
#higher.alkanes

methane = data %>% filter(Category == "Methane") %>% select(-Category)
methane = methane %>% spread(Mechanism, HOx.Prod) %>% mutate(MOZART.rel.diff = (MCM - MOZART)/MCM, RADM2.rel.diff = (MCM - RADM2)/MCM)
#methane

ketones = data %>% filter(Category == "Ketones") %>% select(-Category)
ketones = ketones %>% spread(Mechanism, HOx.Prod) %>% mutate(MOZART.rel.diff = (MCM - MOZART)/MCM, RADM2.rel.diff = (MCM - RADM2)/MCM)
#ketones
