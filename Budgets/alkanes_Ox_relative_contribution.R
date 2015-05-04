library(tidyr)
library(dplyr)

data = read.csv(file = "Ox_production_allocated.csv", header = FALSE)
data = tbl_df(data)
colnames(data) = c("Mechanism", "Speciation", "Type", "Ox")
Total.data = data %>% group_by(Mechanism, Speciation) %>% mutate(Total.Ox = sum(Ox)) %>% select(-Type, -Ox)
Total.data = unique(Total.data)

alkanes = data %>% filter(Type == "Ethane" | Type == "Propane" | Type == "Butanes" | Type == "Pentanes" | Type == "Higher alkanes") %>% select(-Type) 
alkanes = alkanes %>% group_by(Mechanism, Speciation) %>% mutate(Ox = sum(Ox))
alkanes = unique(alkanes)

alkenes = data %>% filter(Type == "Ethene" | Type == "Propene" | Type == "Terpenes" | Type == "Other alkanes, alkynes, dienes") %>% select(-Type) 
alkenes = alkenes %>% group_by(Mechanism, Speciation) %>% mutate(Ox = sum(Ox))
alkenes = unique(alkenes)

aromatics = data %>% filter(Type == "Benzene" | Type == "Toluene" | Type == "Xylenes" | Type == "Trimethylbenzenes" | Type == "Other aromatics") %>% select(-Type) 
aromatics = aromatics %>% group_by(Mechanism, Speciation) %>% mutate(Ox = sum(Ox))
aromatics = unique(aromatics)

oxygenated = data %>% filter(Type == "Acids" | Type == "Alcohols" | Type == "Aldehydes" | Type == "Esters" | Type == "Ethers" | Type == "Ketones") %>% select(-Type) 
oxygenated = oxygenated %>% group_by(Mechanism, Speciation) %>% mutate(Ox = sum(Ox))
oxygenated = unique(oxygenated)

#chlorinated = data %>% filter(Type == "Chlorinated") %>% select(-Type) 
#chlorinated = chlorinated %>% group_by(Mechanism, Speciation) %>% mutate(Ox = sum(Ox))
#chlorinated = unique(chlorinated)

methane = data %>% filter(Type == "Methane") %>% select(-Type) 
methane = methane %>% group_by(Mechanism, Speciation) %>% mutate(Ox = sum(Ox))
methane = unique(methane)

co = data %>% filter(Type == "CO") %>% select(-Type) 
co = co %>% group_by(Mechanism, Speciation) %>% mutate(Ox = sum(Ox))
co = unique(co)

inorganic = data %>% filter(Type == "Inorganic") %>% select(-Type) 
inorganic = inorganic %>% group_by(Mechanism, Speciation) %>% mutate(Ox = sum(Ox))
inorganic = unique(inorganic)

alka.cont = select(alkanes, Mechanism, Speciation) 
alka.cont$Alk.Ox = alkanes$Ox
alka.cont$total.Ox = Total.data$Total.Ox
alka.cont = alka.cont %>% mutate(alka.percent = Alk.Ox / total.Ox * 100) %>% select(-Alk.Ox, -total.Ox)
print.data.frame(alka.cont)

#alke.cont = select(alkenes, Mechanism, Speciation) 
#alke.cont$Alk.Ox = alkenes$Ox
#alke.cont$total.Ox = Total.data$Total.Ox
#alke.cont = alke.cont %>% mutate(alke.percent = Alk.Ox / total.Ox * 100) %>% select(-Alk.Ox, -total.Ox)
#print.data.frame(alke.cont)

arom.cont = select(aromatics, Mechanism, Speciation) 
arom.cont$Alk.Ox = aromatics$Ox
arom.cont$total.Ox = Total.data$Total.Ox
arom.cont = arom.cont %>% mutate(arom.percent = Alk.Ox / total.Ox * 100) %>% select(-Alk.Ox, -total.Ox)
print.data.frame(arom.cont)

oxy.cont = select(oxygenated, Mechanism, Speciation) 
oxy.cont$Alk.Ox = oxygenated$Ox
oxy.cont$total.Ox = Total.data$Total.Ox
oxy.cont = oxy.cont %>% mutate(oxy.percent = Alk.Ox / total.Ox * 100) %>% select(-Alk.Ox, -total.Ox)
print.data.frame(oxy.cont)

#chl.cont = select(chlorinated, Mechanism, Speciation) 
#chl.cont$Alk.Ox = chlorinated$Ox
#chl.cont$total.Ox = Total.data$Total.Ox
#chl.cont = chl.cont %>% mutate(chl.percent = Alk.Ox / total.Ox * 100) %>% select(-Alk.Ox, -total.Ox)
#print.data.frame(chl.cont)

ch4.cont = select(methane, Mechanism, Speciation) 
ch4.cont$Alk.Ox = methane$Ox
ch4.cont$total.Ox = Total.data$Total.Ox
ch4.cont = ch4.cont %>% mutate(ch4.percent = Alk.Ox / total.Ox * 100) %>% select(-Alk.Ox, -total.Ox)
print.data.frame(ch4.cont)

co.cont = select(co, Mechanism, Speciation) 
co.cont$Alk.Ox = co$Ox
co.cont$total.Ox = Total.data$Total.Ox
co.cont = co.cont %>% mutate(co.percent = Alk.Ox / total.Ox * 100) %>% select(-Alk.Ox, -total.Ox)
print.data.frame(co.cont)

inorganic.cont = select(inorganic, Mechanism, Speciation) 
inorganic.cont$Alk.Ox = inorganic$Ox
inorganic.cont$total.Ox = Total.data$Total.Ox
inorganic.cont = inorganic.cont %>% mutate(inorganic.percent = Alk.Ox / total.Ox * 100) %>% select(-Alk.Ox, -total.Ox)
print.data.frame(inorganic.cont)
