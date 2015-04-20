library(dplyr)
library(tidyr)
library(ggplot2)

lm_eqn = function(m) {
  l <- list(a = format(coef(m)[1], digits = 2),
      b = format(abs(coef(m)[2]), digits = 2),
      r2 = format(summary(m)$r.squared, digits = 3));
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  as.character(as.expression(eq));                 
}

mcm.speciations = read.csv("~/Documents/Analysis/2014_Solvent_Sector_Emissions/Speciations/Speciations_Speciated_for_MCM.csv", header = TRUE)
mcm.speciations = select(mcm.speciations, Types, contains('to'))
mcm.speciations = rename(mcm.speciations, TNO = TNO.to.100., IPCC = IPCC.to.100., EMEP = EMEP.to.100., DE94 = DE94.to.100., GR05 = GR05.to.100., GR95 = GR95.to.100., UK98 = UK98.to.100., UK08 = UK08.to.100.)
mcm.speciations = filter(mcm.speciations, Types != "Total")
mcm.speciations$Mechanism = rep("MCM", length(mcm.speciations$Types))
mcm.speciations = gather(mcm.speciations, Speciation, Contribution, -Types, -Mechanism)

Ox.production = read.table("Ox_production_allocated.csv", sep = ",")
Ox.production = Ox.production[, -1]
Ox.production = rename(Ox.production, Mechanism = V2, Speciation = V3, Types = V4, Ox.Production = V5) 
mcm.Ox.production = filter(Ox.production, Mechanism == "MCM")

xylenes.data = filter(mcm.Ox.production, Types == "Xylenes")

xylenes.speciations = mcm.speciations %>% filter(Types == "Xylenes")


xylenes.data$Contribution = c(3.83, 17.9041790, 0.4004004, 0.3000000, 3.5964036, 7.9914087, 4.6810909, 4.7457288)
p = ggplot(xylenes.data, aes(x = Contribution, y = Ox.Production, colour = Speciation))
p = p + geom_point(size = 4)
p + geom_text(aes(x = 9, y = 8e8, label = lm_eqn(lm(Ox.Production ~ Contribution, xylenes.data))), colour = "black", parse = TRUE) + geom_smooth(method = "lm", se = FALSE, colour = "black")
