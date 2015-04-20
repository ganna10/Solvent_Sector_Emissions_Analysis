library(ggplot2)
library(tidyr)
library(Cairo)
library(scales)
library(ggthemes)
library(grid)
library(dplyr)

my.colours = c("Acids" = "#cc6329", "Alcohols" = "#6c254f", "Benzene" = "#898989", "Butanes" = "#77aecc", "Chlorinated" = "#f9c500", "Esters" = "#623812", "Ethane" = "#86b650", "Ethene" = "#f36a71", "Ethers" = "#ba8b01", "Ethyne" = "#dc3522", "Formaldehyde" = "#9bb18d", "Higher alkanes" = "#0e5c28", "Ketones" = "#ef6638", "Aldehydes" = "#8ed6d2", "Other alkenes, alkynes, dienes" = "#58691b", "Other aromatics" = "#b569b3", "Others" = "#2b9eb3", "Pentanes" = "#8c1531", "Propane" = "#e7e85e", "Propene" = "#0c3f78", "Terpenes" = "#ae4901", "Toluene" = "#0352cb", "Trimethylbenzenes" = "#c9a415", "Xylenes" = "#1b695b")

mcm.data = read.csv(file = "/local/home/coates/Documents/Analysis/2014_Solvent_Sector_Emissions/Speciations/Speciations_Speciated_for_MCM.csv", header = TRUE)
mcm.data = select(mcm.data, Types, contains("to"))
mcm.data = rename(mcm.data, TNO = TNO.to.100., IPCC = IPCC.to.100., EMEP = EMEP.to.100., DE94 = DE94.to.100., GR05 = GR05.to.100., GR95 = GR95.to.100., UK98 = UK98.to.100., UK08 = UK08.to.100.)
mcm.data = gather(mcm.data, Speciation, Percentage, -Types)
mcm.data$Percentage = mcm.data$Percentage/100
mcm.data$Mechanism = rep('MCM', length(mcm.data$Percentage))

mozart.data = read.csv(file = "/local/home/coates/Documents/Analysis/2014_Solvent_Sector_Emissions/Speciations/Speciations_Speciated_for_MOZART.csv", header = TRUE)
mozart.data = select(mozart.data, Types, contains("to"))
mozart.data = rename(mozart.data, TNO = TNO.to.100., IPCC = IPCC.to.100., EMEP = EMEP.to.100., DE94 = DE94.to.100., GR05 = GR05.to.100., GR95 = GR95.to.100., UK98 = UK98.to.100., UK08 = UK08.to.100.)
mozart.data = gather(mozart.data, Speciation, Percentage, -Types)
mozart.data$Percentage = mozart.data$Percentage/100
mozart.data$Mechanism = rep('MOZART', length(mozart.data$Percentage))

radm2.data = read.csv(file = "/local/home/coates/Documents/Analysis/2014_Solvent_Sector_Emissions/Speciations/Speciations_Speciated_for_RADM2.csv", header = TRUE)
radm2.data = select(radm2.data, Types, contains("to"))
radm2.data = rename(radm2.data, TNO = TNO.to.100., IPCC = IPCC.to.100., EMEP = EMEP.to.100., DE94 = DE94.to.100., GR05 = GR05.to.100., GR95 = GR95.to.100., UK98 = UK98.to.100., UK08 = UK08.to.100.)
radm2.data = gather(radm2.data, Speciation, Percentage, -Types)
radm2.data$Percentage = radm2.data$Percentage/100
radm2.data$Mechanism = rep('RADM2', length(radm2.data$Percentage))

data = rbind(mcm.data, mozart.data, radm2.data)
data = filter(data, Types != "Total")
data
data$Types = factor(data$Types, levels = c("Acids", "Alcohols", "Aldehydes", "Benzene", "Butanes", "Chlorinated", "Esters", "Ethane", "Ethene", "Ethers", "Higher alkanes", "Ketones", "Other alkenes, alkynes, dienes", "Other aromatics", "Others", "Pentanes", "Propane", "Propene", "Terpenes", "Toluene", "Trimethylbenzenes", "Xylenes")) 

p = ggplot(data, aes(x = Mechanism, y = Percentage, fill = Types))
p = p + geom_bar(stat = "identity", position = "stack")
p = p + facet_wrap( ~ Speciation, nrow = 2)
p = p + scale_y_continuous(expand = c(0,0), labels = percent)
p = p + scale_x_discrete(expand = c(0, 0))
p = p + theme_tufte()
p = p + ggtitle("Solvent Speciations in Different Mechanisms\n")
p = p + theme(plot.title = element_text(size = 24, face = "bold", family = "NimbusSanCond"))
p = p + theme(strip.text = element_text(size = 22, face = "bold"))
p = p + theme(axis.title = element_blank())
p = p + theme(legend.title = element_blank())
p = p + theme(axis.text = element_text(family = "NimbusSanCond", face = "bold"))
p = p + theme(legend.text = element_text(family = "NimbusSanCond", size = 17))
p = p + theme(legend.key.size = unit(8, "mm"))
p = p + theme(axis.ticks.x = element_blank())
p = p + theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 0.8, vjust = 0.8))
p = p + theme(axis.text.y = element_text(size = 17))
p = p + theme(panel.margin.x = unit(7, "mm"))
p = p + scale_fill_manual(values = my.colours, limits = rev(levels(data$Types)))

Cairo(file = "Speciations_all_mechanisms.pdf",  type = "pdf", bg = "transparent", unit = "cm", width = 29, height = 20.5)
print(p)
dev.off()
