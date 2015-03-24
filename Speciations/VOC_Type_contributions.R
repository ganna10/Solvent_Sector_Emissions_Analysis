library(ggplot2)
library(tidyr)
library(Cairo)
library(scales)
library(ggthemes)
library(grid)
library(dplyr)

my.colours = c("Acids" = "#cc6329", "Alcohols" = "#6c254f", "Benzene" = "#898989", "Butanes" = "#77aecc", "Chlorinated" = "#f9c500", "Esters" = "#623812", "Ethane" = "#86b650", "Ethene" = "#f36a71", "Ethers" = "#ba8b01", "Ethyne" = "#dc3522", "Formaldehyde" = "#9bb18d", "Higher alkanes" = "#0e5c28", "Ketones" = "#ef6638", "Aldehydes" = "#8ed6d2", "Other alkenes, alkynes, dienes" = "#58691b", "Other aromatics" = "#b569b3", "Others" = "#2b9eb3", "Pentanes" = "#8c1531", "Propane" = "#e7e85e", "Propene" = "#0c3f78", "Terpenes" = "#ae4901", "Toluene" = "#0352cb", "Trimethylbenzenes" = "#c9a415", "Xylenes" = "#1b695b")

mcm.data = read.csv(file = "/local/home/coates/Documents/Analysis/2014_Solvent_Sector_Emissions/Speciated_Mechanisms_after_Passant-MCM.csv", header = TRUE)
mcm.data = select(mcm.data, Types, contains("Percent"))
mcm.data = rename(mcm.data, TNO = Percent.TNO.Emissions, IPCC = Percent.IPCC.Emissions, EMEP = Percent.EMEP.Emissions, DE94 = Percent.DE94.Emissions, GR05 = Percent.GR05.Emissions, GR95 = Percent.GR95.Emissions, UK98 = Percent.UK98.Emissions, UK08 = Percent.UK08.Emissions)
mcm.data = mcm.data[-c(22,23),]
mcm.data = gather(mcm.data, Speciation, Percentage, -Types)
mcm.data$Mechanism = rep('MCM', length(mcm.data$Percentage))

data = mcm.data
data$Types = factor(data$Types, levels = c("Acids", "Alcohols", "Aldehydes", "Benzene", "Butanes", "Chlorinated", "Esters", "Ethane", "Ethene", "Ethers", "Higher alkanes", "Ketones", "Other alkenes, alkynes, dienes", "Other aromatics", "Pentanes", "Propane", "Propene", "Terpenes", "Toluene", "Trimethylbenzenes", "Xylenes")) 

p = ggplot(data, aes(x = Speciation, y = Percentage, fill = Types))
p = p + geom_bar(stat = "identity")
p = p + scale_y_continuous(expand = c(0,0), labels = percent)
p = p + scale_x_discrete(expand = c(0, 0))
p = p + theme_tufte()
p = p + theme(axis.title = element_blank())
p = p + theme(legend.title = element_blank())
p = p + theme(axis.text = element_text(family = "NimbusSanCond", face = "bold"))
p = p + theme(legend.text = element_text(family = "NimbusSanCond"))
p = p + theme(axis.ticks.x = element_blank())
p = p + theme(axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 0.8))
p = p + scale_fill_manual(values = my.colours, limits = rev(levels(data$Types)))

CairoPDF(file = "MCM_percent_all_VOC.pdf")
print(p)
dev.off()
