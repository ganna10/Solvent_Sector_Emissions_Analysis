library(ggplot2)
library(tidyr)
library(Cairo)
library(scales)
library(ggthemes)
library(grid)
library(dplyr)

my.colours = c("Acids" = "#0352cb", "Alcohols" = "#011e4b", "Benzene" = "#0e5c28", "Butanes" = "#6e254f", "Chlorinated" = "#ef6638", "CO" = "#898989", "Esters" = "#0357d8", "Ethane" = "#a15daf", "Ethene" = "#dfb100", "Ethers" = "#2b9eb3", "Higher alkanes" = "#38103d", "Ketones" = "#02388b", "Aldehydes" = "#0499fc", "Other alkenes, alkynes, dienes" = "#f9c500", "Other aromatics" = "#216105", "Others" = "#898989", "Pentanes" = "#8f2ac9", "Propane" = "#5c1b54", "Propene" = "#796000", "Terpenes" = "#b99300", "Toluene" = "#0b6956", "Trimethylbenzenes" = "#109c00", "Xylenes" = "#0c734b", "CO" = "#6d6537", "Methane" = "#77aecc", "Inorganic" = "#000000")

mcm.data = read.csv(file = "/local/home/coates/Documents/Analysis/2014_Solvent_Sector_Emissions/Speciations/Speciations_Speciated_for_MCM.csv", header = TRUE)
mcm.data = select(mcm.data, Types, contains("to"))
mcm.data = rename(mcm.data, TNO = TNO.to.100., IPCC = IPCC.to.100., EMEP = EMEP.to.100., DE94 = DE94.to.100., GR05 = GR05.to.100., GR95 = GR95.to.100., UK98 = UK98.to.100., UK08 = UK08.to.100.)
mcm.data = gather(mcm.data, Speciation, Percentage, -Types)
mcm.data$Percentage = mcm.data$Percentage/100
mcm.data$Mechanism = rep('MCM', length(mcm.data$Percentage))
#mcm.data = read.csv(file = "/local/home/coates/Documents/Analysis/2014_Solvent_Sector_Emissions/Speciations/Speciated_Mechanisms_after_Passant-MCM.csv", header = TRUE)cm.data$Mechanism = rep('MCM', length(mcm.data$Percentage))
#mcm.data = select(mcm.data, Types, contains("Percent"))
#mcm.data = rename(mcm.data, TNO = Percent.TNO.Emissions, IPCC = Percent.IPCC.Emissions, EMEP = Percent.EMEP.Emissions, DE94 = Percent.DE94.Emissions, GR05 = Percent.GR05.Emissions, GR95 = Percent.GR95.Emissions, UK98 = Percent.UK98.Emissions, UK08 = Percent.UK08.Emissions)
#mcm.data = mcm.data[-c(22,23),]
#mcm.data = gather(mcm.data, Speciation, Percentage, -Types)

data = mcm.data
data = filter(data, Types != "Total")
data$Types = factor(data$Types, levels = c("Ethane", "Propane", "Butanes", "Pentanes", "Higher alkanes", "Ethene", "Propene", "Terpenes", "Other alkenes, alkynes, dienes", "Benzene", "Toluene", "Trimethylbenzenes", "Xylenes", "Other aromatics", "Acids", "Alcohols", "Aldehydes", "Esters", "Ethers", "Ketones", "Chlorinated", "Others"))

p = ggplot(data, aes(x = Speciation, y = Percentage, fill = Types, order = Types))
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
