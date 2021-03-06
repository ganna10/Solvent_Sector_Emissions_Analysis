library(ggplot2)
library(scales)
library(dplyr)
library(ggthemes)
library(Cairo)

Label = c("Public Power", "Residential Combustion", "Industrial Combustion", "Industrial Processes", "Fossil Fuel Processing", "Solvents", "Road Transport", "Non-road Transport", "Waste", "Agriculture")
Ratio = c(.022, .167, .081, .049, .049, .431, .146, .02, .013, .022)
data = data.frame(Label, Ratio)
data = transform(data, Label = reorder(Label, Ratio))

plot = ggplot(data, aes(x = Label, y = Ratio))
plot = plot + geom_bar(stat = "identity")
plot = plot + coord_flip()
plot = plot + theme_tufte()
plot = plot + theme(axis.title = element_blank())
plot = plot + theme(axis.ticks.y = element_blank())
plot = plot + scale_y_continuous(labels = percent, expand = c(0,0))
plot = plot + scale_x_discrete(expand = c(0,0))
plot = plot + theme(axis.text.x = element_text(family = "NimbusSanCond", size = 16.5, face = "bold"))
plot = plot + theme(axis.text.y = element_text(family = "NimbusSanCond", size = 18, face = "bold"))

Cairo(file = "Sector_conributions.pdf", type = "pdf", bg = "transparent", unit = "cm", width = 16.5, height = 20.3)
print(plot)
dev.off() 
