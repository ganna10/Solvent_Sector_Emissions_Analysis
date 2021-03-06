setwd("~/Documents//Analysis//2014_Solvent_Sector_Emissions//Speciations")

# my.colours = c("Acids" = "#0352cb", "Alcohols" = "#011e4b", "Benzene" = "#0e5c28", "Butanes" = "#6e254f", "Chlorinated" = "#ef6638", "CO" = "#898989", "Esters" = "#0357d8", "Ethane" = "#a15daf", "Ethene" = "#dfb100", "Ethers" = "#2b9eb3", "Higher alkanes" = "#38103d", "Ketones" = "#02388b", "Aldehydes" = "#0499fc", "Other alkenes, alkynes, dienes" = "#f9c500", "Other aromatics" = "#216105", "Others" = "#898989", "Pentanes" = "#8f2ac9", "Propane" = "#5c1b54", "Propene" = "#796000", "Terpenes" = "#b99300", "Toluene" = "#0b6956", "Trimethylbenzenes" = "#109c00", "Xylenes" = "#0c734b", "CO" = "#6d6537", "Methane" = "#77aecc", "Inorganic" = "#000000")

data = read.csv(file = "Mapping_initial_species_to_types_all_mechanisms_plot_speciations.csv", header = TRUE)
data = filter(data, Type != "Total")
data$Type = factor(data$Type, levels = c("Ethane", "Propane", "Butanes", "Pentanes", "Higher alkanes", "Ethene", "Propene", "Terpenes", "Other alkenes, alkynes, dienes", "Benzene", "Toluene", "Trimethylbenzenes", "Xylenes", "Other aromatics", "Acids", "Alcohols", "Aldehydes", "Esters", "Ethers", "Ketones", "Chlorinated"))
data$Speciation = factor(data$Speciation, levels = c("TNO", "IPCC", "EMEP", "DE94", "GR95", "GR05", "UK98", "UK08"))
data= data %>% gather(Mechanism, Contribution, -Type, -Speciation)
data = mutate(data, mechanism = factor(Mechanism, labels = c("MCM v3.2", "MOZART-4", "RADM2")))

my.colours = c("Acids" = "grey57", "Alcohols" = "grey89", "Benzene" = "grey71", "Butanes" = "grey37", "Chlorinated" = "grey94", "CO" = "#898989", "Esters" = "grey55", "Ethane" = "grey3", "Ethene" = "grey87", "Ethers" = "grey69", "Higher alkanes" = "grey15", "Ketones" = "grey5", "Aldehydes" = "grey4", "Other alkenes, alkynes, dienes" = "grey28", "Other aromatics" = "grey83", "Others" = "#000000", "Pentanes" = "grey88", "Propane" = "grey72", "Propene" = "grey19", "Terpenes" = "grey96", "Toluene" = "grey40", "Trimethylbenzenes" = "grey79", "Xylenes" = "grey14", "CO" = "#6d6537", "Methane" = "#77aecc", "Inorganic" = "#000000")

p = ggplot(data %>% arrange(Type), aes(x = mechanism, y = Contribution, fill = Type))
p = p + geom_bar(stat = "identity", position = "stack")
p = p + facet_wrap( ~ Speciation, nrow = 2)
p = p + scale_y_continuous(expand = c(0,0), labels = percent)
p = p + scale_x_discrete(expand = c(0, 0))
p = p + theme_tufte()
p = p + theme(strip.text = element_text(face = "bold"))
p = p + theme(axis.title = element_blank())
p = p + theme(legend.title = element_blank())
# p = p + theme(axis.text = element_text(family = "NimbusSanCond", face = "bold"))
# p = p + theme(legend.text = element_text(family = "NimbusSanCond"))
p = p + theme(legend.key.size = unit(8, "mm"))
p = p + theme(axis.ticks.x = element_blank())
p = p + theme(axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 0.8))
p = p + theme(panel.margin.x = unit(7, "mm"))
p = p + scale_fill_manual(values = my.colours, limits = rev(levels(data$Type)), guide = guide_legend(ncol = 1))

CairoPDF(file = "Speciations_all_mechanisms_b_w.pdf", width = 10, height = 7.5)
print(p)
dev.off()
