# compare ozone mixing ratios with different urban area sizes
# Version 0: Jane Coates 19/2/2016

setwd("~/Documents//Analysis//2014_Solvent_Sector_Emissions//Area_Sensitivity")
data <- read.csv("O3_comparison_data.csv")
data <- tbl_df(data)
data

d <- data %>%
  filter(Run == "Area_Sensitivity")

my.colours = c("TNO" = "#000000", "IPCC" = "#2b9eb3", "EMEP" = "#ef6638", "DE94" = "#f9c500", "GR95" = "#6c254f", "GR05" = "#0352cb", "UK98" = "#0e5c28", "UK08" = "#b569b3")

plot <- ggplot(d, aes(x = Time, y = Mixing.Ratio*1e9, colour = Speciation))
plot <- plot + geom_line(size = 1)
plot <- plot + facet_wrap(Run ~ Mechanism)
plot <- plot + plot_theme()
plot <- plot + theme(legend.title = element_blank())
plot <- plot + theme(legend.position = "top")
plot <- plot + ylab("O3 Mixing Ratio (ppbv)")
plot <- plot + scale_colour_manual(values = my.colours)
plot

CairoPDF(file = "Area_Sensitivity_O3_sensitivity.pdf", width = 7, height = 10)
print(plot)
dev.off()

# difference statistics
entire.diffs <- d %>%
  select(-Run) %>%
  mutate(Mixing.Ratio = Mixing.Ratio * 1e9) %>%
  spread(Speciation, Mixing.Ratio, drop = FALSE) %>%
  rowwise() %>%
  group_by(Time, Mechanism) %>%
  mutate(Max.Diff = max(DE94, EMEP, GR05, GR95, IPCC, UK98, UK08) - min(DE94, EMEP, GR05, GR95, IPCC, UK98, UK08)) %>%
  select(Time, Mechanism, Max.Diff) %>%
  rowwise() %>%
  group_by(Mechanism) %>%
  mutate(Std.Dev = sd(Max.Diff), Mean = mean(Max.Diff)) %>%
  select(-Time, -Max.Diff)
alldiffs <- unique.data.frame(entire.diffs)

timed <- d %>%
  select(-Run) %>%
  mutate(Mixing.Ratio = Mixing.Ratio * 1e9) %>%
  spread(Speciation, Mixing.Ratio, drop = FALSE) %>%
  rowwise() %>%
  group_by(Time, Mechanism) %>%
  mutate(Max.Diff = max(DE94, EMEP, GR05, GR95, IPCC, TNO, UK98, UK08) - min(DE94, EMEP, GR05, GR95, IPCC, TNO, UK98, UK08)) %>%
  select(Time, Mechanism, Max.Diff) %>%
  spread(Mechanism, Max.Diff, drop = FALSE)

day1 <- timed[1:72,]
day1$Day <- rep("Day 1", length(day1$Time))
day2 <- timed[73:144,]
day2$Day <- rep("Day 2", length(day2$Time))
day3 <- timed[145:216,]
day3$Day <- rep("Day 3", length(day3$Time))
day4 <- timed[217:288,]
day4$Day <- rep("Day 4", length(day4$Time))
day5 <- timed[289:360,]
day5$Day <- rep("Day 5", length(day5$Time))
day6 <- timed[361:432,]
day6$Day <- rep("Day 6", length(day6$Time))
day7 <- timed[433:504,]
day7$Day <- rep("Day 7", length(day7$Time))

get_mean_sd <- function (df) {
  data <- df %>%
    gather(mechanism, max.diff, -Day, -Time) %>% 
    group_by(mechanism) %>%
    mutate(Std.Dev = sd(max.diff), Mean = mean(max.diff)) %>%
    select(-Time, -max.diff)
  return(unique.data.frame(data))
}
alldays <- rbind(get_mean_sd(day1), get_mean_sd(day2), get_mean_sd(day3), get_mean_sd(day4), get_mean_sd(day5), get_mean_sd(day6), get_mean_sd(day7))
alldays

write.table(alldays, file = "alldays.csv", sep = ",", row.names = FALSE, quote = FALSE)
write.table(alldiffs, file = "alldiffs.csv", sep = ",", row.names = FALSE, quote = FALSE)
