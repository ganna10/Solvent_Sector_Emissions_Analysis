library(tidyr)
library(dplyr)
library(ggplot2)
data = data.frame(Mechanism = c("MCM v3.2", "MOZART-4", "RADM2"), CO = c(1.26e9, 1.34e9, 1.47e9), stdev = c(6.17e7, 1.43e7, 1.35e7))
data = data %>% mutate(ymax = CO + stdev, ymin = CO - stdev)

p = ggplot(data, aes(x = Mechanism, y = CO)) + geom_point()
p + geom_errorbar(aes(ymin = ymin, ymax = ymax))
p = p + geom_errorbar(aes(ymin = ymin, ymax = ymax))
