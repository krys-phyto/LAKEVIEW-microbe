
library(ggplot2)
library(tidyverse)


### Files ###
silva.tax.1 <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/LAKEVIEW-microbe/analysis/mothur/mothur.nr_v132.wang.pick.taxonomy.txt"
silva.tax.2 <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/LAKEVIEW-microbe/analysis/mothur/mothur.nr_v132.wang.pick.tax.summary.txt"

silva.tax.1 <- read.delim(silva.tax.1)
silva.tax.2 <- read.delim(silva.tax.2)

stats <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/LAKEVIEW-microbe/analysis/mothur/final.opti_mcc.groups.ave-std.summary"
stats <- read.delim(stats)

nmds.1 <- "/Users/kjkibler/Library/CloudStorage/GoogleDrive-kjkibler@wisc.edu/My Drive/LAKEVIEW-microbe/analysis/mothur/final.opti_mcc.braycurtis.0.03.lt.ave.nmds.axes"


nmds<-read.table(file=nmds.1, header=T)

nmds %>% ggplot(aes(x = axis1, y = axis2, label = label)) +
  geom_point(aes(color = location)) +
  geom_text() +
  theme_bw()

stats %>% filter(method == "ave") %>% ggplot(aes(x = group, y = invsimpson)) +
  geom_point() +
  theme_bw()

silva.tax.2_filter <- silva.tax.2 %>% filter(taxlevel == 2 |
                                               taxon == "Bacteria")
silva.tax.2_filter <- silva.tax.2_filter %>% filter(total > 1000)

silva.tax.2_filter %>% pivot_longer(c(6:13), names_to = "site", values_to = "reads") %>% 
  ggplot() +
  geom_bar(aes(x = site, y = reads, fill = taxon),
           position = "stack",
           stat = "identity") +
  facet_grid(~ taxlevel, switch = "x") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm"),
        axis.text.x = element_text(angle = 45, hjust = 1))

silva.tax.2_filter_wide <- silva.tax.2_filter %>% pivot_longer(c(6:13), names_to = "site", values_to = "reads")
