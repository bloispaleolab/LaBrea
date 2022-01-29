library(tidyverse)
library(readr)

# plot the isotope data - using ggplot2 ####
rlb_data<- read_csv("data/processed/SIBER/rlb_data2_jlb.csv", trim_ws=TRUE)
rlb_data$Taxon <- as.factor(rlb_data$Taxon)
rlb_data$Time <- as.factor(rlb_data$Time)
rlb_data$Group <- as.factor(rlb_data$Group)

rlb_data<- rlb_data %>%
  arrange(Group) %>%
  mutate(Taxon = fct_relevel(Taxon, 
                             "Bison", "Camelops", "Equus", "Mammut", "Paramylodon", "Microtus", "Neotoma", "Otospermophilus", "Sylvilagus", "Thomomys"), Group = fct_relevel(Group, "Megaherbivore", "Herbivore"), Time = fct_relevel(Time, "Pleistocene", "Holocene"))

labelColors <- c("Bison" = "gray20",
                 "Camelops" = "gray30",
                 "Equus" = "gray40", 
                 "Mammut" = "gray50",
                 "Paramylodon" = "gray60",
                 "Microtus" = "red3",
                 "Neotoma" = "forestgreen",
                 "Otospermophilus" = "royalblue2",
                 "Sylvilagus" = "darkorange",
                 "Thomomys" = "mediumorchid4")

first.plot <- ggplot(data = rlb_data, 
                     aes(x = d13C, y = d15N)) +
                       geom_point(aes(colour = Taxon, 
                                      shape = Time, 
                                      size = Group)) +
                       scale_size_manual(values = c(1, 2)) +
                       scale_color_manual(values=labelColors) +
  ylab(expression(paste(delta^{15}, "N (\u2030) values"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030) values"))) + 
  theme_classic() +
  theme(text = element_text(size=12),
        axis.ticks.length = unit(0.15, "cm"),
        legend.title = element_text(size=12), 
        legend.key.size = unit(0.4, 'cm')) + 
  labs(colour = "Taxon") 

grDevices::cairo_pdf("output/isotope paper final/Figure5_v1_noerrorbars.pdf", width=6, height=4)
  print(first.plot) 
dev.off()

# error bars
sbg <- rlb_data %>% 
  group_by(Taxon, Time) %>% 
  summarise(count = n(),
            mC = mean(d13C), 
            sdC = sd(d13C), 
            mN = mean(d15N), 
            sdN = sd(d15N))

second.plot <- first.plot+
  geom_errorbar(data = sbg, 
                mapping = aes(x = mC, y = mN,
                              ymin = mN - 1.96*sdN, 
                              ymax = mN + 1.96*sdN), 
                width = 0, color="darkgray", alpha=0.5) +
  geom_errorbarh(data = sbg, 
                 mapping = aes(x = mC, y = mN,
                               xmin = mC - 1.96*sdC,
                               xmax = mC + 1.96*sdC),
                 height = 0, color="darkgray", alpha=0.5) + 
  geom_point(data = sbg, aes(x = mC, 
                             y = mN,
                             fill = Taxon), 
             color = "darkgray", shape = 22, size = 4,
             alpha = 0.5, show.legend = FALSE) +
  scale_fill_manual(values=labelColors)

grDevices::cairo_pdf("output/isotope paper final/Figure5_v2_witherrorbars.pdf", width=6, height=4)
  print(second.plot) 
dev.off()
