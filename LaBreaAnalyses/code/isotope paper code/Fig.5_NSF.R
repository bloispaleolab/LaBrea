library(tidyverse)

# plot the isotope data - using ggplot2 ####
rlb_data<- read.csv("data/processed/SIBER/rlb_data2.csv")

first.plot <- ggplot(data = rlb_data, 
                     aes(x = d13C, 
                         y = d15N)) + 
  geom_point(aes(colour = Taxon), size = 3, alpha = 1) +
  scale_colour_manual(labels = c("Bison (n = 20)", "Camelops (n = 13)", "Equus (n = 23)", "Mammut (n = 7)", "Microtus - Ple (n = 1)", "Neotoma - Ple (n = 1)",
  "Otospermophilus - Hol (n = 3)", "Otospermophilus - Ple (n = 24)", "Paramylodon (n = 14)", "Sylvilagus - Hol (n = 2)", "Sylvilagus - Ple (n = 40)", "Thomomys - Ple (n = 1)"), 
  values= c("Bison (n = 20)" = "tan", "Camelops (n = 13)" = "grey", "Equus (n = 23)" = "brown", 
 "Mammut (n = 7)" = "maroon", "Microtus - Ple (n = 1)" = "black", "Neotoma - Ple (n = 1)" = "green",
 "Otospermophilus - Hol (n = 3)" = "darkblue", "Otospermophilus - Ple (n = 24)" = "royalblue2", 
 "Paramylodon (n = 14)" = "yellow", "Sylvilagus - Hol (n = 2)" = "darkorange3", 
  "Sylvilagus - Ple (n = 40)" = "orange", "Thomomys - Ple (n = 1)" = "purple")) +

  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme_classic() +
  theme(text = element_text(size=14),
        axis.ticks.length = unit(0.15, "cm")) + 
  labs(colour = "Taxon") 
 
print(first.plot) 


# error bars
sbg <- rlb_data %>% 
  group_by(Taxon, time_group) %>% 
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
                width = 0) +
  geom_errorbarh(data = sbg, 
                 mapping = aes(x = mC, y = mN,
                               xmin = mC - 1.96*sdC,
                               xmax = mC + 1.96*sdC),
                 height = 0) + 
  geom_point(data = sbg, aes(x = mC, 
                             y = mN,
                             fill = Taxon), 
             color = "black", shape = 22, size = 5,
             alpha = 0.7, show.legend = FALSE) +
  scale_fill_manual(values=c("Bison (n = 20)" = "tan", "Camelops (n = 13)" = "grey", "Equus (n = 23)" = "brown", 
  "Mammut (n = 7)" = "maroon", "Microtus - Ple (n = 1)" = "black", "Neotoma - Ple (n = 1)" = "green",
  "Otospermophilus - Hol (n = 3)" = "darkblue", "Otospermophilus - Ple (n = 24)" = "royalblue2", 
  "Paramylodon (n = 14)" = "yellow", "Sylvilagus - Hol (n = 2)" = "darkorange3", 
  "Sylvilagus - Ple (n = 40)" = "orange", "Thomomys - Ple (n = 1)" = "purple"))

print(second.plot)


