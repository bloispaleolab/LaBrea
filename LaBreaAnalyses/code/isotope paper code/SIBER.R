
library(SIBER)
library(tidyverse)

set.seed(1)

# Create SIBER dataset ----
# load in the isotope dataset
data_raw<- read.csv("data/processed/SIBER/SIBER_raw_final.csv", strip.white=T)
# create sample names for matching isotope file to ages file
samples <- paste("UCIAMS", data_raw$UCIAMS_Number) # this is the final set of samples with isotope data

all_calibrated_ages <- read.csv('output/OxCal/final oxcal models/AllAges_forinput.csv', header=T) # this file stores the calibrated age statistics for each age
all_calibrated_ages$trimmedName <- unlist(lapply(strsplit(all_calibrated_ages$Name, " R_"), '[[', 1))
sample_median_ages <- as.data.frame(cbind(samples, all_calibrated_ages[match(samples, all_calibrated_ages$trimmedName), 'Unmodelled._BP_median']))
colnames(sample_median_ages)[2] <- 'median_age'  
sample_median_ages$median_age <- as.numeric(sample_median_ages$median_age)

data_raw <- cbind(data_raw, median_age=sample_median_ages$median_age)

# turn this into a SIBER data file:
iso1 <- data_raw$del13C_permil
iso2 <- data_raw$del15N_permil
group <- vector(length=nrow(data_raw))
group[which(data_raw$Taxon=="Otospermophilus")] <- 1
group[which(data_raw$Taxon=="Sylvilagus")] <- 2

community <- vector(length=nrow(data_raw))
community[which(data_raw$median_age > 11500)] <- 1
community[which(data_raw$median_age < 11500)] <- 2

data1 <- as.data.frame(cbind(iso1, iso2, group, community))
data2 <- data1
data2$group <- as.factor(data2$group)
data2$community <- as.factor(data2$community)

write.csv(data1, file="data/processed/SIBER/SIBER_data.csv", row.names=F)
rm(list = c('iso1','iso2', 'group', 'community'))

# Read in SIBER dataset  and create SIBER object ----

#data1 <- read.csv(file="data/processed/SIBER/SIBER_data.csv", header=T)
siber.RLB <- createSiberObject(data1) # create the siber object


# default SIBER plots - using SIBER plotting ####

# change ellipse color
palette(c("royalblue2","darkorange"))

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.68, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20") # change color here to be organe vs blue? 

#pdf("output/isotope paper final/Figure2_SIBER.pdf", width=5, height=4)
par(mfrow=c(1,1), mar=c(5,5,4,1)+0.01)
plotSiberObject(siber.RLB,
                ax.pad = 2, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = T, group.hull.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
               
)

legend("bottomleft", legend = c("Otospermophilus Pre-LGM", "Sylvilagus Pre-LGM"),
       col = palette(c("royalblue2","darkorange")), pch = 1, 
       bty = "n", cex = 0.8)

legend("topright", legend = c("Otospermophilus Post-LGM", "Sylvilagus Post-LGM"),
       col = palette(c("royalblue2", "darkorange")), pch = 2,
       bty = "n", cex = 0.8)
dev.off()


# custom SIBER plots - using ggplot2 ####
rlbPalette <- palette(c("royalblue2","darkorange"))
rlb_data <- data1 %>% mutate(group = factor(group), 
                             community = factor(community),
                             d13C = iso1, 
                             d15N = iso2,
                             .keep = "unused") 

first.plot <- ggplot(data = rlb_data, 
                     aes(x = d13C, 
                         y = d15N)) + 
  geom_point(aes(colour = group, shape = community), size = 3) +
  scale_colour_manual(labels = c("Otospermophilus", "Sylvilagus"), 
                      values=rlbPalette) +
  scale_shape_manual(labels = c("Pleistocene", "Holocene"), 
                     values=c(16,17)) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme_classic() +
  theme(text = element_text(size=14),
        axis.ticks.length = unit(0.15, "cm")) + 
  labs(colour = "Taxon", shape="Time Period") 
 
print(first.plot) 


# error bars
sbg <- rlb_data %>% 
  group_by(group, community) %>% 
  summarise(count = n(),
            mC = mean(d13C), 
            sdC = sd(d13C), 
            mN = mean(d15N), 
            sdN = sd(d15N))

second.plot <- first.plot +
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
                             fill = group), 
             color = "black", shape = 22, size = 5,
             alpha = 0.7, show.legend = FALSE) +
  scale_fill_manual(values=rlbPalette)

print(second.plot)


p.ell <- 0.68
ellipse.plot <- first.plot + 
  stat_ellipse(aes(group = interaction(group, community), 
                   fill = group, 
                   color = group), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  scale_fill_manual(values=rlbPalette) 


print(ellipse.plot)

# print Figure 2 plot ----
# will alter in Illustrator afterwards
grDevices::cairo_pdf("output/isotope paper final/Figure2_SIBERplots_Nov2021_JB.pdf", width=8, height=6)
ellipse.plot
dev.off()


# Summary stats ----
# iso1=del13C_permil
# iso2=del15N_permil

# t.test
c.t <- t.test(iso1~group, data=data2)
n.t <- t.test(iso2~group, data=data2)

# cn.aov <- aov(iso1 ~ group + iso2, data = data2)
# summary(cn.aov)
# TukeyHSD(cn.aov, which = 'group')

#SIBER Stats -- Not used, doesn't work except the first few lines !!----

group.ML <- groupMetricsML(siber.RLB)
print(group.ML)
#1.1 = Pleistocene squirrels, #1.2 = Pleistocene rabbits
#2.1 = Holocene squirrels, #2.2 = Holocene rabbits

# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-2

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.RLB, parms, priors)


# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)


# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

# extract the posterior means
mu.post <- extractPosteriorMeans(siber.RLB, ellipses.posterior)

# calculate the corresponding distribution of layman metrics
layman.B <- bayesianLayman(mu.post)

# --------------------------------------
# Visualise the first community
# --------------------------------------
siberDensityPlot(layman.B[[1]], xticklabels = colnames(layman.B[[1]]), 
                 bty="L", ylim = c(0,20))

# add the ML estimates (if you want). Extract the correct means 
# from the appropriate array held within the overall array of means.
comm1.layman.ml <- laymanMetrics(siber.RLB$ML.mu[[1]][1,1,],
                                 siber.RLB$ML.mu[[1]][1,2,]
)
points(1:6, comm1.layman.ml$metrics, col = "red", pch = "x", lwd = 2)



# extract the posterior means
mu.post <- extractPosteriorMeans(siber.example, ellipses.posterior)

# calculate the corresponding distribution of layman metrics
layman.B <- bayesianLayman(mu.post)


# --------------------------------------
# Visualise the second community
# --------------------------------------
siberDensityPlot(layman.B[[2]], xticklabels = colnames(layman.B[[2]]), 
                 bty="L", ylim = c(0,20))

# add the ML estimates. (if you want) Extract the correct means 
# from the appropriate array held within the overall array of means.
comm2.layman.ml <- laymanMetrics(siber.RLB$ML.mu[[2]][1,1,],
                                 siber.RLB$ML.mu[[2]][1,2,]
)
points(1:6, comm2.layman.ml$metrics, col = "red", pch = "x", lwd = 2)



# --------------------------------------
# Alternatively, pull out TA from both and aggregate them into a 
# single matrix using cbind() and plot them together on one graph.
# --------------------------------------

# go back to a 1x1 panel plot
par(mfrow=c(1,1))

siberDensityPlot(cbind(layman.B[[1]][,"TA"], layman.B[[2]][,"TA"]),
                 xticklabels = c("Community 1", "Community 2"), 
                 bty="L", ylim = c(0,20),
                 las = 1,
                 ylab = "TA - Convex Hull Area",
                 xlab = "")


