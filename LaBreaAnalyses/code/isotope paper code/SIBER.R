library(tidyverse)
library(gridExtra)
library(SIBER)
# Create final dataset ----
# originally input through SIBER, hence the file names. No longer using SIBER

# load in the isotope dataset
data_raw<- read.csv("data/processed/SIBER/SIBER_raw_final.csv", strip.white=T)

# create sample names for matching isotope file to ages file
samples <- paste("UCIAMS", data_raw$UCIAMS_Number) # this is the final set of samples with isotope data

# match with the sample ages
all_calibrated_ages <- read.csv('output/OxCal/final oxcal models/AllAges_forinput.csv', header=T) # this file stores the calibrated age statistics for each age
all_calibrated_ages$trimmedName <- unlist(lapply(strsplit(all_calibrated_ages$Name, " R_"), '[[', 1))
sample_median_ages <- as.data.frame(cbind(samples, all_calibrated_ages[match(samples, all_calibrated_ages$trimmedName), 'Unmodelled._BP_median']))
colnames(sample_median_ages)[2] <- 'median_age'  
sample_median_ages$median_age <- as.numeric(sample_median_ages$median_age)

# add ages to the file
data_raw <- cbind(data_raw, median_age=sample_median_ages$median_age)

# Create age groups
time_group <- vector(length=nrow(data_raw))
time_group[which(data_raw$median_age > 11500)] <- "Pleistocene"
time_group[which(data_raw$median_age < 11500)] <- "Holocene"

# add age groups to raw data and create primary data file (data1)
data1 <- as.data.frame(cbind(data_raw, time_group))

#write.csv(data1, file="data/processed/SIBER/SIBER_data.csv", row.names=F)


# plot the isotope data - using ggplot2 ####
rlbPalette <- palette(c("royalblue2","darkorange"))
rlb_data <- data1 %>% mutate(Taxon = factor(Taxon), 
                             time_group = factor(time_group),
                             d13C = del13C_permil, 
                             d15N = del15N_permil,
                             .keep = "unused") 

first.plot <- ggplot(data = rlb_data, 
                     aes(x = d13C, 
                         y = d15N)) + 
  geom_point(aes(colour = Taxon, shape = time_group), size = 3) +
  scale_colour_manual(labels = c("Otospermophilus", "Sylvilagus"), 
                      values=rlbPalette) +
  scale_shape_manual(labels = c("Holocene", "Pleistocene"), 
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
  group_by(Taxon, time_group) %>% 
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
                             fill = Taxon), 
             color = "black", shape = 22, size = 5,
             alpha = 0.7, show.legend = FALSE) +
  scale_fill_manual(values=rlbPalette)

print(second.plot)


p.ell <- 0.68
ellipse.plot <- first.plot + 
  stat_ellipse(aes(Taxon = interaction(Taxon, time_group), 
                   fill = Taxon, 
                   color = Taxon), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  scale_fill_manual(values=rlbPalette) 

print(ellipse.plot)

# test out new plots ----


# add probability ellipses for Pleistocene and Holocene data 
# (taxa differentiated)

p.ell <- 0.68
alldata.taxa.ellipse.plot <- first.plot + 
  stat_ellipse(data = rlb_data, 
               aes(x = d13C, 
                   y = d15N, color=Taxon), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  stat_ellipse(data = rlb_data[-c(1:5),], 
               aes(x = d13C, 
                   y = d15N, color=Taxon), 
               linetype = 2,
               alpha = 0, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  scale_fill_manual(values=rlbPalette) 

grDevices::cairo_pdf("output/isotope paper final/Figure2_SIBERplots_Dec2021_JLB_bytaxa_time.pdf", width=7, height=4)
print(alldata.taxa.ellipse.plot)
dev.off()

# add probability ellipses for Pleistocene and Holocene data 
# (combining taxa)

p.ell <- 0.68
alldata.time.ellipse.plot <- first.plot + 
  stat_ellipse(data = rlb_data[-c(1:5)], 
               aes(x = d13C, 
                   y = d15N), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  stat_ellipse(data = rlb_data[c(1:5),], 
               aes(x = d13C, 
                   y = d15N), 
               linetype = 2,
               alpha = 0.5, 
               level = p.ell,
               type = "norm",
               geom = "polygon")  

grDevices::cairo_pdf("output/isotope paper final/Figure2_SIBERplots_Dec2021_JLB_bytime.pdf", width=7, height=4)
print(alldata.time.ellipse.plot)
dev.off()

# merged - Final Figure 2 plot ----
# will alter in Illustrator afterwards

new.first.plot <- ggplot(data = rlb_data, 
                     aes(x = d13C, 
                         y = d15N)) + 
  lims(x = c(-23, -17), y= c(0,11)) +
  geom_point(aes(colour = Taxon, shape = time_group), size = 3) +
  scale_colour_manual(labels = c("Otospermophilus", "Sylvilagus"), 
                      values=rlbPalette) +
  scale_shape_manual(labels = c("Holocene", "Pleistocene"), 
                     values=c(16,17)) +
  ylab(expression({delta}^15*N~'value ('~'\u2030'~', AIR)')) +
  xlab(expression({delta}^13*C~'value ('~'\u2030'~', VPDB)')) + 
  theme_classic() +
  theme(text = element_text(size=14),
        axis.ticks.length = unit(0.15, "cm")) + 
  labs(colour = "Taxon", shape="Time Period") 

new.alldata.taxa.ellipse.plot <- new.first.plot + 
  stat_ellipse(data = rlb_data, 
               aes(x = d13C, 
                   y = d15N, color=Taxon), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  stat_ellipse(data = rlb_data[-c(1:5),], 
               aes(x = d13C, 
                   y = d15N, color=Taxon), 
               linetype = 2,
               alpha = 0, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  scale_fill_manual(values=rlbPalette) 

new.alldata.time.ellipse.plot <- new.first.plot + 
  stat_ellipse(data = rlb_data[-c(1:5)], 
               aes(x = d13C, 
                   y = d15N), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  stat_ellipse(data = rlb_data[c(1:5),], 
               aes(x = d13C, 
                   y = d15N), 
               linetype = 2,
               alpha = 0.5, 
               level = p.ell,
               type = "norm",
               geom = "polygon") 

figure <- ggarrange(new.alldata.time.ellipse.plot, 
                    new.alldata.taxa.ellipse.plot, 
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)

grDevices::cairo_pdf("output/isotope paper final/Figure2_SIBERplots_Jan2022_JLB_both.pdf", width=12, height=4)
print(figure)
dev.off()



# Summary stats ----

cn.manova <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon, data1)
summary(cn.manova)

# factor two of the grouping variables the run t-tests
data2 <- data1
data2$Taxon <- as.factor(data2$Taxon)
data2$time_group <- as.factor(data2$time_group)

# t.test
c.t <- t.test(del13C_permil~Taxon, data=data2)
c.t
n.t <- t.test(del15N_permil~Taxon, data=data2)
n.t

# Just pleistocene data
data3 <- data2[which(data2$time_group=="Pleistocene"),]
cn.manova.p <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon, data3)
summary(cn.manova.p)

# t.test
c.t.p <- t.test(del13C_permil~Taxon, data=data3)
c.t.p
n.t.p <- t.test(del15N_permil~Taxon, data=data3)
n.t.p

# try adding sample age to the manova
cnt.manova <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon + median_age, data1)
summary(cnt.manova)

cnt.manova.p <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon + median_age, data3)
summary(cnt.manova.p)

# try adding sample age to the manova
cnt.manova.interact <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon * median_age, data1)
summary(cnt.manova.interact)

cnt.manova.p.interact <- manova(cbind(del13C_permil, del15N_permil) ~ Taxon * median_age, data3)
summary(cnt.manova.p.interact)

# try adding a color gradient
# Note: this doesn't add valuable info and is  confusing
# plot the isotope data - using ggplot2 ####
rlbPalette <- palette(c("royalblue2","darkorange"))
rlb_data <- data1 %>% mutate(Taxon = factor(Taxon), 
                             time_group = factor(time_group),
                             d13C = del13C_permil, 
                             d15N = del15N_permil,
                             .keep = "unused") 

first.plot <- ggplot(data = rlb_data, 
                     aes(x = d13C, 
                         y = d15N)) + 
  geom_point(aes(colour = Taxon, shape = time_group, alpha=median_age), size = 3) +
  scale_colour_manual(labels = c("Otospermophilus", "Sylvilagus"), 
                      values=rlbPalette) +
  scale_shape_manual(labels = c("Holocene", "Pleistocene"), 
                     values=c(16,17)) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme_classic() +
  theme(text = element_text(size=14),
        axis.ticks.length = unit(0.15, "cm")) + 
  labs(colour = "Taxon", shape="Time Period", alpha="Median Age") 

print(first.plot) 

p.ell <- 0.68
ellipse.plot <- first.plot + 
  stat_ellipse(aes(Taxon = interaction(Taxon, time_group), 
                   fill = Taxon, 
                   color = Taxon), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  scale_fill_manual(values=rlbPalette) 

print(ellipse.plot)

# print Figure 2 plot ----
# will alter in Illustrator afterwards
grDevices::cairo_pdf("output/isotope paper final/Figure2_SIBERplots_Nov2021_JB_agegradient.pdf", width=8, height=6)
ellipse.plot
dev.off()



###
## Everything following is old, unused code ----
###

# old anovas
cn.aov <- aov(group~iso1+iso2, data=data1)
summary(cn.aov)

# cn.aov <- aov(iso1 ~ group + iso2, data = data2)
# summary(cn.aov)
# TukeyHSD(cn.aov, which = 'group')

#SIBER Stats ----

## Create data frame for SIBER analysis
# First create a df with both groups (taxa) and communities (Holo vs Pleisto). The sample sizes are too low for Holo
siber.data <- data1
siber.data <- siber.data %>% 
  rename(iso1 = del13C_permil,
        iso2 = del15N_permil,
        community = time_group,
        group = Taxon) %>%
    select(iso1, iso2, group, community)

siber.obj <- createSiberObject(siber.data)

# create new object that groups Pleisto and Holo into one community
siber.data.1comm <- siber.data
siber.data.1comm$community <- "Quaternary"
siber.1comm.obj <- createSiberObject(siber.data.1comm)

# Then create separate Pleisto dataframe
siber.data.pleisto <- data1 %>% 
  rename(iso1 = del13C_permil,
         iso2 = del15N_permil,
         community = time_group,
         group = Taxon) %>%
  filter(median_age > 11500) %>%
  select(iso1, iso2, group, community)
siber.pleisto.obj <- createSiberObject(siber.data.pleisto)

community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")

# plot data with all four separate groups and communities
par(mfrow=c(1,1))
plotSiberObject(siber.obj,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

group.ML.siber.obj <- groupMetricsML(siber.obj)
print(group.ML.siber.obj)

# plot data with two groups and 1 Quaternary community
par(mfrow=c(1,1))
plotSiberObject(siber.1comm.obj,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

group.ML.siber.1comm.obj <- groupMetricsML(siber.1comm.obj)
print(group.ML.siber.1comm.obj)


# plot data with two groups and 1 Pleisto community
par(mfrow=c(1,1))
plotSiberObject(siber.pleisto.obj,
                ax.pad = 2, 
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = T, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030')
)

group.ML.siber.pleisto.obj <- groupMetricsML(siber.pleisto.obj)
print(group.ML.siber.pleisto.obj)

# fit Bayesian multivariate normal distributions to each group in the dataset. do this for all four separately ----
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
ellipses.posterior <- siberMVN(siber.obj, parms, priors)


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

