
library(SIBER)

set.seed(1)

# load in the dataset
data1<- read.csv("data/processed/SIBER/SIBER_data.csv")
#
# create the siber object
siber.RLB <- createSiberObject(data1)

# change ellipse color
palette(c("royalblue2","orange"))


# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.68, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")


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

legend("topleft", legend = c("Otospermophilus Pre-LGM", "Sylvilagus Pre-LGM"),
       col = palette(c("royalblue2","orange")), pch = 1, 
       bty = "n", cex = 0.8)

legend("topright", legend = c("Otospermophilus Post-LGM", "Sylvilagus Post-LGM"),
       col = palette(c("royalblue2", "orange")), pch = 2,
       bty = "n", cex = 0.8)






# Summary stats
group.ML <- groupMetricsML(siber.RLB)
print(group.ML)
#1.1 = Old audubonii, #1.2 = Old beecheyi
#2.1 = young audubonii, #2.2 = young beecheyi

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


