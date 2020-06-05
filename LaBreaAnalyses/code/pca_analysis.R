
# Perform PCA
tax_scales <- c("known", "cf", "estimated")
spat_scales <- c("local", "regional", "continental")

# adjust these values based on what scale you want to look at
i <- 1
j <- 1  

# read in bioclim and other trait data
trait_means <- read.csv(file=paste0("data/processed/assemblage_trait_mean_files/trait_weighted_means_", tax_scales[i], "_", spat_scales[j], ".csv"))

# pull out just the median climate trait values
medianT <- trait_means[grep("median", trait_means$X),]
rownames(medianT) <- medianT$X
medianT <- medianT[,-1] # remove rownames
medianT <- medianT[, -which(colnames(medianT) == "Box_14")] # remove Box 14
medianT <- t(medianT)

# remove seasonality (bio_4, bio_15) and just look at temp/precip variables
colsToRemove <- match(c('median_bio_4', 'median_bio_15'), colnames(medianT))
medianT_trimmed <- medianT[,-colsToRemove]

# perform a PCA on the median trait values
pca_median <- prcomp(medianT_trimmed, center=T, scale=T)
biplot(pca_median)
pc_scores <- pca_median$x
pca_median$rotation # shows the loadings

# read in the d180 box-level climate data
d180 <- read.csv("data/processed/Interpolation/Box_summary.csv")
d180 <- d180[,-1] # NATE DELETE THIS LINE OF CODE ONCE YOUVE FIXED 

# discard Box 14 data
d180 <- d180[-which(d180$Deposit==14),]

# get pc_scores in same order as d180
pc_scores <- pc_scores[c(1,3,2),]

merged <- cbind(d180, pc_scores)
merged$Deposit <- as.factor(merged$Deposit)

# plot pca results
par(mfrow=c(4,2), mar=c(4,4,1,1)+0.1)
plot(P23_d18O_med~Deposit, data=merged)
plot(-PC1~Deposit, data=merged)
plot(P23_d18O_mean~Deposit, data=merged)
plot(-PC1~Deposit, data=merged)
plot(P23_d18O_med~Deposit, data=merged)
plot(PC2~Deposit, data=merged)
plot(P23_d18O_mean~Deposit, data=merged)
plot(PC2~Deposit, data=merged)

# plot original trait mean results
trait_means_arranged <- trait_means
# remove box 14
trait_means_arranged <- trait_means_arranged[,-grep("14", colnames(trait_means_arranged))]

rownames(trait_means_arranged) <- trait_means$X
trait_means_arranged <- trait_means_arranged[,-1]
trait_means_arranged <- t(trait_means_arranged)

# reorder to match deposits
trait_means_arranged <- trait_means_arranged[c(1,3,2),]
trait_means_arranged <- as.data.frame(trait_means_arranged)

# plot just temp (bio 1)
par(mfrow=c(2,2), mar=c(4,4,1,1)+0.1)
plot(P23_d18O_med~Deposit, data=merged)
plot(-trait_means_arranged$median_bio_1~merged$Deposit)

plot(P23_d18O_med~Deposit, data=merged)
plot(-PC1~Deposit, data=merged)
