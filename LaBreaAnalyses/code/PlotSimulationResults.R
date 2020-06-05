tax_scales <- c("known", "cf", "estimated")
spat_scales <- c("local", "regional", "continental")

# read in the d180 box-level climate data
d180 <- read.csv("data/processed/Interpolation/Box_summary.csv")
d180 <- d180[,-1] # NATE DELETE THIS LINE OF CODE ONCE YOUVE FIXED 

# discard Box 14 data
d180 <- d180[-which(d180$Deposit==14),]
# reorder by time
d180$Deposit <- as.factor(d180$Deposit)
d180 <- d180[c(1,3,2),]

# adjust these values based on what scale you want to look at
i <- 1
j <- 1  

pdf(file="output/sensitivy_results_bioclim4.pdf")

par(mfrow=c(3,3), mar=c(4,4,1,1)+0.1)

for (i in 1:length(tax_scales)){
  for (j in 1:length(spat_scales)){

    # read in bioclim and other trait data
    trait_means <- read.csv(file=paste0("data/processed/assemblage_trait_mean_files/trait_weighted_means_", tax_scales[i], "_", spat_scales[j], ".csv"))
    
    # remove box 14
    trait_means_arranged <- trait_means[,-grep("14", colnames(trait_means))]
    
    rownames(trait_means_arranged) <- trait_means$X
    trait_means_arranged <- trait_means_arranged[,-1]
    trait_means_arranged <- t(trait_means_arranged)
    
    # merge with climate data
    merged <- cbind(d180, trait_means_arranged)
    
    
    # plot just temp (bio 1)
    plot(median_bio_1 ~ P23_d18O_med, data=merged, pch=16)
    with(merged, text(median_bio_1 ~ P23_d18O_med, labels=Deposit, pos=c(1, 4, 2)))
    mtext(paste0("taxonomic scale = ", tax_scales[i], "; spatial scale = ", spat_scales[j]), side=3, line=0, cex=0.5)
  }
}

dev.off()