library(scales)

tax_scales <- c("known", "cf", "estimated")
spat_scales <- c("local", "regional", "continental")

# read in the d180 box-level climate data
d180 <- read.csv("data/processed/Interpolation/Box_summary.csv")

# discard Box 14 data
d180 <- d180[-which(d180$Deposit==14),]
# reorder by time
d180$Deposit <- as.factor(d180$Deposit)
d180 <- d180[c(1,3,2),]

# adjust these values based on what scale you want to look at
i <- 1
j <- 1  

xvars <- c("P23_d18O_mean", "median_ndvi")
yvars <- c("median_bio_1", "median_bio_12", "AdultBodyMass_g")

for (k in 1:length(xvars)){
  for (m in 1:length(yvars)){
    
  pdf(file=paste0("output/new sensitivity/sensitivy_results_", yvars[m], "-", xvars[k], ".pdf"))
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
        cols <- c(match(c(xvars[k], yvars[m]), colnames(merged)), which(colnames(merged)=="Deposit"))
        merged_trim <- merged[,cols]
        
        # plot just temp (bio 1)
        plot(merged_trim[,2] ~ merged_trim[,1], pch=16, 
             xlim=scales::expand_range(range(merged_trim[,1]), mul = 0.1), 
             ylim=scales::expand_range(range(merged_trim[,2]), mul = 0.1),
             xlab="", ylab="")
        with(merged_trim, text(merged_trim[,2] ~ merged_trim[,1], labels=Deposit, pos=4))
        mtext(paste0("taxonomic scale = ", tax_scales[i], "; spatial scale = ", spat_scales[j]), side=3, line=0, cex=0.5)
        mtext(colnames(merged_trim)[1], side=1, line=2.5, cex=0.9)
        mtext(colnames(merged_trim)[2], side=2, line=2.5, cex=0.9)
      }
      }
  dev.off()
  }
}