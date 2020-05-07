library(tidyverse)
library(scales)

# read in master mammal_nisp and trait data ----
mammals_nisp <- read.delim("data/processed/mammal_nisp.txt", sep="\t", stringsAsFactors = F)
mammals_nisp[is.na(mammals_nisp)] <- 0 # convert the NAs to 0

taxonomy_traits <- read.delim("data/processed/TaxonomyTraits_Reg_est_Bioclim.txt", sep="\t", stringsAsFactors = F)

# merge trait data with nisp ----
dat <- cbind(mammals_nisp, taxonomy_traits[match(mammals_nisp$RevisedName, taxonomy_traits$RevisedName),-c(1:5)])
colsForBoxes <- grep("Box", colnames(dat)) 

# Estimate Sylvilagus spp contributions ----
# calculate Sylvilagus ratios and determin addition factors
temp_bach <- dat[which(dat$RevisedName == "Sylvilagus bachmani"),]
temp_audu <- dat[which(dat$RevisedName == "Sylvilagus audubonii"),]
temp_df <- colSums(rbind(temp_bach, temp_audu)[,colsForBoxes], na.rm=T)
prop_bach <- temp_bach[,colsForBoxes]/temp_df
prop_audu <- temp_audu[,colsForBoxes]/temp_df

temp_sylv_sp <- dat[which(dat$RevisedName == "Sylvilagus sp"),]
total_new_bach <- round(prop_bach*temp_sylv_sp[,colsForBoxes],0)
total_new_audu <- round(prop_audu*temp_sylv_sp[,colsForBoxes],0) 

# Add on these genus Sylvilagus proportions to individual species then delete the Sylvilagus sp. row
dat[which(dat$RevisedName == "Sylvilagus bachmani"),colsForBoxes] <- dat[which(dat$RevisedName == "Sylvilagus bachmani"),colsForBoxes] + total_new_bach
dat[which(dat$RevisedName == "Sylvilagus audubonii"),colsForBoxes] <- dat[which(dat$RevisedName == "Sylvilagus audubonii"),colsForBoxes] + total_new_audu
dat <- dat[-which(dat$RevisedName == "Sylvilagus sp"),]

# calculate weighted trait averages ----

# for each trait, first remove the taxa with no trait data
traits <- c("BodyMass_Mean", "Precip_Mean", "Temp_Mean", "NDVI_Mean", "Temp_SD",	"Temp_max",	"Temp_min"
)
trait_weighted_mean <- as.data.frame(matrix(NA, nrow=length(traits), ncol=length(colsForBoxes)))
colnames(trait_weighted_mean) <- colnames(dat)[colsForBoxes]
rownames(trait_weighted_mean) <- traits

for (i in 1:length(traits)){
  dat_temp <- dat[-which(is.na(dat[,traits[i]])),]
  trait_nisp <- colSums(dat_temp[,colsForBoxes])
  trait_weighted_mean[i,] <- colSums(dat_temp[,traits[i]]*dat_temp[,colsForBoxes])/trait_nisp
}


