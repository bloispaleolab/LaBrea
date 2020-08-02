library(tidyverse)
library(scales)

# read in master file that stores names/assemblages for trait analyses
dat <- read.delim("data/processed/Master_trait_taxon_names.txt", sep="\t", stringsAsFactors = F)
dat <- dat[,-which(colnames(dat)=='include')]

# redo the "include" to match new "include" file
include <- read.csv(file="data/processed/includes.file.csv", header=T)

dat <- cbind(dat, include[match(dat$RevisedName, include$RevisedName),'Include'])
colnames(dat)[ncol(dat)] <- "include"

# read in master mammal_nisp and process to account for Sylvilagus ----

mammals_nisp <- read.delim("data/processed/mammal_nisp.txt", sep="\t", stringsAsFactors = F)
mammals_nisp[is.na(mammals_nisp)] <- 0 # convert the NAs to 0
colsForBoxes <- grep("Box", colnames(mammals_nisp)) 

# Estimate Sylvilagus spp contributions and add onto NISP for the named species ----

# calculate Sylvilagus ratios and determin addition factors
temp_bach <- mammals_nisp[which(mammals_nisp$RevisedName == "Sylvilagus bachmani"),]
temp_audu <- mammals_nisp[which(mammals_nisp$RevisedName == "Sylvilagus audubonii"),]
temp_df <- colSums(rbind(temp_bach, temp_audu)[,colsForBoxes], na.rm=T)
prop_bach <- temp_bach[,colsForBoxes]/temp_df
prop_audu <- temp_audu[,colsForBoxes]/temp_df

temp_sylv_sp <- mammals_nisp[which(mammals_nisp$RevisedName == "Sylvilagus sp"),4:7]
total_new_bach <- round(prop_bach*temp_sylv_sp,0)
total_new_audu <- round(prop_audu*temp_sylv_sp,0) 

# Add on these genus Sylvilagus proportions to individual species then delete the Sylvilagus sp. row
mammals_nisp[which(mammals_nisp$RevisedName == "Sylvilagus bachmani"),colsForBoxes] <- mammals_nisp[which(mammals_nisp$RevisedName == "Sylvilagus bachmani"),colsForBoxes] + total_new_bach
mammals_nisp[which(mammals_nisp$RevisedName == "Sylvilagus audubonii"),colsForBoxes] <- mammals_nisp[which(mammals_nisp$RevisedName == "Sylvilagus audubonii"),colsForBoxes] + total_new_audu
mammals_nisp <- mammals_nisp[-which(mammals_nisp$RevisedName == "Sylvilagus sp"),]


# set the scales for the sensitivity analysis
tax_scales <- c("known", "cf", "estimated")
spat_scales <- c("local", "regional", "continental")

# create the final trait files for each assemblage
for (i in 1:length(tax_scales)){
  for (j in 1:length(spat_scales)){
    
    # determine the final taxon list
    if (i == 1){tax_rows <- which(dat$BestTaxData == tax_scales[i])}
    if (i == 2){tax_rows <- c(which(dat$BestTaxData == tax_scales[1]), which(dat$BestTaxData == tax_scales[i]))}
    if (i == 3){tax_rows <- seq(1:nrow(dat))}
   
    if (j == 1){spat_rows <- which(dat$BestSpatData == spat_scales[j])}
    if (j == 2){spat_rows <- c(which(dat$BestSpatData == spat_scales[1]), which(dat$BestSpatData == spat_scales[j]))}
    if (j == 3){spat_rows <- seq(1:nrow(dat))}
    
    taxon_list <- unique(dat$Species[intersect(tax_rows, spat_rows)])
    taxon_list <- taxon_list[which(!is.na(taxon_list))]
    taxon_list <- taxon_list[which(dat$include[match(taxon_list, dat$Species)]=="y")]

    # now that we have a taxon list, create a blank dataframe to store the trait data
    trait_dat <- as.data.frame(matrix(nrow=length(taxon_list), ncol=55))
    
    #we want to associate the bioclim data for each species
    for (s in 1:length(taxon_list)){
    
    # read in bioclim and ndvis summary file, arrange, and store in trait dat
      bioclim <- read.csv(paste0("data/processed/range_wide_bioclim_stats/range_wide_bioclim_stats_", taxon_list[s], ".csv"), header=T)
      colnames(bioclim)[1] <- "statistic"
      bioclim_long <- reshape(bioclim, idvar = "statistic", varying=colnames(bioclim)[2:ncol(bioclim)], direction = "long")
      colnames(bioclim_long)[2:3] <- c("bioclim_variable", "value")    
      bioclim_long$bioclim_variable <- gsub("1_5m_", "", bioclim_long$bioclim_variable)  
      bioclim_colnames <- gsub(".1_5m", "", rownames(bioclim_long))
      
      trait_dat[s, 1:length(bioclim_colnames)]<- bioclim_long$value
      colnames(trait_dat)[1:length(bioclim_colnames)] <- bioclim_colnames
      
      ndvi <- read.csv(paste0("data/processed/range_wide_ndvi_stats/range_wide_ndvi_stats_", taxon_list[s], ".csv"), header=T)
      ndvi_transposed <- t(ndvi[,2])
      colnames(ndvi_transposed) <- paste0(ndvi[,1], "_ndvi")
      
      trait_dat[s, 49:54]<- ndvi_transposed
      colnames(trait_dat)[49:54] <- colnames(ndvi_transposed)
      
      }

    # read in Pantheria trait data and store
    pan <- read.delim("data/raw/PanTHERIA_1-0_WR05_Aug2008.txt", sep="\t", header=T)
    
    # pull out PantheriaID for each taxon in taxon list
    pantheria_list <- dat$PantheriaID[match(taxon_list, dat$Species)]
    trait_dat[,ncol(trait_dat)] <- pan[match(pantheria_list, pan$MSW05_Binomial),'X5.1_AdultBodyMass_g']
    colnames(trait_dat)[ncol(trait_dat)] <- 'AdultBodyMass_g'
    trait_dat <- cbind(taxon_list, trait_dat)
    
    write.csv(trait_dat, file=paste0("data/processed/community_trait_files/traits_", tax_scales[i], "_", spat_scales[j], ".csv"), row.names=F)

    # pull out the nisp data for all the taxa in the focal assemblage ----
    
    # match back the Species name with RevisedName
    RevisedName_list <- dat$RevisedName[match(taxon_list, dat$Species)]
    matched_nisp <- mammals_nisp[match(RevisedName_list, mammals_nisp$RevisedName),]
    total_nisp <- colSums(matched_nisp[,colsForBoxes])
    
    # calculate weighted trait averages ----
    
    # for each trait, first remove the taxa with no trait data
    traits <- colnames(trait_dat)[-1]
    
    trait_weighted_mean <- as.data.frame(matrix(NA, nrow=length(traits), ncol=length(colsForBoxes)))
    colnames(trait_weighted_mean) <- colnames(matched_nisp)[colsForBoxes]
    rownames(trait_weighted_mean) <- traits
    
    for (t in 1:length(traits)){
      trait_weighted_mean[t,] <- colSums(matched_nisp[,colsForBoxes]*trait_dat[,t+1])/total_nisp
    }
    
    rownames(trait_weighted_mean) <- traits
    
    write.csv(trait_weighted_mean, file=paste0("data/processed/assemblage_trait_mean_files/trait_weighted_means_", tax_scales[i], "_", spat_scales[j], ".csv"), row.names=T)
    
  }
}


