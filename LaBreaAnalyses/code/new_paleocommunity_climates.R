## Infer paleo-community climates ----
library(tidyverse)
library(readr)

## Read in data and assign names ----

# read in La Brea mammal data
nisp <- read_delim(file="data/processed/mammal_nisp.txt", delim="\t")

# read in contemporary climate ranges
load(file="data/gbif_occurrences/final/species_clim_stats.RData") #species_clim_stats
load(file="data/gbif_occurrences/final/species_occs_final.RData") #species_occs_final
load(file="data/gbif_occurrences/final/species_clim_final.RData") #species_clim_final

NatureServeNameList <- c('Mustela frenata',
                 'Mephitis mephitis',
                 'Spilogale gracilis',
                 'Microtus californicus',
                 'Peromyscus californicus',
                 'Peromyscus maniculatus',
                 'Neotoma macrotis',
                 'Reithrodontomys megalotis',
                 'Onychomys torridus',
                 'Sylvilagus audubonii',
                 'Sylvilagus bachmani',
                 'Otospermophilus beecheyi',
                 'Neotamias merriami',
                 'Dipodomys agilis',
                 'Chaetodipus californicus',
                 'Thomomys bottae',
                 'Scapanus latimanus',
                 'Sorex ornatus')

# load species name matches
matches <- read_delim(file="data/processed/Master_trait_taxon_names.txt", delim="\t")

# Arrange NISP data ----
# allocate genus nisp based on ratios
nisp_revised <- nisp[-grep("Peromyscus sp", nisp$RevisedName),]
nisp_revised <- nisp_revised[-grep("Sylvilagus sp", nisp_revised$RevisedName),]

# Peromyscus
# calculate ratio, californicus/maniculatus:
nisp_Pc <- nisp[grep("Peromyscus californicus", nisp$RevisedName),]
nisp_Pm <- nisp[grep("Peromyscus maniculatus", nisp$RevisedName),]
nisp_Pero <- nisp[grep("Peromyscus sp", nisp$RevisedName),]
percent_Pc <- nisp_Pc[,4:7]/(nisp_Pc[,4:7]+nisp_Pm[,4:7])
percent_Pm <- nisp_Pm[,4:7]/(nisp_Pc[,4:7]+nisp_Pm[,4:7])

nisp_Pero_Pc <- percent_Pc*nisp_Pero[,4:7]
nisp_Pero_Pm <- percent_Pm*nisp_Pero[,4:7]
  
nisp_Pc_total <- nisp_Pc[,4:7]+nisp_Pero_Pc
nisp_Pm_total <- nisp_Pm[,4:7]+nisp_Pero_Pm

# replace original values with new
nisp_revised[grep("Peromyscus californicus", nisp_revised$RevisedName),4:7] <- nisp_Pc_total
nisp_revised[grep("Peromyscus maniculatus", nisp_revised$RevisedName),4:7] <- nisp_Pm_total

# Sylvilagus
# calculate ratio, audubonii/bachmani:
nisp_Sa <- nisp[grep("Sylvilagus audubonii", nisp$RevisedName),]
nisp_Sb <- nisp[grep("Sylvilagus bachmani", nisp$RevisedName),]
nisp_Syl <- nisp[grep("Sylvilagus sp", nisp$RevisedName),]
percent_Sa <- nisp_Sa[,4:7]/(nisp_Sa[,4:7]+nisp_Sb[,4:7])
percent_Sb <- nisp_Sb[,4:7]/(nisp_Sa[,4:7]+nisp_Sb[,4:7])

nisp_Syl_Sa <- percent_Sa*nisp_Syl[,4:7]
nisp_Syl_Sb <- percent_Sb*nisp_Syl[,4:7]

nisp_Sa_total <- nisp_Sa[,4:7]+nisp_Syl_Sa
nisp_Sb_total <- nisp_Sb[,4:7]+nisp_Syl_Sb

# replace original values with new
nisp_revised[grep("Sylvilagus audubonii", nisp_revised$RevisedName),4:7] <- nisp_Sa_total
nisp_revised[grep("Sylvilagus bachmani", nisp_revised$RevisedName),4:7] <- nisp_Sb_total

# calculate community climates ----
boxes <- colnames(nisp_revised[,4:7])


for (i in 1:length(boxes)){
  focalBox <- boxes[i]
  colNumber <- match(focalBox, colnames(nisp_revised))

  # find NatureServeName from Revised Name
  speciesListIndexes <- which(!is.na(nisp_revised[,colNumber]))

  rangeMapNames <- matches$NatureServeName[match(nisp_revised$RevisedName[speciesListIndexes], matches$RevisedName)]
  
  # find names in nature serve species name list
  unique(match(rangeMapNames, NatureServeNameList))
  
  # match to species_clim_stats and calculate community climate - PA----
  
  # pull out medians, mins, and maxes based on presence-absence data
  medians <- do.call('rbind', lapply(species_clim_stats,function(x) x[which(row.names(x)=="Median"),]))[unique(match(rangeMapNames, NatureServeNameList)),]
  rownames(medians) <- NatureServeNameList[unique(match(rangeMapNames, NatureServeNameList))]
  community.mean.medians <- apply(medians, 2, mean)
  
  mins <- do.call('rbind', lapply(species_clim_stats,function(x) x[which(row.names(x)=="Min."),]))[unique(match(rangeMapNames, NatureServeNameList)),]
  rownames(mins) <- NatureServeNameList[unique(match(rangeMapNames, NatureServeNameList))]
  community.mean.mins <- apply(mins, 2, mean)
  
  maxes <- do.call('rbind', lapply(species_clim_stats,function(x) x[which(row.names(x)=="Max."),]))[unique(match(rangeMapNames, NatureServeNameList)),]
  rownames(maxes) <- NatureServeNameList[unique(match(rangeMapNames, NatureServeNameList))]
  community.mean.maxes <- apply(maxes, 2, mean)

# match to species_clim_stats and calculate community climate - abundance weighted----
  
}