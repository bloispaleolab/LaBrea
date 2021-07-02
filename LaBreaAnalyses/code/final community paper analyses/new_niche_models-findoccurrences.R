# species occurrence data ----
# retrieve and clean
# install.packages(c("ggmap", "mapdata"))
# library(wallace)
#install.packages("tidyverse")
library(spocc)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(maptools)
library(readr)
library(dplyr)

speciesName <- c('Mustela frenata',
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

# query selected database for occurrence records
for (i in 1:length(speciesName)){
  results <- spocc::occ(query = speciesName[i], from = "gbif", limit = 10000, has_coords = TRUE, gbifopts = list(basisOfRecord = 'PRESERVED_SPECIMEN'))
  # retrieve data table from spocc object
  results.data <- results[["gbif"]]$data[[paste(strsplit(speciesName[i], split=' ')[[1]], collapse='_')]]
  # remove rows with duplicate coordinates
  occs.dups <- duplicated(results.data[c('longitude', 'latitude')])
  occs <- results.data[!occs.dups,]
  # make sure latitude and longitude are numeric (sometimes they are characters)
  occs$latitude <- as.numeric(occs$latitude)
  occs$longitude <- as.numeric(occs$longitude)
  
  # remove networkKeys
  occs <-  select(occs, -networkKeys)
  
  
  # remove the lat/longs at 0,0
  if(any(c(which(occs$longitude==0), which(occs$latitude==0)))){
    occs <- occs[-unique(which(occs$longitude==0), which(occs$latitude==0)),]
  }

  # remove points on the Hawaiian Islands
  if(any(unique(occs$islandGroup[!is.na(occs$islandGroup)]) == "Hawaiian Islands")){
    occs <- occs[-which(occs$islandGroup == "Hawaiian Islands"),]
  }
  
  # give all records a unique ID
  occs$occID <- row.names(occs)
  
  # save the occurrences for later point cleanup
  write_delim(occs, file=paste0("data/gbif_occurrences/", speciesName[i], "_initialOccurrences.txt"), delim="\t")

}


# more intense point cleaning
data(wrld_simpl)
pdf(file="data/gbif_occurrences/final_niche_occurrences.pdf")
for (i in 1:length(speciesName)){

  occs <- read_delim(paste0("data/gbif_occurrences/", speciesName[i], "_initialOccurrences.txt"), delim="\t")

  # plot(wrld_simpl, axes=TRUE, col="light yellow", main=paste0(speciesName[i], "-initial"))
  # points(occs$longitude, occs$latitude, pch=16, col="red")

  # Mustela frenata - DONE (i=1)
  if (speciesName[i] == "Mustela frenata"){
    occs <- occs[-which(occs$latitude < -40),]
    occs <- occs[-which(occs$longitude > -20),]
    xlim =c(-130, -50)
    ylim= c(-20,80)
  }
  # Mephitis mephitis - DONE (i=2)
  if (speciesName[i] == "Mephitis mephitis"){
    occs <- occs[-which(occs$latitude < 20),]
    occs <- occs[-which(occs$longitude > -20),]
    xlim =c(-130, -55)
    ylim= c(20,60)
  }
  # Spilogale gracilis - DONE (i=3)
  if (speciesName[i] == "Spilogale gracilis"){
    occs <- occs[-which(occs$longitude > 0),]
    xlim =c(-130, -80)
    ylim= c(20,60)
  }
  # Microtus californicus - DONE (i=4)
    # removed points at tip of baja and in mainland mexico and US that are too far east.
    # a few points seem outside the distribution, but better to accept that uncertainty
  if (speciesName[i] == "Microtus californicus"){
    occs <- occs[-which(occs$longitude > -113),]
    xlim =c(-130, -110)
    ylim= c(20,50)
  }
  # Peromyscus californicus - DONE (i=5) 
    # a few points seem outside the distribution, but better to accept that uncertainty
  if (speciesName[i] == "Peromyscus californicus"){
    occs <- occs[-which(occs$latitude > 45),]
    occs <- occs[-which(occs$longitude > -10),]
    xlim =c(-130, -110)
    ylim= c(20,50)
  }
  # Peromyscus maniculatus - DONE (i=6)
  if (speciesName[i] == "Peromyscus californicus"){
    xlim=c(-140, -60)
    ylim= c(20,60)
  }
  # Neotoma macrotis - DONE (i=7)
  if (speciesName[i] == "Neotoma macrotis"){
    xlim=c(-130, -100)
    ylim= c(20,50)
  }
  # Reithrodontomys megalotis - DONE (i=8)
  if (speciesName[i] == "Reithrodontomys megalotis"){
    occs <- occs[-which(occs$longitude > -20),]
    xlim=c(-140, -70)
    ylim= c(10,60)
  }
  # Onychomys torridus - DONE (i=9)
    # based on IUCN, "native to Mexico and the states of Arizona, California, Nevada, New Mexico, and Utah in the United States."
    # so removed the occurrences in states below, but kept OR, CO, TX because flanking states and so boundaries uncertain
  if (speciesName[i] == "Onychomys torridus"){
    occs <- occs[-which(occs$latitude < 10),]
    occs <- occs[-which(occs$stateProvince == "New York" | occs$stateProvince == "Minnesota" | occs$stateProvince == "Montana"| is.na(occs$stateProvince)),]
    xlim=c(-140, -60)
    ylim= c(10,60)
  }
  # Sylvilagus audubonii- DONE (i=10)
  if (speciesName[i] == "Sylvilagus audubonii"){
    occs <- occs[-which(occs$longitude > -90),]
    xlim=c(-140, -90)
    ylim= c(10,60)
  }
  # Sylvilagus bachmani - DONE (i=11)
  if (speciesName[i] == "Sylvilagus bachmani"){
    occs <- occs[-which(occs$longitude > -100),]
    xlim=c(-130, -100)
    ylim= c(15,55)
  }
  # Otospermophilus beecheyi - DONE (i=12)
  if (speciesName[i] == "Otospermophilus beecheyi"){
    occs <- occs[-which(occs$longitude > -110),]
    occs <- occs[-which(occs$longitude < -150),]
    xlim=c(-130, -100)
    ylim= c(20,60)
  }
  # Neotamias merriami - DONE (i=13)
  if (speciesName[i] == "Neotamias merriami"){
    xlim=c(-130, -100)
    ylim= c(25,45)
  }
  # Dipodomys agilis - DONE (i=14)
  if (speciesName[i] == "Dipodomys agilis"){
    occs <- occs[-which(occs$latitude > 40),]
    occs <- occs[-which(occs$longitude > -110),]
    xlim= c(-130, -100)
    ylim= c(20,50)
  }
  # Chaetodipus californicus - DONE (i=15)
  if (speciesName[i] == "Chaetodipus californicus"){
    occs <- occs[-which(occs$longitude > -100),]
    xlim=c(-130, -110)
    ylim= c(25,45)
  }
  # Thomomys bottae - DONE (i=16)
  if (speciesName[i] == "Thomomys bottae"){
    occs <- occs[-which(occs$longitude > -90),]
    xlim=c(-130, -90)
    ylim= c(20,50)
    }
  # Scapanus latimanus - DONE (i=17) 
  if (speciesName[i] == "Scapanus latimanus"){
    occs <- occs[-which(occs$longitude > -113),]
    xlim=c(-130, -110)
    ylim= c(20,50)
  }
  # Sorex ornatus - DONE (i=18)
  if (speciesName[i] == "Sorex ornatus"){
    xlim=c(-130, -100)
    ylim= c(20,50)
    #note: occurrence at tip of Baja is ok
  }
  
  plot(wrld_simpl, xlim=xlim, ylim=ylim, axes=TRUE, col="light yellow", main=paste0(speciesName[i], "-final"))
  points(occs$longitude, occs$latitude, pch=16, col="red")
  
write_delim(occs, file=paste0("data/gbif_occurrences/", speciesName[i], "_finalOccurrences.txt"), delim="\t")
}
dev.off()