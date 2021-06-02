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

for (i in 1:length(speciesName)){

  occs <- read_delim(paste0("data/gbif_occurrences/", speciesName[i], "_initialOccurrences.txt"), delim="\t")

  plot(wrld_simpl, xlim=c(-180,-40), ylim=c(0,80), axes=TRUE, col="light yellow", main=paste0(speciesName[i], "-initial"))
  points(occs$longitude, occs$latitude, pch=16, col="red")
  
  # Mustela frenata
  if (speciesName[i] == "Mustela frenata"){
    occs <- occs[-which(occs$latitude < -40),]
  }
  # Mephitis mephitis
  if (speciesName[i] == "Mephitis mephitis"){
    occs <- occs[-which(occs$latitude < 0),]
  }
  # Microtus californicus
  if (speciesName[i] == "Microtus californicus"){
    occs <- occs[-which(occs$longitude > -50),]
  }
  # Peromyscus californicus
  if (speciesName[i] == "Peromyscus californicus"){
    occs <- occs[-which(occs$latitude > 45),]
  }

  # Onychomys torridus
  if (speciesName[i] == "Onychomys torridus"){
    occs <- occs[-which(occs$latitude < 10),]
  }
  
  # Sylvilagus audubonii
  if (speciesName[i] == "Sylvilagus audubonii"){
    occs <- occs[-which(occs$longitude > -90),]
  }

  # Sylvilagus bachmani
  if (speciesName[i] == "Sylvilagus bachmani"){
    occs <- occs[-which(occs$longitude > -100),]
  }

  # Otospermophilus beecheyi
  if (speciesName[i] == "Otospermophilus beecheyi"){
    occs <- occs[-which(occs$longitude > -110),]
    occs <- occs[-which(occs$longitude < -150),]
  }
  
  # Dipodomys agilis
  if (speciesName[i] == "Dipodomys agilis"){
    occs <- occs[-which(occs$latitude < 20),]
  }
  
  # Chaetodipus californicus
  if (speciesName[i] == "Chaetodipus californicus"){
    occs <- occs[-which(occs$latitude < 20),]
  }
  
  # Thomomys bottae
  if (speciesName[i] == "Thomomys bottae"){
    occs <- occs[-which(occs$longitude > -90),]
  }
  
  # Sorex ornatus
  if (speciesName[i] == "Sorex ornatus"){
    #occs[which(occs$latitude < 25),]$'issues'
    #occs <- occs[-which(occs$latitude < 25),'issues'] #maybe these are ok?
  }
  
    plot(wrld_simpl, xlim=c(-180,-40), ylim=c(0,80), axes=TRUE, col="light yellow", main=paste0(speciesName[i], "-final"))
  points(occs$longitude, occs$latitude, pch=16, col="red")
  
}
