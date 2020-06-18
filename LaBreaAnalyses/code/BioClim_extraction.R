#install.packages("remotes")
#remotes::install_github("rebeccalpowell/grassmapr")1
library(grassmapr)
library(remotes)
library(raster)
library(sp)
library(rgdal)
library(rgeos)
library(plyr)

## BIOCLIM ####
# Read in the relevant bioclim raster files 
vars <- c(1, 4:6, 12:15)
var_names <- c("temp_mean", "temp_SD", "temp_max", "temp_min", "precip_total_annual", "precip_monthly_max", "precip_monthly_min", "precip_monthly_SD")

bioclim_files_all <- list.files("data/raw/BioClim/")
bioclim_files <- bioclim_files_all[match(paste0("wc2.1_5m_bio_", vars, ".tif"), bioclim_files_all)]

bioclim_stack <- stack(paste0("data/raw/BioClim/", bioclim_files))

## Read in the species range maps, overlay on bioclim, and extract variables, calculate species range statistics ####

# Read in master names file

lookup <- read.delim("data/processed/Master_trait_taxon_names.txt", header=T, sep="\t")

# generate list of unique species to match
uniqueSp <- levels(lookup$Species)
  
for (i in 1:length(uniqueSp)){
  rowForMatching <- which(lookup$Species == uniqueSp[i])[1]
  
  if (!is.na(lookup$NatureServeID[rowForMatching])){
    # generate blank file to store species-range statistics
    range_summary <- as.data.frame(matrix(ncol = length(bioclim_files), nrow=6))
    rownames(range_summary) <- c("mean", "median", "var", "sd", "max", "min")
    
    # Read in the range shapefiles and simplify
    range <- readOGR(dsn=path.expand(paste0("data/raw/SpeciesRangeShapefiles/", 
                                            lookup$NatureServeID[rowForMatching], ".shp")))
    range <- gSimplify(range, tol=0.01, topologyPreserve=FALSE)
    #plot(range, add=T)
  
    # extract bioclim variables from the range
    bioclim_data <- extract(bioclim_stack, range)
    bioclim_data <- ldply(bioclim_data, rbind)
    # calculate various summary statistics on the range-wide climate data
    range_summary[1,] <- apply(bioclim_data, 2, mean, na.rm=TRUE)#Added na.rm=TRUE, this is what messed up M. frenata
    range_summary[2,] <- apply(bioclim_data, 2, median, na.rm=TRUE)
    range_summary[3,] <- apply(bioclim_data, 2, var, na.rm=TRUE)
    range_summary[4,] <- apply(bioclim_data, 2, sd, na.rm=TRUE)
    range_summary[5,] <- apply(bioclim_data, 2, max, na.rm=TRUE)
    range_summary[6,] <- apply(bioclim_data, 2, min, na.rm=TRUE)
    colnames(range_summary) <- colnames(bioclim_data) 
    
    # save the file
    write.csv(range_summary, file=paste0("data/processed/range_wide_bioclim_stats/range_wide_bioclim_stats_", uniqueSp[i], ".csv"), row.names=T)
  
    }
}


# ndvi ----
# Calculate mean NDVI across all months in a year and save for future use
ndvi_mean <- calc(ndvi_NA, mean)
writeRaster(ndvi_mean, filename="ndvi.asc", format = "ascii", overwrite=TRUE)

## Read in the species range maps, overlay on ndvi, and extract variables, calculate species range statistics ####

# Read in master names file

lookup <- read.delim("data/processed/Master_trait_taxon_names.txt", header=T, sep="\t")

# generate list of unique species to match
uniqueSp <- levels(lookup$Species)

for (i in 1:length(uniqueSp)){
  rowForMatching <- which(lookup$Species == uniqueSp[i])[1]
  
  if (!is.na(lookup$NatureServeID[rowForMatching])){
    # generate blank file to store species-range statistics
    range_summary <- as.data.frame(matrix(ncol = 1, nrow=6))
    rownames(range_summary) <- c("mean", "median", "var", "sd", "max", "min")
    
    # Read in the range shapefiles and simplify
    range <- readOGR(dsn=path.expand(paste0("data/raw/SpeciesRangeShapefiles/", 
                                            lookup$NatureServeID[rowForMatching], ".shp")))
    range <- gSimplify(range, tol=0.01, topologyPreserve=FALSE)
    #plot(range, add=T)
    
    # extract ndvi data from the range
    ndvi_data <- extract(ndvi_mean, range)
    ndvi_data <- unlist(ndvi_data)
    # calculate various summary statistics on the range-wide climate data
    range_summary[1,] <- mean(ndvi_data, na.rm=TRUE)#Added na.rm=TRUE, this is what messed up M. frenata
    range_summary[2,] <- median(ndvi_data, na.rm=TRUE)
    range_summary[3,] <- var(ndvi_data, na.rm=TRUE)
    range_summary[4,] <- sd(ndvi_data, na.rm=TRUE)
    range_summary[5,] <- max(ndvi_data, na.rm=TRUE)
    range_summary[6,] <- min(ndvi_data, na.rm=TRUE)
    colnames(range_summary) <- "ndvi"
    
    # save the file
    write.csv(range_summary, file=paste0("data/processed/range_wide_ndvi_stats/range_wide_ndvi_stats_", uniqueSp[i], ".csv"), row.names=T)
    
  }
}