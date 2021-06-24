# load libraries ####
library(readr)
library(raster)
library(dplyr)
library(maptools)
library(rgeos)
library(spThin)
library(rgdal)

# load special functions ####
## function to create a minimum convex polygon 
# simpleMCP written by Jamie Kass spring 2014 
# Rob Boria vetted code, spring 2014, during his second chapter analyses by comparing code to MCP created in ArcGIS
simpleMCP <- function (xy) {
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2])
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  p <- SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1)))
}

# set species names ####
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


# arrange worldclim climate data ####

#load climate variables 
env_w <- list.files(path="../../GIS_data/wc2.1_2.5m_bio/", pattern='tif', full.names=TRUE)
wc_names <- gsub("../../GIS_data/wc2.1_2.5m_bio//wc2.1_2.5m_", "", env_w)
wc_names <- gsub(".tif", "", wc_names)
wc_names <- wc_names[c(1,12:19, 2:11)] #reorder variables
env_w <- env_w[c(1,12:19, 2:11)]
WC <- stack(env_w)
projection(WC) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
names(WC) <- wc_names
WC

#Clip WC to NA 
data(wrld_simpl)

#load North America shapefile
# NorthAmerica <- readOGR(dsn="../../GIS_data/stanford-cq068zf3261-shapefile-North_America", layer = "cq068zf3261")

#clip layers by full extent 
North_AM <- wrld_simpl[c(which(wrld_simpl@data$NAME=="United States"),
                         which(wrld_simpl@data$NAME=="Canada"), 
                         which(wrld_simpl@data$NAME=="Mexico")),]

# NA_layers <- crop(WC, NorthAmerica)
- 
  NA_layers <- crop(WC, North_AM)
plot(NA_layers[[1]], xlim=c(-180, 0))

North_Am2 <- crop(WC, extent(-180, -40, 0, 90))
plot(North_Am2[[1]])

writeRaster(North_Am2, file="data/processed/NA_climate_layers.tif", format="GTiff", overwrite=TRUE)

# spatial thinning ----

for (i in 1:length(speciesName)){
  
  # load in occurrences
  occs <- read_delim(file=paste0("data/gbif_occurrences/", speciesName[i], "_finalOccurrences.txt"), delim="\t")
  occs2 <- occs %>% select(occID, name, longitude, latitude)
  
  # check point plotting
  plot(WC, 1, main=speciesName[i])
  points(occs2$longitude, occs2$latitude, pch=16, col="red", cex=0.5)
  
  # Generate MCP for each species   
  occs_mcp <- simpleMCP(occs2[3:4])
  plot(occs_mcp, add=T)
  # Generate a buffered MCP (here 3.0 degrees)
  occs_buf <- gBuffer(occs_mcp, width=3.0)
  plot(occs_buf, add=T)
  
  #clip layers by full extent 
  layers <- crop(NA_layers, occs_buf)
  occs_pres <- mask(layers, occs_buf)
  
  #spatial filter in spThin
  occs_thin <- thin(loc.data = occs2, 
                    lat.col = "latitude", long.col = "longitude", 
                    spec.col = "occID", 
                    thin.par = 50, reps = 100, 
                    locs.thinned.list.return = TRUE, 
                    write.files = TRUE, 
                    max.files = 5,
                    out.dir = "data/gbif_occurrences/thin_50/", out.base = paste0(speciesName[i], "_t"), 
                    write.log.file = TRUE,
                    log.file = paste0("data/gbif_occurrences/thin_50/", speciesName[i],"_thin_log.txt"))
  
  #load in spatial thin dataset
  occs_thin_final <- read.csv(paste0("data/gbif_occurrences/thin_50/", speciesName[i], "_t_thin1.csv"))
  
  #plot 
  plot(occs_pres[[1]])
  points(occs2$longitude, occs2$latitude, pch=16, col="black", cex=0.5)
  points(occs_thin_final$longitude, occs_thin_final$latitude, pch=16, col="red", cex=.5)
  
}

# Associate species occurrences with climate variables ----
# read in climate layers
NA_layers <- stack("data/processed/NA_climate_layers.tif")
species_clim_stats <- list()
species_occs_final <- list()
species_clim_final <- list()

for (i in 1:length(speciesName)){
  
  #load in spatial thin dataset
  occs_thin_final <- read.csv(paste0("data/gbif_occurrences/thin_50/", speciesName[i], "_t_thin1.csv"))
  
  # Generate MCP for the species   
  occs_mcp <- simpleMCP(occs_thin_final[2:3])

  # plot(NA_layers[[1]])
  # plot(occs_mcp, add=T)
  # points(occs_thin_final$longitude, occs_thin_final$latitude, pch=16, col="black", cex=.5)
  
  # extract climate at the points
  clim_occs <- extract(NA_layers, occs_thin_final[, c('longitude', 'latitude')])
  
  # points to be removed:
    # misc other points in South America
  SA_points <- intersect(which(occs_thin_final$longitude >-80), which(occs_thin_final$latitude <15))
    # Add clim NA values
  remove <- union(which(is.na(clim_occs[,1])), SA_points)

  if (length(remove)>0){
   clim_occs_trim <- clim_occs[-remove,]
  # remove NA values from original occurrences
  occs_thin_final_NA <- occs_thin_final[-remove,]
  # Check for NA values. Should be FALSE
  any(is.na(clim_occs_trim))
  }else{
    clim_occs_trim <- clim_occs
    occs_thin_final_NA <- occs_thin_final
  }
  
  # check points 
  plot(NA_layers[[1]])
  plot(occs_mcp, add=T)
  points(occs_thin_final_NA$longitude, occs_thin_final_NA$latitude, pch=16, col="black", cex=.5)
  
 # calculate climate statistics
  species_clim_stats[[i]] <- apply(clim_occs_trim, 2, summary)

  # save final occurrences
  species_occs_final[[i]] <- occs_thin_final_NA
  
  # save final climate dataframe
  species_clim_final[[i]] <- clim_occs_trim
  
}

save(species_clim_stats, file="data/gbif_occurrences/final/species_clim_stats.RData")
save(species_occs_final, file="data/gbif_occurrences/final/species_occs_final.RData")
save(species_clim_final, file="data/gbif_occurrences/final/species_clim_final.RData")