# load libraries ####
library(readr)
library(raster)
library(dplyr)
library(maptools)
library(rgeos)
library(spThin)

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

#clip layers by full extent 
North_AM <- wrld_simpl[c(which(wrld_simpl@data$NAME=="United States"),
                         which(wrld_simpl@data$NAME=="Canada"), 
                         which(wrld_simpl@data$NAME=="Mexico")),]
NA_layers <- crop(WC, North_AM)
plot(NA_layers[[1]], xlim=c(-180, 0))

# climate niche model analyses ----

# set species names
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