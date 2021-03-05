#### Peromyscus maniculatus modeling ##########
## Code developed by Robert Boria, modified by Jessica Blois #
options(java.parameters = "-Xmx8000m")
library("raster")
library("maptools")
library("spThin")
library("dismo")
library("rJava")
library("ENMeval")
library("rgeos")
library("rgdal")
library("rgbif")

## Find species occurrences and coordinates ####
# Note: still need to add this code to directly query gbif


# load in peromyscus maniculatus lineage 1 data 
pero_man <- read.csv("GIS/PM_all_coords.csv", header=T)

#load evironemntal variables 
#setwd("GIS/env/")
env_w <- list.files(path="../../../../Downloads/wc2/", pattern='tif', full.names=TRUE)
WC <- stack(env_w)
projection(WC) <- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
WC
plot(WC, 1)
points(pero_man$LONG, pero_man$LAT, pch=16, col="red", cex=0.5)

#load in linages shapefile 
# PM_lineages <-readOGR(dsn = "GIS/lineages/", layer = "lineages")

#plot(PM_lineages, add=T, col="black")

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

PM_all <- pero_man

#Clip WC to NA 
data(wrld_simpl)

# plot the data

North_AM <- wrld_simpl[c(which(wrld_simpl@data$NAME=="United States"),
                         which(wrld_simpl@data$NAME=="Canada"), 
                         which(wrld_simpl@data$NAME=="Mexico")),]

#clip layers by full extent 
NA_layers <- crop(WC, North_AM)
plot(NA_layers[[1]])

# Generate MCP for each species   
PM_mcp <- simpleMCP(PM_all[3:2])
plot(PM_mcp, add=T)
# Generate a buffered MCP (here 3.0 degrees)
PM_buf <- gBuffer(PM_mcp, width=3.0)
plot(PM_buf, add=T)

#clip layers by full extent 
layers <- crop(NA_layers, PM_buf)
PM_pres <- mask(layers, PM_buf)

plot(PM_pres[[1]])
#plot(north_am, 1)
points(pero_man$LONG, pero_man$LAT, pch=16, col="red", cex=0.5)


#spatial filter in spThin
PM_thin <- thin(loc.data = PM_all, 
                lat.col = "LAT", long.col = "LONG", 
                spec.col = "id", 
                thin.par = 50, reps = 100, 
                locs.thinned.list.return = TRUE, 
                write.files = TRUE, 
                max.files = 5,
                out.dir = "GIS/thin_50/", out.base = "PM_f", 
                write.log.file = TRUE,
                log.file = "GIS/thin_50/PM_thin_log.txt")

#load in spatial thin dataset
PM_F <- read.csv("GIS/thin_50/PM_f_thin1.csv")

#plot 
plot(PM_pres[[1]])
points(pero_man$LONG, pero_man$LAT, pch=16, col="black", cex=0.5)
points(PM_F$LONG, PM_F$LAT, pch=16, col="red", cex=.5)

#lineage one model testing 
tune.args <- list(fc = c("L", "LQ", "H", "LQH", 'LQHT'), rm = 1:6)
model_sel_50 <- ENMevaluate(PM_F[2:3], PM_pres, algorithm = "maxent.jar",
                            tune.args = tune.args, partitions = "block", overlap=F, 
                            updateProgress = F, numCores = 4)
## look at the results 
PM_result <- model_sel@results

write.csv(PM_result, "../output/Dissertation/ENMs/PM_results.csv", quote=FALSE, row.names=FALSE)

#see below
#PM_model Thinning parameter 25km  
#M_predict
#PM_model_50 Thinning parameter 50km  
#M_predict_50 
##Generating final models for P. maniculatus
## Using the top perfroming model (for now)  ##
#projected to SR; #PM_ best settings were linear Quadratic; RM = 6 
args=c("noaddsamplestobackground","noautofeature", "noproduct","nohinge","nothreshold","betamultiplier=6.0", "responsecurves=true") 
pred.args_N <- c("outputformat=Cloglog", "doclamp=TRUE")

PM_model_50 <- maxent(PM_pres, PM_F[2:3], args = args)
PM_predict_50 <- predict(PM_pres, PM_model_50, args=pred.args_N) 
plot(PM_predict_50)
points(PM_F$LONG, PM_F$LAT, pch=16, col="black", cex=.25)
plot(PM_predict) 
points(PM_F$LONG, PM_F$LAT, pch=16, col="red", cex=.25)


plot(wrld_simpl[c(which(wrld_simpl@data$NAME=="United States"), #plot the outline of US, CA, and MEX
                  which(wrld_simpl@data$NAME=="Canada"), 
                  which(wrld_simpl@data$NAME=="Mexico")),], 
     xlim=c(-180, 0), #this controls the x (longitude) boundaries
     cex=0.25)   #this controls how large the points are (ranges from 0 to 1)
plot(PM_predict, add=T) 
points(PM_F$LONG, PM_F$LAT, pch=16, col="black", cex=.25)
