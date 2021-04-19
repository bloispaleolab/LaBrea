# New climate niche model analyses ----

# set species names

# download specimens

# load in peromyscus maniculatus lineage 1 data 
pero_man <- read.csv("GIS/PM_all_coords.csv", header=T)

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
plot(WC, 1)
points(pero_man$LONG, pero_man$LAT, pch=16, col="red", cex=0.5)

