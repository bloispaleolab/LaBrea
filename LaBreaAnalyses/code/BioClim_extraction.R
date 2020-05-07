#install.packages("remotes")
#remotes::install_github("rebeccalpowell/grassmapr")1
library(grassmapr)
library(remotes)
library(raster)
library(sp)
library(rgdal)
library(rgeos)


temp <- raster("data/raw/BioClim/wc2.1_5m_bio_1.tif")#mean temp
plot(temp)

temp <- raster("data/raw/BioClim/wc2.1_5m_bio_4.tif") #temp SD
plot(temp)

temp <- raster("data/raw/BioClim/wc2.1_5m_bio_5.tif") #temp max
plot(temp)

temp <- raster("data/raw/BioClim/wc2.1_5m_bio_6.tif") #temp min
plot(temp)

precip <- raster("data/raw/BioClim/wc2.1_5m_bio_12.tif") # annual precip
plot(precip)

precip <- raster("data/raw/BioClim/wc2.1_5m_bio_15.tif") #precip SD
plot(precip)

precip <- raster("data/raw/BioClim/wc2.1_5m_bio_13.tif") #precip max
plot(precip)

precip <- raster("data/raw/BioClim/wc2.1_5m_bio_14.tif") #precip min
plot(precip)


#Microtus californicus
Mcali_range <- readOGR(dsn=path.expand("data/raw/SpeciesRangeShapefiles/micr_cali_pl.shp"))
Mcali_range <- gSimplify(Mcali_range, tol=0.01, topologyPreserve=FALSE)
plot(Mcali_range, add=T)

#temp
Mcali_temp <- extract(temp, Mcali_range)

Mcali_temp_all <- unlist(Mcali_temp)
Mcali_temp_all
mean(Mcali_temp_all)#na.rm=TRUE
max(Mcali_temp_all)
min(Mcali_temp_all)

#precip
Mcali_precip <- extract(precip, Mcali_range)

Mcali_precip_all <- unlist(Mcali_precip)
Mcali_precip_all
mean(Mcali_precip_all)#na.rm=TRUE
max(Mcali_precip_all)
min(Mcali_precip_all)

#Peromyscus californicus
Pcali_range <- readOGR(dsn=path.expand("data/raw/SpeciesRangeShapefiles/pero_cali_pl.shp"))
Pcali_range <- gSimplify(Pcali_range, tol=0.01, topologyPreserve=FALSE)
plot(Pcali_range, add=T)

Pcali_temp <- extract(temp, Pcali_range)

Pcali_temp_all <- unlist(Pcali_temp)
Pcali_temp_all
mean(Pcali_temp_all)#na.rm=TRUE
max(Pcali_temp_all)
min(Pcali_temp_all)

#precip
Pcali_precip <- extract(precip, Pcali_range)

Pcali_precip_all <- unlist(Pcali_precip)
Pcali_precip_all
mean(Pcali_precip_all)#na.rm=TRUE
max(Pcali_precip_all)
min(Pcali_precip_all)


#Neotoma macrotis
Nmacr_range <- readOGR(dsn=path.expand("data/raw/SpeciesRangeShapefiles/neot_macr_pl.shp"))
Nmacr_range <- gSimplify(Nmacr_range, tol=0.01, topologyPreserve=FALSE)
plot(Nmacr_range, add=T)

Nmacr_temp <- extract(temp, Nmacr_range)

Nmacr_temp_all <- unlist(Nmacr_temp)
Nmacr_temp_all
mean(Nmacr_temp_all)#na.rm=TRUE
max(Nmacr_temp_all)
min(Nmacr_temp_all)

#precip
Nmacr_precip <- extract(precip, Nmacr_range)

Nmacr_precip_all <- unlist(Nmacr_precip)
Nmacr_precip_all
mean(Nmacr_precip_all)#na.rm=TRUE
max(Nmacr_precip_all)
min(Nmacr_precip_all)

#Reithrodontomys megalotis
Rmega_range <- readOGR(dsn=path.expand("data/raw/SpeciesRangeShapefiles/reit_mega_pl.shp"))
Rmega_range <- gSimplify(Rmega_range, tol=0.01, topologyPreserve=FALSE)
plot(Rmega_range, add=T)

Rmega_temp <- extract(temp, Rmega_range)

Rmega_temp_all <- unlist(Rmega_temp)
Rmega_temp_all
mean(Rmega_temp_all)#na.rm=TRUE
max(Rmega_temp_all)
min(Rmega_temp_all)
#Precip
Rmega_precip <- extract(precip, Rmega_range)

Rmega_precip_all <- unlist(Rmega_precip)
Rmega_precip_all
mean(Rmega_precip_all)#na.rm=TRUE
max(Rmega_precip_all)
min(Rmega_precip_all)

#Onychomys torridus
Otorr_range <- readOGR(dsn=path.expand("data/raw/SpeciesRangeShapefiles/Onyc_torr_pl.shp"))
Otorr_range <- gSimplify(Otorr_range, tol=0.01, topologyPreserve=FALSE)
plot(Otorr_range, add=T)

Otorr_temp <- extract(temp, Otorr_range)

Otorr_temp_all <- unlist(Otorr_temp)
Otorr_temp_all
mean(Otorr_temp_all)#na.rm=TRUE
max(Otorr_temp_all)
min(Otorr_temp_all)
#Precip
Otorr_precip <- extract(precip, Otorr_range)

Otorr_precip_all <- unlist(Otorr_precip)
Otorr_precip_all
mean(Otorr_precip_all)#na.rm=TRUE
max(Otorr_precip_all)
min(Otorr_precip_all)

#Sylvilagus audubonii
Saudu_range <- readOGR(dsn=path.expand("data/raw/SpeciesRangeShapefiles/Sylv_audu_pl.shp"))
Saudu_range <- gSimplify(Saudu_range, tol=0.01, topologyPreserve=FALSE)
plot(Saudu_range, add=T)

Saudu_temp <- extract(temp, Saudu_range)

Saudu_temp_all <- unlist(Saudu_temp)
Saudu_temp_all
mean(Saudu_temp_all)#na.rm=TRUE
max(Saudu_temp_all)
min(Saudu_temp_all)
#Precip
Saudu_precip <- extract(precip, Saudu_range)

Saudu_precip_all <- unlist(Saudu_precip)
Saudu_precip_all
mean(Saudu_precip_all)#na.rm=TRUE
max(Saudu_precip_all)
min(Saudu_precip_all)

#Sylvilagus bachmani
Sbach_range <- readOGR(dsn=path.expand("data/raw/SpeciesRangeShapefiles/Sylv_bach_pl.shp"))
Sbach_range <- gSimplify(Sbach_range, tol=0.01, topologyPreserve=FALSE)
plot(Sbach_range, add=T)

Sbach_temp <- extract(temp, Sbach_range)

Sbach_temp_all <- unlist(Sbach_temp)
Sbach_temp_all
mean(Sbach_temp_all)#na.rm=TRUE
max(Sbach_temp_all)
min(Sbach_temp_all)
#Precip
Sbach_precip <- extract(precip, Sbach_range)

Sbach_precip_all <- unlist(Sbach_precip)
Sbach_precip_all
mean(Sbach_precip_all)#na.rm=TRUE
max(Sbach_precip_all)
min(Sbach_precip_all)

#Otospermophilus beecheyi
Obeec_range <- readOGR(dsn=path.expand("data/raw/SpeciesRangeShapefiles/sper_beec_pl.shp"))
Obeec_range <- gSimplify(Obeec_range, tol=0.01, topologyPreserve=FALSE)
plot(Obeec_range, add=T)

Obeec_temp <- extract(temp, Obeec_range)

Obeec_temp_all <- unlist(Obeec_temp)
Obeec_temp_all
mean(Obeec_temp_all)#na.rm=TRUE
max(Obeec_temp_all)
min(Obeec_temp_all)
#Precip
Obeec_precip <- extract(precip, Obeec_range)

Obeec_precip_all <- unlist(Obeec_precip)
Obeec_precip_all
mean(Obeec_precip_all)#na.rm=TRUE
max(Obeec_precip_all)
min(Obeec_precip_all)

#Neotamias merriami
Nmerr_range <- readOGR(dsn=path.expand("data/raw/SpeciesRangeShapefiles/neot_merr_pl.shp"))
Nmerr_range <- gSimplify(Nmerr_range, tol=0.01, topologyPreserve=FALSE)
plot(Nmerr_range, add=T)

Nmerr_temp <- extract(temp, Nmerr_range)

Nmerr_temp_all <- unlist(Nmerr_temp)
Nmerr_temp_all
mean(Nmerr_temp_all)#na.rm=TRUE
max(Nmerr_temp_all)
min(Nmerr_temp_all)
#Precip
Nmerr_precip <- extract(precip, Nmerr_range)

Nmerr_precip_all <- unlist(Nmerr_precip)
Nmerr_precip_all
mean(Nmerr_precip_all)#na.rm=TRUE
max(Nmerr_precip_all)
min(Nmerr_precip_all)

#Dipodomys agilis
Dagil_range <- readOGR(dsn=path.expand("data/raw/SpeciesRangeShapefiles/Dipo_agil_pl.shp"))
Dagil_range <- gSimplify(Dagil_range, tol=0.01, topologyPreserve=FALSE)
plot(Dagil_range, add=T)

Dagil_temp <- extract(temp, Dagil_range)

Dagil_temp_all <- unlist(Dagil_temp)
Dagil_temp_all
mean(Dagil_temp_all)#na.rm=TRUE
max(Dagil_temp_all)
min(Dagil_temp_all)
#Precip
Dagil_precip <- extract(precip, Dagil_range)

Dagil_precip_all <- unlist(Dagil_precip)
Dagil_precip_all
mean(Dagil_precip_all)#na.rm=TRUE
max(Dagil_precip_all)
min(Dagil_precip_all)

#Chaeodipus californicus
Ccali_range <- readOGR(dsn=path.expand("data/raw/SpeciesRangeShapefiles/chae_cali_pl.shp"))
Ccali_range <- gSimplify(Ccali_range, tol=0.01, topologyPreserve=FALSE)
plot(Ccali_range, add=T)

Ccali_temp <- extract(temp, Ccali_range)

Ccali_temp_all <- unlist(Ccali_temp)
Ccali_temp_all
mean(Ccali_temp_all)#na.rm=TRUE
max(Ccali_temp_all)
min(Ccali_temp_all)
#Precip
Ccali_precip <- extract(precip, Ccali_range)

Ccali_precip_all <- unlist(Ccali_precip)
Ccali_precip_all
mean(Ccali_precip_all)#na.rm=TRUE
max(Ccali_precip_all)
min(Ccali_precip_all)

#Scapanus latimanus
Scap_range <- readOGR(dsn=path.expand("data/raw/SpeciesRangeShapefiles/scap_lati_pl.shp"))
Scap_range <- gSimplify(Scap_range, tol=0.01, topologyPreserve=FALSE)
plot(Scap_range, add=T)

Scap_temp <- extract(temp, Scap_range)

Scap_temp_all <- unlist(Scap_temp)
Scap_temp_all
mean(Scap_temp_all)#na.rm=TRUE
max(Scap_temp_all)
min(Scap_temp_all)
#Precip
Scap_precip <- extract(precip, Scap_range)

Scap_precip_all <- unlist(Scap_precip)
Scap_precip_all
mean(Scap_precip_all)#na.rm=TRUE
max(Scap_precip_all)
min(Scap_precip_all)

#Sorex ornatus
Sorna_range <- readOGR(dsn=path.expand("data/raw/SpeciesRangeShapefiles/sore_orna_pl.shp"))
Sorna_range <- gSimplify(Sorna_range, tol=0.01, topologyPreserve=FALSE)
plot(Sorna_range, add=T)

Sorna_temp <- extract(temp, Sorna_range)

Sorna_temp_all <- unlist(Sorna_temp)
Sorna_temp_all
mean(Sorna_temp_all)#na.rm=TRUE
max(Sorna_temp_all)
min(Sorna_temp_all)
#Precip
Sorna_precip <- extract(precip, Sorna_range)

Sorna_precip_all <- unlist(Sorna_precip)
Sorna_precip_all
mean(Sorna_precip_all)#na.rm=TRUE
max(Sorna_precip_all)
min(Sorna_precip_all)

#Thomomys bottae
Tbott_range <- readOGR(dsn=path.expand("data/raw/SpeciesRangeShapefiles/thom_bott_pl.shp"))
Tbott_range <- gSimplify(Tbott_range, tol=0.01, topologyPreserve=FALSE)
plot(Tbott_range, add=T)

Tbott_temp <- extract(temp, Tbott_range)

Tbott_temp_all <- unlist(Tbott_temp)
Tbott_temp_all
mean(Tbott_temp_all)#na.rm=TRUE
max(Tbott_temp_all)
min(Tbott_temp_all)
#Precip
Tbott_precip <- extract(precip, Tbott_range)

Tbott_precip_all <- unlist(Tbott_precip)
Tbott_precip_all
mean(Tbott_precip_all)#na.rm=TRUE
max(Tbott_precip_all)
min(Tbott_precip_all)

#Spilogale gracilis
Sgrac_range <- readOGR(dsn=path.expand("data/raw/SpeciesRangeShapefiles/spil_grac_pl.shp"))
Sgrac_range <- gSimplify(Sgrac_range, tol=0.01, topologyPreserve=FALSE)
plot(Sgrac_range, add=T)

Sgrac_temp <- extract(temp, Sgrac_range)

Sgrac_temp_all <- unlist(Sgrac_temp)
Sgrac_temp_all
mean(Sgrac_temp_all)#na.rm=TRUE
max(Sgrac_temp_all)
min(Sgrac_temp_all)
#Precip
Sgrac_precip <- extract(precip, Sgrac_range)

Sgrac_precip_all <- unlist(Sgrac_precip)
Sgrac_precip_all
mean(Sgrac_precip_all)#na.rm=TRUE
max(Sgrac_precip_all)
min(Sgrac_precip_all)
