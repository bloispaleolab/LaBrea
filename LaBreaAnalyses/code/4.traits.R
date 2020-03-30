library(tidyverse)
library(scales)

# calculate NISP ----

# read in master mammal data 
mammals_orig <- read.delim("data/processed/master_mammal_file.txt", sep="\t", stringsAsFactors = F)
mammals <- mammals_orig[-which(mammals_orig$misc == "y"),]

# Estimate Sylvilagus spp contributions

#which(mammals$prelim_taxon_name == "Sylvilagus audubonii" & mammals$box == "7b")
#aud<-16
#which(mammals$prelim_taxon_name == "Sylvilagus bachmani" & mammals$box == "7b")
#bac<-1
#aud_ratio=aud/(aud+bac)
#bac_ratio=bac/(aud+bac)
#Syl_sp<-which(mammals$prelim_taxon_name == "Sylvilagus sp" & mammals$box == "7b")
#Syl_aud<-76*aud_ratio
#Syl_bac<-76*bac_ratio

# this will provide list of row indexes with all Sylv. sp frm a box
# convert percentage to a whole number based on total number of row indices
# replace !st x% of the rows with one spp
# replace the remainder of the rows with the other spp. (not whole row, just prelim_taxon_name)

taxonomy_traits <- read.delim("data/processed/TaxonomyTraitsFile_Reg_est.txt", sep="\t", stringsAsFactors = F)


taxon <- as.data.frame(matrix(ncol=9, nrow=nrow(mammals)))
colnames(taxon) <- c('RevisedName', "Order", "Family", "Genus", "Species", 
  "BodyMass_Mean", "Precip_Mean", "Temp_Mean", "NDVI_Mean") 



for (i in 1:nrow(taxon)){
  taxon[i,] <- taxonomy_traits[which(taxonomy_traits$'OriginalName' == mammals$prelim_taxon_name[i]),
  c('RevisedName', "Order", "Family", "Genus", "Species", "BodyMass_Mean", "Precip_Mean", 
  "Temp_Mean", "NDVI_Mean")] 
}


mammals_trim <- select(mammals, "prelim_taxon_name", "box")
mammals_trim <- cbind(mammals_trim, taxon)
mammals_trim$Family <- as.factor(mammals_trim$Family)
mammals_trim$Order <- as.factor(mammals_trim$Order)
mammals_trim$Genus <- as.factor(mammals_trim$Genus)
mammals_trim$box <- as.factor(mammals_trim$box)

mammals_filtered <- mammals_trim %>% 
  filter(!is.na(Family)) %>%
  filter(Order != "Artiodactyla") %>%
  filter(Genus != "Canis") %>%
  filter(Genus != "Taxidea") %>%
  filter(Genus != "Urocyon") %>%
  filter(Genus != "Lepus") %>%
  filter(!is.na(Genus))
mammals_filtered$Order <- droplevels(mammals_filtered$Order)
mammals_filtered$Family <- droplevels(mammals_filtered$Family)
mammals_filtered$Genus <- droplevels(mammals_filtered$Genus)

mammals_filtered <- mammals_filtered %>%   
  group_by(box, Precip_Mean) %>%
  summarise(n=n()) #"BodyMass_Mean", "Precip_Mean", "Temp_Mean", "NDVI_Mean"

### Apply estimated Sylvilagus species contributions
mammals_filtered$n[3] <- mammals_filtered$n[3]+115 #box 1 audubonii
mammals_filtered$n[9] <- mammals_filtered$n[9]+25 #box 1 bachmani
mammals_filtered$n[15] <- mammals_filtered$n[15]+204 #box 13 audubonii
mammals_filtered$n[24] <- mammals_filtered$n[24]+99 #box 14 audubonii
mammals_filtered$n[31] <- mammals_filtered$n[31]+50 #box 14 bachmani
mammals_filtered$n[36] <- mammals_filtered$n[36]+72 #box 7b audubonii
mammals_filtered$n[43] <- mammals_filtered$n[43]+4  #box 7b bachmani

###Stats###

Box1<-mammals_filtered[1:11,]
Box13<-mammals_filtered[12:21,]
Box14<-mammals_filtered[22:33,]
Box7b<-mammals_filtered[34:45,]

#Box1
n_Box1<-sum(Box1$n[1:10])
Box1_weighted<-Box1$Precip_Mean*Box1$n
Box1_sum<-sum(Box1_weighted, na.rm=TRUE)
Box1_avg<-Box1_sum/n_Box1
Box1_avg
#Box13
n_Box13<-sum(Box13$n[1:9])
Box13_weighted<-Box13$Precip_Mean*Box13$n
Box13_sum<-sum(Box13_weighted, na.rm=TRUE)
Box13_avg<-Box13_sum/n_Box13
Box13_avg
#Box14
n_Box14<-sum(Box14$n[1:11])
Box14_weighted<-Box14$Precip_Mean*Box14$n
Box14_sum<-sum(Box14_weighted, na.rm=TRUE)
Box14_avg<-Box14_sum/n_Box14
Box14_avg
#Box7b
n_Box7b<-sum(Box7b$n[1:11])
Box7b_weighted<-Box7b$Precip_Mean*Box7b$n
Box7b_sum<-sum(Box7b_weighted, na.rm=TRUE)
Box7b_avg<-Box7b_sum/n_Box7b
Box7b_avg
