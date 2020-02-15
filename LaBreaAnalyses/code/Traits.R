library(tidyverse)
library(scales)

# calculate NISP ----

# read in master mammal data 
mammals <- read.delim("data/processed/master_mammal_NISP.txt", sep="\t", stringsAsFactors = F)
taxonomy_traits <- read.delim("data/processed/TaxonomyTraitsFile.txt", sep="\t", stringsAsFactors = F)


taxon <- as.data.frame(matrix(ncol=8, nrow=nrow(mammals)))
  colnames(taxon) <- c('RevisedName', "Order", "Family", "Genus", "Precipitation.Affinity", 
  "Temperature.Affinity", "Diet.Group", "Diet.Breadth") 

for (i in 1:nrow(taxon)){
  taxon[i,] <- taxonomy_traits[which(taxonomy_traits$'OriginalName' == mammals$prelim_taxon_name[i]),
  c('RevisedName', "Order", "Family", "Genus", "Precipitation.Affinity", "Temperature.Affinity", 
  "Diet.Group", "Diet.Breadth")] 
} 



taxon[which(taxon$Genus=="Lepus"),'Genus'] <- "Sylvilagus"


mammals_trim <- select(mammals, "prelim_taxon_name", "box")
mammals_trim <- cbind(mammals_trim, taxon)
mammals_trim$Family <- as.factor(mammals_trim$Family)
mammals_trim$Order <- as.factor(mammals_trim$Order)
mammals_trim$Genus <- as.factor(mammals_trim$Genus)
mammals_trim$box <- as.factor(mammals_trim$box)


mammals_filtered <- mammals_trim %>% 
  ## Filter out specimens with no trait data
  filter(!is.na(Precipitation.Affinity)) %>%
  filter(!is.na(Family)) %>%
  filter(Order != "Artiodactyla") %>%
  filter(Genus != "Canis") %>%
  filter(Genus != "Taxidea") %>%
  filter(!is.na(Genus))
mammals_filtered$Order <- droplevels(mammals_filtered$Order)
mammals_filtered$Family <- droplevels(mammals_filtered$Family)
mammals_filtered$Genus <- droplevels(mammals_filtered$Genus)


mammals_filtered_test <- mammals_filtered %>%   
  group_by(box, Temperature.Affinity) %>%
  summarise(n=n())

#Dont run
#Genus_Colors <- NULL
#start <- 0
#stop <- 0
#for (i in 1:length(levels(mammals_filtered$Family))){
  #Genera <- unique(mammals_filtered[which(mammals_filtered$Family==levels(mammals_filtered$Family)[i]),'Genus'])
  #nGenera <- nrow(Genera)
  
  #if (nGenera <= 3){
    #stop <- start+20
  #}else{
    #stop <- start+50 
  #}
  
  #show_col(hue_pal(h = c(start, stop))(nGenera))
  #Genus_Colors <- rbind(Genus_Colors, cbind(Genera, hue_pal(h = c(start, stop))(nGenera)))
  #start <- stop+20
#}
#colnames(Genus_Colors)[2] <- "color"

# reorder Genus levels to match taxonomy in Genus_Colors
#mammals_filtered$Genus <- factor(mammals_filtered$Genus,levels(mammals_filtered$Genus)[match(Genus_Colors[,'Genus'], levels(mammals_filtered$Genus))])


bp<- ggplot(mammals_filtered_test, aes(x="", y=n, fill= Temperature.Affinity)) +
  #scale_fill_manual(values=rgb(t(col2rgb(Genus_Colors$color, alpha = FALSE)/255))) +
  geom_bar(stat = "identity", position = position_fill()) + 
  coord_polar("y", start=0) + 
  facet_grid(.~ box) +
  theme_light()

pdf(file="output/mammal_nisp_temp.pdf", width=15, height=7.5, encoding="MacRoman") 
bp  
dev.off()


#stats

(library(vegan))
(library(reshape2))
options(scipen=999)

Trans <- dcast(mammals_filtered_test, box ~ Temperature.Affinity, value.var = "n")
Trans <- Trans[, -1]
Trans <- as.table(as.matrix(Trans))

#For missing temp
Trans[3,1]=0


dimnames(Trans) <- list(box = c("1", "13","14", "7b"),
                        Temperature.Affinity = c("Cool","Neutral", "Warm"))

Xsq<-chisq.test(Trans)
Xsq

#Doesnt work
adonis(formula=Trans[, -1]~box, data=Trans[,-1], method=dist)








