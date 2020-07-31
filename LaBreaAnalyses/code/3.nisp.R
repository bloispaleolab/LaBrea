library(tidyverse)
library(scales)

# arrange taxonomy a bit more
taxonomy <- read.delim("data/raw/TaxonomyMatchingFile.txt", sep="\t", stringsAsFactors = F)
includes.file <- unique(taxonomy$RevisedName)
includes.file <- as.data.frame(includes.file)
colnames(includes.file) <- "RevisedName"
write.csv(includes.file,file="data/processed/includes.file.csv", row.names=F)


# calculate NISP ----

# read in master mammal data 
mammals_orig <- read.delim("data/processed/master_mammal_file.txt", sep="\t", stringsAsFactors = F)

# read in trait taxon names, then match the names to pull out the "include" column
dat <- read.delim("data/processed/Master_trait_taxon_names.txt", sep="\t", stringsAsFactors = F)


# get rid of hancock collection stuff for now
mammals <- mammals_orig[-which(mammals_orig$misc == "y"),]

taxon <- as.data.frame(matrix(ncol=4, nrow=nrow(mammals)))
colnames(taxon) <- c('RevisedName', "Order", "Family", "Genus") 

# for each row (specimen), retrieve the revised taxonomic name, as well as the Order, Family, and Genus info
for (i in 1:nrow(taxon)){
  taxon[i,] <- taxonomy[which(taxonomy$'OriginalName' == mammals$prelim_taxon_name[i]),c('RevisedName', "Order", "Family", "Genus")] 
}

# Code check - make sure all the specimens are matched with new names for analysis
if (any(is.na(taxon$RevisedName))==F){
  print("All rows match a name in our taxon matching file! PROCEED")
}else{
  # examine which specimens in the original mammal data don't match up with our revised taxonomy name:
  mammals[which(is.na(taxon$RevisedName)),]
  # If these mismatches are OK, then PROCEED
  # Otherwise, troubleshoot and fix
}

mammals_trim <- select(mammals, "prelim_taxon_name", "box") # trim down original mammal data to only prelim_taxon_name and box
mammals_trim <- cbind(mammals_trim, taxon) # bind on the updated taxonomy
mammals_trim$Family <- as.factor(mammals_trim$Family) #convert to factor for later filtering
mammals_trim$Order <- as.factor(mammals_trim$Order) #convert to factor for later filtering
mammals_trim$Genus <- as.factor(mammals_trim$Genus) #convert to factor for later filtering
mammals_trim$box <- as.factor(mammals_trim$box) #convert to factor for later filtering

# remove the taxa we don't want to include, this should be verified/changed based on user decisions
mammals_filtered <- mammals_trim %>% 
  filter(!is.na(Family)) %>%
  filter(Order != "Artiodactyla") %>%
  filter(Genus != "Canis") %>%
  filter(Genus != "Taxidea") %>%
  filter(Genus != "Urocyon") %>%
  filter(Genus != "Lepus") %>% ## JESSICA AND NATE check: Why did we include Lepus here? Because it is only to genus?
  filter(!is.na(Genus))

# remove the factor levels that were dropped
mammals_filtered$Order <- droplevels(mammals_filtered$Order)
mammals_filtered$Family <- droplevels(mammals_filtered$Family)
mammals_filtered$Genus <- droplevels(mammals_filtered$Genus)

# count up the specimens to create NISP
mammals_filtered <- mammals_filtered %>%   
  group_by(box, RevisedName, Family, Genus) %>%
  summarise(n=n())

# covert the data from "long" to "wide" format
mammals_nisp <- spread(mammals_filtered, box, n)
mammals_nisp <- mammals_nisp[order(mammals_nisp$Family, mammals_nisp$Genus),]
colnames(mammals_nisp)[4:7] <- paste0("Box_", colnames(mammals_nisp)[4:7])
write.table(mammals_nisp, file="data/processed/mammal_nisp.txt", sep="\t", row.names=F)

### PLOT the NISP data ###
Genus_Colors <- NULL
start <- 0
stop <- 0
for (i in 1:length(levels(mammals_filtered$Family))){
  Genera <- unique(mammals_filtered[which(mammals_filtered$Family==levels(mammals_filtered$Family)[i]),'Genus'])
  nGenera <- nrow(Genera)
  
    if (nGenera <= 3){
      stop <- start+20
    }else{
      stop <- start+50 
    }
  
    show_col(hue_pal(h = c(start, stop))(nGenera))
    Genus_Colors <- rbind(Genus_Colors, cbind(Genera, hue_pal(h = c(start, stop))(nGenera)))
  start <- stop+20
}
colnames(Genus_Colors)[2] <- "color"

# reorder Genus levels to match taxonomy in Genus_Colors
mammals_filtered$Genus <- factor(mammals_filtered$Genus,levels(mammals_filtered$Genus)[match(Genus_Colors[,'Genus'], levels(mammals_filtered$Genus))])

mammals_filtered$box <- factor(mammals_filtered$box, levels = c("7b", "14", "13", "1"))

bp<- ggplot(mammals_filtered, aes(x="", y=n, fill= Genus)) +
  scale_fill_manual(values=rgb(t(col2rgb(Genus_Colors$color, alpha = FALSE)/255))) +
  geom_bar(stat = "identity", position = position_fill()) +   
  facet_grid(.~ box) +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_text(
          size = 18, face = "bold"
        ))

pdf(file="output/mammal_nisp_2.pdf", width=12, height=6, encoding="MacRoman") 
  bp  
dev.off()
