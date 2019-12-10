library(tidyverse)
library(scales)

# calculate NISP ----

# read in master mammal data 
mammals_orig <- read.delim("data/processed/master_mammal_file.txt", sep="\t", stringsAsFactors = F)

# get rid of hancock collection stuff for now
mammals <- mammals_orig[-which(mammals_orig$box == "HC"),]

taxonomy <- read.delim("data/raw/TaxonomyMatchingFile.txt", sep="\t", stringsAsFactors = F)

taxon <- as.data.frame(matrix(ncol=4, nrow=nrow(mammals)))
colnames(taxon) <- c('RevisedName', "Order", "Family", "Genus") 

for (i in 1:nrow(taxon)){
  taxon[i,] <- taxonomy[which(taxonomy$'OriginalName' == mammals$prelim_taxon_name[i]),c('RevisedName', "Order", "Family", "Genus")] 
}

mammals_trim <- select(mammals, "prelim_taxon_name", "box")
mammals_trim <- cbind(mammals_trim, taxon)
mammals_trim$Family <- as.factor(mammals_trim$Family)
mammals_trim$Order <- as.factor(mammals_trim$Order)
mammals_trim$Genus <- as.factor(mammals_trim$Genus)
mammals_trim$box <- as.factor(mammals_trim$box)

mammals_filtered <- mammals_trim %>% 
  filter(!is.na(Family)) %>%
  filter(Order != "Carnivora") %>%
  filter(Order != "Artiodactyla") %>%
  filter(!is.na(Genus))
mammals_filtered$Order <- droplevels(mammals_filtered$Order)
mammals_filtered$Family <- droplevels(mammals_filtered$Family)
mammals_filtered$Genus <- droplevels(mammals_filtered$Genus)

mammals_filtered <- mammals_filtered %>%   
  group_by(box, RevisedName, Family, Genus) %>%
  summarise(n=n())

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

mammals_filtered$box <- factor(mammals_filtered$box, levels = c("14", "7b", "13", "1"))

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

pdf(file="output/mammal_nisp.pdf", width=12, height=6, encoding="MacRoman") 
  bp  
dev.off()


# rarefaction curves
library(iNEXT)
mamm <- list(mammals_filtered$'box', mammals_filtered$'n')
str(mamm)
names(mamm) <- c("box", "abundance")
out <- iNEXT(mamm, q=0, datatype="abundance")
spider.girdled <- iNEXT(spider$Girdled, q=0, datatype="abundance")


# Sample-size-based R/E curves, separating by "site""
ggiNEXT(out, type=1, facet.var=box)
