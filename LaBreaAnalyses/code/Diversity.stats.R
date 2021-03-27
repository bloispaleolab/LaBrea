## install the latest version from github
#install_github('JohnsonHsieh/iNEXT')
library(devtools)
library(iNEXT)
library(ggplot2)
library(vegan)
library(BiodiversityR)
library("dplyr")

#Example
# data(BCI)
# H <- diversity(BCI)
# simp <- diversity(BCI, "simpson")
# invsimp <- diversity(BCI, "inv")
# r.2 <- rarefy(BCI, 2)
# alpha <- fisher.alpha(BCI)
# pairs(cbind(H, simp, invsimp, r.2, alpha), pch="+", col="blue")
# ## Species richness (S) and Pielou's evenness (J):
# S <- specnumber(BCI) ## rowSums(BCI > 0) does the same...
# J <- H/log(S)

#Ages: Box 1 = ~37,970 Cal BP; Box 13 = ~39,880 Cal BP; Box 7b = 43,170 Cal BP;
#Box 14 = >46,000 Cal BP

#Box 14 = warmest mean, but rapidly cooling; Box 7b = 2nd coldest mean temp, relative stability;
#Box 13 = coldest mean, rapid cooling event; Box 1 = second warmest mean, rapid warming event)

# read.nisp and make sure included taxa are correct
nisp <- read.delim(file="data/processed/mammal_nisp.txt", sep="\t", header=T)

# add the "include" to match new "include" file
include <- read.csv(file="data/processed/includes.file.csv", header=T)

nisp <- cbind(nisp, include[match(nisp$RevisedName, include$RevisedName),'Include'])
colnames(nisp)[ncol(nisp)] <- "include"

nisp <- nisp[which(nisp$include=="y"),] # note all data in nisp are included
nisp[is.na(nisp)] <- 0
# sum data by Genus
genus_nisp <- nisp %>%                                            # Specify data frame
  group_by(Genus) %>%                               # Specify group indicator
  summarise_at(vars(Box_1, Box_13, Box_14, Box_7b), # Specify column
               list(sum))                    # Specify function

# arrange in iNEXT format
div.data.genus <- t(genus_nisp[,2:5])
colnames(div.data.genus) <- genus_nisp$Genus
div.data.genus[which(is.na(div.data.genus))] <- 0

Div.data <- div.data.genus

##loadData - Nate's original code
# Div.data <- read.csv("data/processed/P23_diversity_data_formatted.csv", header = T, row.names = 1)

##Richness
apply(Div.data>0,1,sum)

##Richness corrected for sampling effort

#Menhinick's index (D = n/sqrtN)
n<-apply(Div.data>0,1,sum)
N <- apply(Div.data,1,sum)
n/sqrt(N)

#Margalef's index (D = n-1/lnN)
n<-apply(Div.data>0,1,sum)
N <- apply(Div.data,1,sum)
(n-1)/log(N)

#Abundances
apply(Div.data,1,sum)

#Rarified abundances
rarefy(Div.data, sample=100, MARGIN=1)

#Shannon_index
diversity(Div.data, index="shannon")

#Hills's numbers (evenness)
S <- apply(Div.data>0,1,sum)
exp(diversity(Div.data, index="shannon"))/S

#rank abundance curves
mod <- rad.lognormal(Div.data)
mod
plot(mod)
mod <- radfit(Div.data)
## Standard plot overlaid for all models
## Pre-emption model is a line
plot(mod)
## log for both axes: Zipf model is a line
plot(mod, log = "xy")


#example
# data(dune.env)
# data(dune)
# RankAbun.1 <- rankabundance(dune)
# RankAbun.1
# rankabunplot(RankAbun.1, scale='abundance', addit=FALSE, specnames=c(1,2,3))
# rankabunplot(RankAbun.1, scale='logabun', addit=FALSE, specnames=c(1:30), 
#              srt=45, ylim=c(1,100))
# rankabuncomp(dune, y=dune.env, factor='Management', 
#              scale='proportion', legend=FALSE)

#My data
Dep1 <- t(as.matrix(Div.data[1, ]))
Dep13 <- t(as.matrix(Div.data[2, ]))
Dep14 <- t(as.matrix(Div.data[3, ]))
Dep7b <- t(as.matrix(Div.data[4, ]))


RankAbun.1 <- rankabundance(Dep13)#enter deposit of interest
RankAbun.1
#rank abundance
pdf("Output/diversity/Dep13_RA.pdf", width=5, height=5)
rankabunplot(RankAbun.1, scale='abundance', addit=FALSE, specnames=c(1,2,3,4))
dev.off()
#log abundance
rankabunplot(RankAbun.1, scale='logabun', addit=FALSE, specnames=c(1:15), 
             srt=45, xlim=c(1, 14), ylim=c(1,1000))

#### i_NEXT ###

## install the latest version from github

iNEXT.data <- read.csv("data/processed/P23_div_iNEXT_spp.csv", row.names = 1) # Note: this was Nate's original code, replaced below by nisp output

iNEXT.data <- nisp[,c('Box_1', 'Box_13', 'Box_14', 'Box_7b')]
# iNEXT.data <- apply(iNEXT.data, 2, function(x) replace_na(x, 0))
rownames(iNEXT.data) <- nisp$'RevisedName'

out <- iNEXT(iNEXT.data, q=c(1), datatype="abundance", conf=0.95, endpoint=1000)
# q=0 = taxon richness, q=1 = shannon diversity, q=2 = Simpson diversity 
# Hill numbers can be viewed all togther or individually
pdf("Output/diversity/inext_spp.pdf", width=7, height=4)
div_plot<-ggiNEXT(out, type=1, se=TRUE) 
div_plot + labs(y = "Species diversity")
dev.off()
# type can be changed from 1-3
# type=1 = Sample-size-based R/E curve
# type=2 = Sample completeness curve
# type=3 = Coverage-based R/E curve

#facet.var = "none", "order", "site" or "both" 
# removed this from the ggiNEXT code: facet.var="order" 
# don't need this

# just plot the data without Box 14

iNEXT.data.sub <- read.csv("data/processed/P23_div_iNEXT_gen.csv", row.names = 1) # Note: this was Nate's original code, replaced below by nisp output

iNEXT.data.sub <- iNEXT.data.sub[,c('Deposit.1', 'Deposit.13', 'Deposit.7B')]
#iNEXT.data.sub <- nisp[,c('Box_1', 'Box_13', 'Box_7b')]
# iNEXT.data <- apply(iNEXT.data, 2, function(x) replace_na(x, 0))
#rownames(iNEXT.data.sub) <- nisp$'RevisedName'

out.sub <- iNEXT(iNEXT.data.sub, q=c(1), datatype="abundance", conf=0.95, endpoint=1000)

pdf("Output/diversity/inext_spp_no14.pdf", width=7, height=4)
div_plot<-ggiNEXT(out.sub, type=1, se=TRUE, color.var = "site") 
div_plot +   
  scale_colour_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=14),
        axis.title=element_text(size=18,face="bold")) +
  labs(y = "Species diversity")
dev.off()
