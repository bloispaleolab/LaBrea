## install the latest version from github
install_github('JohnsonHsieh/iNEXT')
library(iNEXT)
library(ggplot2)
library(vegan)
library(BiodiversityR)
library(devtools)

#Example
data(BCI)
H <- diversity(BCI)
simp <- diversity(BCI, "simpson")
invsimp <- diversity(BCI, "inv")
r.2 <- rarefy(BCI, 2)
alpha <- fisher.alpha(BCI)
pairs(cbind(H, simp, invsimp, r.2, alpha), pch="+", col="blue")
## Species richness (S) and Pielou's evenness (J):
S <- specnumber(BCI) ## rowSums(BCI > 0) does the same...
J <- H/log(S)

#Ages: Box 1 = ~37,970 Cal BP; Box 13 = ~39,880 Cal BP; Box 7b = 43,170 Cal BP;
#Box 14 = >46,000 Cal BP

#Box 14 = warmest mean, but rapidly cooling; Box 7b = 2nd coldest mean temp, relative stability;
#Box 13 = coldest mean, rapid cooling event; Box 1 = second warmest mean, rapid warming event)

##loadData
Div.data <- read.csv("data/processed/P23_diversity_data_formatted.csv", header = T, row.names = 1)

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
data(dune.env)
data(dune)
RankAbun.1 <- rankabundance(dune)
RankAbun.1
rankabunplot(RankAbun.1, scale='abundance', addit=FALSE, specnames=c(1,2,3))
rankabunplot(RankAbun.1, scale='logabun', addit=FALSE, specnames=c(1:30), 
             srt=45, ylim=c(1,100))
rankabuncomp(dune, y=dune.env, factor='Management', 
             scale='proportion', legend=FALSE)
#My data
Dep1 <-Div.data[1, ]
Dep13 <-Div.data[2, ]
Dep7b <-Div.data[3, ]
Dep14 <-Div.data[4, ]

RankAbun.1 <- rankabundance(Dep14)#enter deposit of interest
RankAbun.1
#rank abundance
rankabunplot(RankAbun.1, scale='abundance', addit=FALSE, specnames=c(1,2,3,4))
#log abundance
rankabunplot(RankAbun.1, scale='logabun', addit=FALSE, specnames=c(1:15), 
             srt=45, xlim=c(1, 14), ylim=c(1,1000))

#### i_NEXT ###

## install the latest version from github

iNEXT.data <- read.csv("data/processed/P23_diversity_iNEXT.csv", row.names = 1)

out <- iNEXT(iNEXT.data, q=c(0,1,2), datatype="abundance", conf=0.95, endpoint=1000)
# q=0 = taxon richness, q=1 = shannon diversity, q=2 = Simpson diversity 
# Hill numbers can be viewed all togther or individually

ggiNEXT(out, type=1, se=TRUE, facet.var="order") 
# type can be changed from 1-3
# type=1 = Sample-size-based R/E curve
# type=2 = Sample completeness curve
# type=3 = Coverage-based R/E curve

#facet.var = "none", "order", "site" or "both"
