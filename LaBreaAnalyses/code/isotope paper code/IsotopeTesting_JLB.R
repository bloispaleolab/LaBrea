library(dplyr)
library(stringr)
library(ggplot2)

# random function so I can extract p-values from lm
Regressionp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
  
# Section 1: Import and Clean data ----
## import isotope data ----
isoDat2 <- read.csv('data/processed/SIBER/SIBER_raw.csv', header=T, strip.white=T)

# create sample names for matching isotope file to ages file
samples <- paste("UCIAMS", isoDat2$UCIAMS_Number) # this is the final set of samples with isotope data

## read in calibrated ages ----
allAges<- read.csv('output/OxCal/final oxcal models/AllAges_ages_probs.csv', header=T) # this file stores the raw probabilities for each age across the distribution
all_calibrated_ages <- read.csv('output/OxCal/final oxcal models/AllAges_forinput.csv', header=T) # this file stores the calibrated age statistics for each age
all_calibrated_ages$trimmedName <- unlist(lapply(strsplit(all_calibrated_ages$Name, " R_"), '[[', 1))
sample_median_ages <- as.data.frame(cbind(samples, all_calibrated_ages[match(samples, all_calibrated_ages$trimmedName), 'Unmodelled._BP_median']))
colnames(sample_median_ages)[2] <- 'median_age'  
sample_median_ages$median_age <- as.numeric(sample_median_ages$median_age)
  
## data cleaning on allAges file----
allAges<- allAges[allAges$probability != 0, ] # remove age estimates with 0 probability

#remove repetitive information and duplicate samples 
allAges<- allAges[!(allAges$name=="R Combine 1 LACMP23-33228" | 
                     allAges$name=="R Combine 2 LACMHC-142773" |
                     allAges$name=="R Combine 3 LACMHC-142779" |
                     allAges$name=="R Combine 4 LACMP23-35541" |
                     allAges$name=="R Combine 5 LACMP23-40642" |
                     allAges$name=="UCIAMS 198206" | allAges$name=="UCIAMS 223585" |
                     allAges$name=="UCIAMS 217076" | allAges$name=="UCIAMS 217077" |  
                     allAges$name=="UCIAMS 223587"),]

#remove specimens that are not rabbits and squirrels
allAges<- allAges[!(allAges$name=="UCIAMS 198199" | 
                     allAges$name=="UCIAMS 198200" | allAges$name=="UCIAMS 198208" |
                     allAges$name=="UCIAMS 216782" | allAges$name=="UCIAMS 223520" |
                     allAges$name=="UCIAMS 198297"),]

#remove specimens without stable isotope values & unmeasureable dates
allAges<- allAges[!(allAges$name=="UCIAMS 223515" | 
                     allAges$name=="UCIAMS 223519" | allAges$name=="UCIAMS 223522"),]

# Section 2: Match d180 to all age estimates ----
# (both median_age as well as the full distribution of ages) for every sample 

if (all(allAges$value<0)==T){ # this should be true, then take the absolute values to match to climate
  xout <- abs(allAges$value)
}

# first match hendy data to age data
  hendyDat<- read.delim("data/raw/climate/hendy2002data.txt")
  
  # adding d18O to allAges to match to the full distribution of age estimates for each sample
  Hendy_extracted <- approx(x=hendyDat$HendyAge, y=hendyDat$pach.d18O, method="linear", xout=xout)
  d18O_hendy <-Hendy_extracted$y
  allAges_d18O<-cbind(allAges, d18O_hendy) # this matches d180 to all the age estimates, for all probabilities
  
  # matching d18O to the median ages
  d18O_hendy_medianage <- approx(x=hendyDat$HendyAge, y=hendyDat$pach.d18O, 
                           method="linear", xout=sample_median_ages$median_age)$y

  # next match ngrip data to age data
  ngripDat<- read.delim("data/raw/climate/ngrip.txt")
  ngripDat$Age <- ngripDat$Age-50 #ngrip dates are "before 2000", so need to subtract 50 years
  
  # adding d18O to allAges to match to the full distribution of age estimates for each sample
  ngrip_extracted <- approx(x=ngripDat$Age, y=ngripDat$d18O, method="linear", xout=xout)
  d18O_ngrip <-ngrip_extracted$y
  allAges_d18O <- cbind(allAges_d18O, d18O_ngrip) # this matches d180 to all the age estimates, for all probabilities
  
  # matching d18O to the median ages
  d18O_ngrip_medianage <- approx(x=ngripDat$Age, y=ngripDat$d18O, 
                           method="linear", xout=sample_median_ages$median_age)$y

# clean up matched d18O data
  # d18O_hendy # not needed - bound to allAges_d18O
  # d18O_ngrip # not needed - bound to allAges_d18O
  rm(d18O_hendy)
  rm(d18O_ngrip)
  
  # d18O_hendy_medianage - matched to median age within sample_median_ages
  # d18O_ngrip_medianage - matched to median age within sample_median_ages
  sample_median_ages <- cbind(sample_median_ages, d18O_hendy_medianage, d18O_ngrip_medianage)
  colnames(sample_median_ages)[3:4] <- c("d18O_hendy", "d18O_ngrip")
  rm(d18O_hendy_medianage)
  rm(d18O_ngrip_medianage)
  
  ### at end of Section 2, have 2 main files for later use: ###

  # allAges_d18O (stores both ngrip and hendy)
  # sample_median_ages (stores both ngrip and hendy)
  
# Section 3: additional data cleaning ----
# remove any unmatched specimen from the files
matches <- match(samples, allAges_d18O$name)

# added ifelse statement to catch instances where there are or are not matches
if (any(is.na(matches))){
  isoDat2_final <- isoDat2[-which(is.na(matches)),]
  samples_final <- samples[-which(is.na(matches))]
}else{
  isoDat2_final <- isoDat2
  samples_final <- samples
}
# end result dataframes:
# isoDat2_final - final C/N isotope data
# samples_final - final set of samples
# allAges_d18O - ages, probabilities, and d18O for the full distribution of calibrated age estimates associated with each sample
# sample_median_ages - median ages and d180 estimates at those median ages for each sample

# Section 4: Age uncertainty: calculate weighted age and weighted mean d18O for each specimen ----
# For each specimen, calculate the weighted d180 and weighted mean age, based on the probability of each calibrated age from the age models

specimen_wage <- vector(mode="numeric", length=length(samples_final))
specimen_wd18O_hendy <- vector(mode="numeric", length=length(samples_final))
specimen_wd18O_ngrip <- vector(mode="numeric", length=length(samples_final))

for (i in 1:length(samples_final)){
  # pull out specimen ID
  spec <- samples_final[i] 
  # pull out all calibrated ages and probabilities for that single specimen
  specd18O <- allAges_d18O[which(allAges_d18O$name == spec),]
  # calculate weighted mean d18O
  specimen_wd18O_hendy[i] <- weighted.mean(specd18O$d18O_hendy, specd18O$probability)
  specimen_wd18O_ngrip[i] <- weighted.mean(specd18O$d18O_ngrip, specd18O$probability)
  # calculate weighted mean age
  specimen_wage[i] <- weighted.mean(specd18O$value, specd18O$probability)
}

# match all age and d18O estimates with the original C&N isotope data
matchedDF_weighted <- cbind(isoDat2_final[,c('UCIAMS_Number', 'Taxon', 'del15N_permil', 'del13C_permil', 'X14C_age_BP')], specimen_wage, specimen_wd18O_hendy, specimen_wd18O_ngrip)

# make sure still in same order
match(gsub("UCIAMS ", "", sample_median_ages$samples), matchedDF_weighted$UCIAMS_Number)

matchedDF_all <- cbind(matchedDF_weighted, sample_median_ages[,2:4]) 
# matchedDF_all thus is the full final dataset:
# UCIAMS_Number
# Taxon
# del15N_permil      
# del13C_permil
# X14C_age_BP
# specimen_wage # sensitivity age
# specimen_wd18O_hendy # sensitivity d18O estimate - hendy
# specimen_wd18O_ngrip # sensitivity d18O estimate - ngrip
# median_age # primary age
# d18O_hendy #primary d18O estimate at median age - hendy
# d18O_ngrip #primary d18O estimate at median age - hendy

######### Section 5: Final Models and Figures ##############
# using the final dataframe: matchedDF_all
# we should use the median d18O values for final analysis, and present the weighted ones in supplemental, along with sensitivity analysis

# Methods potential text:
# For carbon and nitrogen, we fit a linear model that originally included oxygen, taxon, and the interaction between the two variables as independent variables. We then performed stepwise regression to determine a final model. 
# Results potential text
# Stepwise regression indicated that there was no significant interaction between 18O and taxon for either carbon or nitrogen stable isotope values. For carbon, variation in 13C was significantly associated with both 18O and taxon (stats from summary(carbon.lm.final)). For nitrogen, neither 18O nor taxon explained significant variation in 15N, though taxon as a single variable was marginally significant (stats from summary(nitrogen.lm.taxon)).

final.model.stats <- as.data.frame(matrix(data=NA, nrow=4, ncol=5))
colnames(final.model.stats) <- c("element", "climate", "final_model", "p", "adjR2")

## start with Hendy climate data ----
### models - Carbon ----

# Start simple and build complexity

# climate-only model
carbon.lm.clim<-lm(del13C_permil~d18O_hendy, matchedDF_all)
summary(carbon.lm.clim)
plot(del13C_permil ~ d18O_hendy, data=matchedDF_all, pch=16)
abline(carbon.lm.clim)

# taxon-only model
carbon.lm.taxon<-lm(del13C_permil~Taxon, data=matchedDF_all)
summary(carbon.lm.taxon)

# model with no interaction term 
# NOTE: this is what the final model is in the end, so should be the same as carbon.lm.final
carbon.lm.all.additive<-lm(del13C_permil~d18O_hendy + Taxon, data=matchedDF_all)
summary(carbon.lm.all.additive)

# model with interaction term included
carbon.lm.all.interaction<-lm(del13C_permil~d18O_hendy * Taxon, 
                              data=matchedDF_all)
summary(carbon.lm.all.interaction)

##NOTE - no significant interactions, total amount of variation explained is lower than with additive model

#stepwise regression 
carbon.lm.final <- step(lm(del13C_permil~d18O_hendy * Taxon, data=matchedDF_all), direction="both")
summary(carbon.lm.final)
# --> this shows that the d180 + Taxon model (no interaction) is the best final model. 

# t-test of residuals from climate only model
# residuals of Otospermophilus and Sylvilagus are signif different
c.t.clim <- t.test(carbon.lm.clim$residuals ~ matchedDF_all$Taxon[which(!is.na(matchedDF_all$d18O_hendy))])
c.t.clim 

### models - Nitrogen ----

# Start simple and build complexity

# climate-only model
nitrogen.lm.clim<-lm(del15N_permil ~ d18O_hendy, data=matchedDF_all)
summary(nitrogen.lm.clim)

# taxon-only model
nitrogen.lm.taxon<-lm(del15N_permil~Taxon, data=matchedDF_all)
summary(nitrogen.lm.taxon)

# model with no interaction term 
nitrogen.lm.all.additive<-lm(del15N_permil ~ d18O_hendy + Taxon, data=matchedDF_all)
summary(nitrogen.lm.all.additive)

# model with interaction term included
nitrogen.lm.all.interaction<-lm(del15N_permil ~ d18O_hendy * Taxon, data=matchedDF_all)
summary(nitrogen.lm.all.interaction)

#stepwise regression --> this shows that there is not a great final model
nitrogen.lm.final <- step(lm(del15N_permil~d18O_hendy*Taxon, data=matchedDF_all), direction="both")
summary(nitrogen.lm.final)
# --> the best final model is the additive model, which is the same as carbon.

# t-test of residuals from climate-only model
# residuals of Otospermophilus and Sylvilagus are NOT signif different
n.t.clim <- t.test(nitrogen.lm.clim$residuals ~ matchedDF_all$Taxon[which(!is.na(matchedDF_all$d18O_hendy))])
n.t.clim

### Do we need to add age as a co-variate in the models? ----
# Does this answer the question: "does isotope vary through time?"
carbon.lm.all.age<-step(lm(del13C_permil~d18O_hendy * Taxon * median_age, data=matchedDF_all), direction="both")
summary(carbon.lm.all.age)
# explain a bit more variation with age included: 0.4877 vs 0.428

nitrogen.lm.all.age<-step(lm(del15N_permil ~ d18O_hendy * Taxon * median_age, data=matchedDF_all), direction="both")
summary(nitrogen.lm.all.age)
# explain a bit more variation with age included: 0.1702 vs 0.08476

final.model.stats[1,] <- c("carbon", "hendy", "carbon.lm.all.additive", 
                           Regressionp(carbon.lm.all.additive), summary(carbon.lm.all.additive)$adj)
final.model.stats[2,] <- c("nitrogen", "hendy", "nitrogen.lm.all.additive", 
                           Regressionp(nitrogen.lm.all.additive), summary(nitrogen.lm.all.additive)$adj)

# rename final models before changing climate data
carbon.lm.final.hendy <- carbon.lm.final
carbon.lm.clim.hendy <- carbon.lm.clim
c.t.clim.hendy <- c.t.clim
nitrogen.lm.final.hendy <- nitrogen.lm.final
nitrogen.lm.clim.hendy <- nitrogen.lm.clim
n.t.clim.hendy <- n.t.clim

# remove variables that are being generated anew in ngrip section
rm("c.t.clim", "carbon.lm.all.additive", "carbon.lm.all.age", "carbon.lm.all.interaction", "carbon.lm.clim", "carbon.lm.final", "carbon.lm.taxon",
   "n.t.clim", "nitrogen.lm.all.additive", "nitrogen.lm.all.age", "nitrogen.lm.all.interaction", "nitrogen.lm.clim", "nitrogen.lm.final", "nitrogen.lm.taxon")

## next with ngrip climate data ----
### models - Carbon ----

# Start simple and build complexity

# climate-only model
carbon.lm.clim<-lm(del13C_permil~d18O_ngrip, matchedDF_all)
summary(carbon.lm.clim)
plot(del13C_permil ~ d18O_ngrip, data=matchedDF_all, pch=16)
abline(carbon.lm.clim)

# taxon-only model
carbon.lm.taxon<-lm(del13C_permil~Taxon, data=matchedDF_all)
summary(carbon.lm.taxon)

# model with no interaction term 
# NOTE: this is what the final model is in the end, so should be the same as carbon.lm.final
carbon.lm.all.additive<-lm(del13C_permil~d18O_ngrip + Taxon, data=matchedDF_all)
summary(carbon.lm.all.additive)

# model with interaction term included
carbon.lm.all.interaction<-lm(del13C_permil~d18O_ngrip * Taxon, 
                              data=matchedDF_all)
summary(carbon.lm.all.interaction)

##NOTE - no significant interactions, total amount of variation explained is about same as additive model

#stepwise regression 
carbon.lm.final <- step(lm(del13C_permil~d18O_ngrip * Taxon, data=matchedDF_all), direction="both")
summary(carbon.lm.final)
# --> this shows that the d180 + Taxon model (no interaction, additive model) is the best final model. 

# t-test of residuals from climate only model
# residuals of Otospermophilus and Sylvilagus are signif different
c.t.clim <- t.test(carbon.lm.clim$residuals ~ matchedDF_all$Taxon[which(!is.na(matchedDF_all$d18O_ngrip))])
c.t.clim 

### models - Nitrogen ----

# Start simple and build complexity

# climate-only model
nitrogen.lm.clim<-lm(del15N_permil ~ d18O_ngrip, data=matchedDF_all)
summary(nitrogen.lm.clim)

# taxon-only model
nitrogen.lm.taxon<-lm(del15N_permil~Taxon, data=matchedDF_all)
summary(nitrogen.lm.taxon)

# model with no interaction term 
nitrogen.lm.all.additive<-lm(del15N_permil ~ d18O_ngrip + Taxon, data=matchedDF_all)
summary(nitrogen.lm.all.additive)

# model with interaction term included
nitrogen.lm.all.interaction<-lm(del15N_permil ~ d18O_ngrip * Taxon, data=matchedDF_all)
summary(nitrogen.lm.all.interaction)

#stepwise regression --> this shows that there is not a great final model
nitrogen.lm.final <- step(lm(del15N_permil~d18O_ngrip*Taxon, data=matchedDF_all), direction="both")
summary(nitrogen.lm.final)
# --> the best final model is the additive model, which is the same as carbon.

# t-test of residuals from climate-only model
# residuals of Otospermophilus and Sylvilagus are NOT signif different
n.t.clim <- t.test(nitrogen.lm.clim$residuals ~ matchedDF_all$Taxon[which(!is.na(matchedDF_all$d18O_ngrip))])
n.t.clim

### Do we need to add age as a co-variate in the models? ----
# Does this answer the question: "does isotope vary through time?"
carbon.lm.all.age<-step(lm(del13C_permil~d18O_ngrip * Taxon * median_age, data=matchedDF_all), direction="both")
summary(carbon.lm.all.age)
# explain a bit more variation with age included: 0.4902 vs 0.4063

nitrogen.lm.all.age<-step(lm(del15N_permil ~ d18O_ngrip * Taxon * median_age, data=matchedDF_all), direction="both")
summary(nitrogen.lm.all.age)
# explain a bit more variation with age included: 0.148 vs 0.08543

final.model.stats[3,] <- c("carbon", "ngrip", "carbon.lm.all.additive", 
                           Regressionp(carbon.lm.all.additive), summary(carbon.lm.all.additive)$adj)
final.model.stats[4,] <- c("nitrogen", "ngrip", "nitrogen.lm.all.additive", 
                           Regressionp(nitrogen.lm.all.additive), summary(nitrogen.lm.all.additive)$adj)

# rename final models 
carbon.lm.final.ngrip <- carbon.lm.final
carbon.lm.clim.ngrip <- carbon.lm.clim
c.t.clim.ngrip <- c.t.clim
nitrogen.lm.final.ngrip <- nitrogen.lm.final
nitrogen.lm.clim.ngrip <- nitrogen.lm.clim
n.t.clim.ngrip <- n.t.clim

# remove variables 
rm("c.t.clim", "carbon.lm.all.additive", "carbon.lm.all.age", "carbon.lm.all.interaction", "carbon.lm.clim", "carbon.lm.final", "carbon.lm.taxon",
   "n.t.clim", "nitrogen.lm.all.additive", "nitrogen.lm.all.age", "nitrogen.lm.all.interaction", "nitrogen.lm.clim", "nitrogen.lm.final", "nitrogen.lm.taxon")


## Summarize results ----
print(final.model.stats)
summary(carbon.lm.final.ngrip)
summary(carbon.lm.final.hendy)
# for carbon, the hendy data better explains variation in carbon isotope values. this is seen for the full model, as well as the greater partial variation explained by d18O. In the ngrip models, the strength of Taxon relative to d18O increases in explaining variation in carbon isotopes.
summary(nitrogen.lm.final.ngrip)
summary(nitrogen.lm.final.hendy)
# for nitrogen, ngrip and hendy dat produce virtually identical and very weak models, both overall as well as in terms of the partial effect of d18O on carbon isotope variation.

# Based on this, move forward with Hendy dataset

## Final plots ####

  ### Figure 3 ----
  
  #grDevices::pdf("output/Figure3_lm_carbon_nitrogen_all_July2021_NF.pdf", width=8, height=8)
  grDevices::cairo_pdf("output/isotope paper final/Figure3_lm_carbon_nitrogen_all_Nov2021_JB.pdf", width=8, height=8)
  
  layout(matrix(seq(1:6), ncol=3, nrow=2, byrow=F), widths=c(2.5,2.5,1))
  par(mar=c(4,4,1,1), cex.axis=1, bty="l")
  
  
  # carbon
  plot(del13C_permil~d18O_hendy, data=matchedDF_all, pch=16, type="n", xlab="", ylab="")
  points(del13C_permil~d18O_hendy, 
         data=matchedDF_all[which(matchedDF_all$Taxon == "Sylvilagus"),], 
         pch=16, col="darkorange")
  points(del13C_permil~d18O_hendy, 
         data=matchedDF_all[which(matchedDF_all$Taxon == "Otospermophilus"),], 
         pch=16, col="royalblue2")
  abline(carbon.lm.final.hendy, lty=2)
  abline(carbon.lm.clim.hendy, lty=1)
  mtext(expression({delta}^18*O~'\u2030'), side=1, line=2.25) 
  mtext(expression({delta}^13*C~'\u2030'), side=2, line=2)
  
  boxplot(carbon.lm.clim.hendy$residuals ~ matchedDF_all$Taxon[which(!is.na(matchedDF_all$d18O_hendy))], 
          xlab="", ylab="")
  stripchart(carbon.lm.clim.hendy$residuals ~ matchedDF_all$Taxon[which(!is.na(matchedDF_all$d18O_hendy))], vertical=TRUE, add=TRUE, method="stack", col=c("royalblue2","darkorange"), pch=16)
  mtext("Taxon", side=1, line=2.25)
  mtext(expression('Residuals (Carbon Climate-only Model)'), side=2, line=2.25, cex=0.8)
  # mtext(expression('Residuals ('~{delta}^13*C~'\u2030'~' ~ '~{delta}^18*O~'\u2030'~')'), side=2, line=2.25, cex=0.8)
  legend("topright", legend=paste0("t=", round(c.t.clim.hendy$statistic,2), "; df=", round(c.t.clim.hendy$parameter,2), "; p=", round(c.t.clim.hendy$p.value,2)), bty = "n", cex = 0.8)
  
  # nitrogen
  plot(del15N_permil~d18O_hendy, data=matchedDF_all, pch=16, 
       xlab = "", ylab = "", type="n")
  points(del15N_permil~d18O_hendy, 
         data=matchedDF_all[which(matchedDF_all$Taxon == "Sylvilagus"),], 
         pch=16, col="darkorange")
  points(del15N_permil~d18O_hendy, 
         data=matchedDF_all[which(matchedDF_all$Taxon == "Otospermophilus"),], 
         pch=16, col="royalblue2")
  abline(nitrogen.lm.final.hendy, lty=2)
  abline(nitrogen.lm.clim.hendy, lty=1)
  mtext(expression({delta}^18*O~'\u2030'), side=1, line=2.25)
  
  mtext(expression({delta}^15*N~'\u2030'), side=2, line=2.25)
  
  boxplot(nitrogen.lm.clim.hendy$residuals ~ matchedDF_all$Taxon[which(!is.na(matchedDF_all$d18O_hendy))],
          xlab="", ylab="")
  stripchart(nitrogen.lm.clim.hendy$residuals ~ matchedDF_all$Taxon[which(!is.na(matchedDF_all$d18O_hendy))], vertical=TRUE, add=TRUE, method="stack", col=c("royalblue2","darkorange"), pch=16)
  mtext("Taxon", side=1, line=2.25)
  mtext(expression('Residuals (Nitrogen Climate-only Model)'), side=2, line=2.25, cex=0.8)
  # mtext(expression('Residuals ('~{delta}^15*N~'\u2030'~' ~ '~{delta}^18*O~'\u2030'~')'), side=2, line=2.25, cex=0.8)
  legend("topright", legend=paste0("t=", round(n.t.clim.hendy$statistic,2), "; df=", round(n.t.clim.hendy$parameter,2), "; p=", round(n.t.clim.hendy$p.value,2)), bty = "n", cex = 0.8)
  
  # taxon legend
  # Draw an empty plot
  plot(5, 5, 
       type="n", axes=FALSE, ann=FALSE, 
       xlim=c(0, 10), ylim = c(0,10))
  legend(xpd=T, x=-5, y=5, 
         legend = c("Otospermophilus", "Sylvilagus"),
         title="Taxon",
         col = c("royalblue2","darkorange"), pch = 16,
         bty = "n", cex = 0.8)
  
  # model legend
  plot(5, 5, 
       type="n", axes=FALSE, ann=FALSE, 
       xlim=c(0, 10), ylim = c(0,10))
  legend(xpd=T, x=-5, y=5,
         legend = c("Climate+Taxon", "Climate-only"),
         title="Model",
         lty = c(2,1),
         bty = "n", cex = 0.8)
  
  dev.off()
  
  # I haven't put A-D labels on them yet.
  # Figure caption text.
  # Figure 3. The relationship between isotope niche and climate. A) & C) show the relationship between 13C or 15N, respectively, and 18O.  In both panels, the dashed line indicates the fitted relationship between 13C or 15N and 18O from the final model which  includes Taxon as an independent variable. The solid line indicates the fitted relationship between 13C or 15N and 18O from a linear model that just includes 18O as the independent variable. B) & D) indicate the residuals from the climate-only model, plotted by taxon.
  
  ### Figure 4 ----
  #plot carbon and climate thru time
  
  # order matchedDF_all
  matchedDF_all <- matchedDF_all[order(matchedDF_all$median_age),]
  
  #grDevices::pdf(file="output/Figure4_carbon_climate_time_updated_NF.pdf", height=6, width=8)
  grDevices::cairo_pdf(file="output/isotope paper final/Figure4_carbon_climate_time_updated_Jan2022_JB.pdf", height=6, width=8)
  
  layout(matrix(seq(1:2), ncol=1, nrow=2), heights=c(0.75,1))
  par(mar=c(4,5,0,5), cex.axis=1, bty="l", xpd=F)
  
  # d13C plot
  plot(del13C_permil~median_age, data=matchedDF_all, 
       type="n", axes=FALSE, ann=FALSE, xaxs="i", yaxs="r", 
       xlim=c(55000, 0))
  axis(4)
  mtext(expression({delta}^13*C~'\u2030'), side=4, line=2.25)
  lines(del13C_permil~median_age, data=matchedDF_all, lty=1, col="black")
  points(del13C_permil~median_age, data=matchedDF_all[which(matchedDF_all$Taxon=="Otospermophilus"),], pch=16, col="royalblue2", cex=0.5)
  points(del13C_permil~median_age, data=matchedDF_all[which(matchedDF_all$Taxon=="Sylvilagus"),], pch=16, col="darkorange", cex=0.5)
  
  #d18O plot - Hendy
  plot(pach.d18O~HendyAge, dat=hendyDat,
       type="l",
       xlab="Years before present",
       ylab = expression({delta}^18*O~'\u2030'),
       bty="n",
       xlim=c(55000, 0), ylim=c(3, 0),
       lab=c(12, 8, 7), xaxs="i", yaxs="r", col="lightgray")
  x1.05<- loess(hendyDat$pach.d18O~hendyDat$HendyAge, span=0.05)
  lines(x1.05$fitted~x1.05$x, type="l", col="darkgray")

  points(d18O_hendy~median_age, data=matchedDF_all[which(matchedDF_all$Taxon=="Sylvilagus"),], pch=16, col="darkorange", cex=0.5)
  points(d18O_hendy~median_age, data=matchedDF_all[which(matchedDF_all$Taxon=="Otospermophilus"),], pch=16, col="royalblue2", cex=0.5)
  dev.off()



# Sensitivity analysis  ----
# sensitivity analysis #1: compare weighted and median age models, using the final models from the primary analysis

C_model_weighted <- lm(del13C_permil ~ specimen_wd18O + Taxon, data=matchedDF_all)
summary(C_model_weighted)
C_model_median <- lm(del13C_permil ~ d18O_hendy + Taxon, data=matchedDF_all)
summary(C_model_median) # this is the final model from the primary analyses

N_model_weighted <- lm(del15N_permil ~ specimen_wd18O + Taxon, data=matchedDF_all)
summary(N_model_weighted)
N_model_median <- lm(del15N_permil ~ d18O_hendy + Taxon, data=matchedDF_all)
summary(N_model_median) # this is the final model from the primary analyses

# compare weighted d180 vs d18O at median age
# No substantial difference. N still almost signif, C still highly signif. 
# Median is more significant than weighted, for what that's worth. Maybe only matters when we go to match isotopes with a specific age? weighted age and median age can be quite different


# How much of a difference does the variation in age make?
N=100 # Note: some specimens do not have 100 age estimates.
N_model_res <- as.data.frame(matrix(data=NA, nrow=N, ncol=6))
C_model_res<- as.data.frame(matrix(data=NA, nrow=N, ncol=6))
colnames(N_model_res) <- colnames(C_model_res) <- c('Fstat', 'lm_pVal', 'coeff', 'AdjR2', 'cor', 'cor_pVal')

for (k in 1:N){
  
  # for each specimen, sample a d18O value
  sampledd18O <- vector(mode="numeric", length=length(samples_final))
  
  for (i in 1:length(samples_final)){
    # pull out specimen ID
    spec <- samples_final[i] 
    # pull out all calibrated ages and probabilities
    specd18O <- allAges_d18O[which(allAges_d18O$name == spec),]
    # sample a single age, with sampling weighted per the probability
    sampledd18O[i] <- sample(specd18O$d18O, 1, prob=specd18O$probability)
  }
  
  # match with isotope data
  matchedDF_sensitivity <- cbind(isoDat2_final[,c('UCIAMS_Number', 'Taxon', 'del15N_permil', 'del13C_permil', 'X14C_age_BP')], sampledd18O)
  
  N_model_sensitivity <- lm(del15N_permil ~ sampledd18O, data=matchedDF_sensitivity)
  C_model_sensitivity <- lm(del13C_permil ~ sampledd18O, data=matchedDF_sensitivity)
  
  par(mfrow=c(1,2))
  plot(del15N_permil ~ sampledd18O, data=matchedDF_sensitivity, pch=16)
  abline(N_model_sensitivity)
  plot(del13C_permil ~ sampledd18O, data=matchedDF_sensitivity, pch=16)
  abline(C_model_sensitivity)
  
  N_cor.test_sensitivity <- cor.test(matchedDF_sensitivity$del15N_permil, matchedDF_sensitivity$sampledd18O)
  C_cor.test_sensitivity <- cor.test(matchedDF_sensitivity$del13C_permil, matchedDF_sensitivity$sampledd18O)
  
  N_model_res[k, 1]<- summary(N_model_sensitivity)$fstatistic[1]
  N_model_res[k, 2]<- summary(N_model_sensitivity)$coefficients[8]
  N_model_res[k, 3]<- summary(N_model_sensitivity)$coefficients[2]
  N_model_res[k, 4]<- summary(N_model_sensitivity)$adj.r.squared
  N_model_res[k, 5]<- N_cor.test_sensitivity$estimate
  N_model_res[k, 6]<- N_cor.test_sensitivity$p.value
  
  C_model_res[k, 1]<- summary(C_model_sensitivity)$fstatistic[1]
  C_model_res[k, 2]<- summary(C_model_sensitivity)$coefficients[8]
  C_model_res[k, 3]<- summary(C_model_sensitivity)$coefficients[2]
  C_model_res[k, 4]<- summary(C_model_sensitivity)$adj.r.squared
  C_model_res[k, 5]<- C_cor.test_sensitivity$estimate
  C_model_res[k, 6]<- C_cor.test_sensitivity$p.value
  
}

#Jessica - Can we run many iterations on the model above (e.g., 1000)
#and use the min of mins and max of maxes below to better quantify 
#the sensitivity of these estimates and the full range of R2 values?
# sure: that's what the script above does (though only for N=100), then the apply function below shows a variety of summary stats

# Calculate summary statistics on each column
apply(N_model_res, 2, summary)
apply(C_model_res, 2, summary)

# Supplemental figure for appendix - sensitivity analysis ----
## Compare final estimate with sensitivity analysis ----

# new supplemental figure, showing the test results for both Adj. R2 and the fitted coefficient for d18O.
# Note the following models to use:
# Carbon median ages univariate climate: carbon.lm.clim
# Nitrogen median ages univariate climate: nitrogen.lm.clim

# Carbon weighted age univariate climate
carbon.lm.clim.weighted <- lm(del13C_permil ~ specimen_wd18O, data=matchedDF_all)
# Nitrogen weighted age univariate climate
nitrogen.lm.clim.weighted <- lm(del15N_permil ~ specimen_wd18O, data=matchedDF_all)

# grDevices::pdf("output/SuppFigureX_climate_sensitivity_July2021_NF.pdf", width=8, height=6)
#Jessica - use line below instead of above
grDevices::cairo_pdf("output/SuppFigureX_climate_sensitivity_R2andcoeff.pdf", width=8, height=6)

par(mfrow=c(2,2)) 

# plot Carbon R2
y1<- hist(C_model_res$AdjR2, xlab=expression('Adjusted R'^2), main=expression('Model:'~{delta}^13*C~'\u2030'~'~'~{delta}^18*O~'\u2030'))
segments(summary(carbon.lm.clim.weighted)$adj.r.squared, 0, 
         summary(carbon.lm.clim.weighted)$adj.r.squared, max(y1$counts), 
         col="red", lwd=1)
segments(summary(carbon.lm.clim)$adj.r.squared, 0, 
         summary(carbon.lm.clim)$adj.r.squared, max(y1$counts), 
         col="blue", lwd=2)
legend(xpd=T, x=0.25, y=35, bty="n", legend=c("median d18O", "weighted d18O"), 
       col=c("blue", "red"), lwd=c(2,1), cex=0.5)

# plot Carbon d18O fitted coefficient
y2 <- hist(C_model_res$coeff, xlab=expression({delta}^18*O~'fitted coefficient'), main="")
segments(summary(carbon.lm.clim.weighted)$coefficients[2], 0, 
         summary(carbon.lm.clim.weighted)$coefficients[2], max(y2$counts), 
         col="red", lwd=1)
segments(summary(carbon.lm.clim)$coefficients[2], 0, 
         summary(carbon.lm.clim)$coefficients[2], max(y2$counts), 
         col="blue", lwd=2)


# plot Nitrogen R2
x1<- hist(N_model_res$AdjR2, xlab=expression('Adjusted R'^2), main=expression('Model:'~{delta}^15*N~'\u2030'~'~'~{delta}^18*O~'\u2030'))
segments(summary(nitrogen.lm.clim.weighted)$adj.r.squared, 0, 
         summary(nitrogen.lm.clim.weighted)$adj.r.squared, max(x1$counts), 
         col="red", lwd=1)
segments(summary(nitrogen.lm.clim)$adj.r.squared, 0, 
         summary(nitrogen.lm.clim)$adj.r.squared, max(x1$counts), 
         col="blue", lwd=2)

# plot Nitrogen d18O fitted coefficient
x2 <- hist(N_model_res$coeff, xlab=expression({delta}^18*O~'fitted coefficient'), main="")
segments(summary(nitrogen.lm.clim.weighted)$coefficients[2], 0, 
         summary(nitrogen.lm.clim.weighted)$coefficients[2], max(x2$counts), 
         col="red", lwd=1)
segments(summary(nitrogen.lm.clim)$coefficients[2], 0, 
         summary(nitrogen.lm.clim)$coefficients[2], max(x2$counts), 
         col="blue", lwd=2)
dev.off()

# old supplemental figure, showing the correlation test results
grDevices::pdf("output/SuppFigureX_climate_sensitivity_July2021_NF.pdf", width=8, height=6)
#Jessica - use commented out line below instead of above
#grDevices::cairo_pdf("output/SuppFigureX_climate_sensitivity_July2021_JB.pdf", width=8, height=6)
par(mfrow=c(1,2)) 
x<- hist(N_model_res$cor, xlab=expression('Correlation:'~{delta}^15*N~'\u2030'~'~'~{delta}^18*O~'\u2030'), main="")
segments(N_cor.test_weighted$estimate, 0, N_cor.test_weighted$estimate, max(x$counts), 
         col="red", lwd=2)
segments(N_cor.test_median$estimate, 0, N_cor.test_median$estimate, max(x$counts), 
         col="blue", lwd=2)
legend(xpd=T, x=0.25, y=35, bty="n", legend=c("median d18O", "weighted d18O"), 
       col=c("blue", "red"), lwd=c(2,1), cex=0.5)

y<- hist(C_model_res$cor, xlab=expression('Correlation:'~{delta}^13*C~'\u2030'~'~'~{delta}^18*O~'\u2030'), main="")
segments(C_cor.test_weighted$estimate, 0, C_cor.test_weighted$estimate, max(y$counts), 
         col="red", lwd=1)
segments(C_cor.test_median$estimate, 0, C_cor.test_median$estimate, max(y$counts), 
         col="blue", lwd=2)
dev.off()


# Final additional plot - age of Squirrels vs rabbits ----
age_ttest<- t.test(specimen_medianage~Taxon, data=matchedDF_all)

pdf(file="output/SuppFigY_taxon_age_distribution.pdf", height=6, width=6)
boxplot(specimen_medianage~Taxon, data=matchedDF_all, 
        ylim=c(55000, 0), 
        xlab="", ylab="Median Specimen Age (cal years BP)")
stripchart(specimen_medianage~Taxon, data=matchedDF_all, vertical=TRUE, ylim=c(55000, 0), add=TRUE, method="stack", col=c("royalblue2","darkorange"), pch=16)
legend("topright", legend=paste0("t=", round(age_ttest$statistic,2), "; df=", round(age_ttest$parameter,2), "; p=", round(age_ttest$p.value,2)), bty = "n", cex = 0.8)
dev.off()