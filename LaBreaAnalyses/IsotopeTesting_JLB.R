library(dplyr)
library(stringr)

# Import and Clean data ----
## import isotope data ----
# isoDat <- read.delim('data/processed/master_isotopes_file.txt', header=T, sep="\t") # this is not used in script
isoDat2 <- read.csv('data/processed/SIBER/SIBER_raw.csv', header=T)
# match isotope file to ages file
samples <- paste("UCIAMS", isoDat2$UCIAMS_Number) # this is the final set of samples with isotope data

## read in climate Data ----
climDat<- read.delim("data/raw/climate/hendy2002data.txt")

## read in calibrated ages ----
allAges<- read.csv('output/OxCal/final oxcal models/AllAges_ages_probs.csv', header=T)

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

#remove specimens without stable isotope values
allAges<- allAges[!(allAges$name=="UCIAMS 223515" | 
                     allAges$name=="UCIAMS 223519" | allAges$name=="UCIAMS 223522"),]

#remove specimens with unmeasureable dates
allAges<- allAges[!(allAges$name=="UCIAMS 223515" | 
                     allAges$name=="UCIAMS 223519" | allAges$name=="UCIAMS 223522"),]


# Match d180 to all age estimates for every sample ----
# end result is a new dataframe, adding d18O to allAges
all(allAges$value<0) # this should be true, then take the absolute values to match to climate
xout <- abs(allAges$value)
Hendy_extracted <- approx(x=climDat$HendyAge, y=climDat$pach.d18O, method="linear", xout=xout)
d18O<-Hendy_extracted$y
allAges_d18O<-cbind(allAges, d18O) # this matches d180 to all the age estimates, for all probabilities

# remove any unmatched specimen from the files
matches <- match(samples, allAges_d18O$name)

# added ifelse statement to catch instances where there are or are not matches
if (length(which(is.na(matches)))>0){
  isoDat2_final <- isoDat2[-which(is.na(matches)),]
  samples_final <- samples[-which(is.na(matches))]
}else{
  isoDat2_final <- isoDat2
  samples_final <- samples
}
# end result dataframes:
# isoDat2_final - final C/N isotope data
# samples_final - final set of samples
# allAges_d18O - ages, probabilities, and d18O for each sample


# Primary analysis: calculate weighted mean d18O and weighted age for each specimen ----
# For each specimen, calculate the weighted d180 and weighted mean age, based on the probability of each calibrated age from the age models, then do the linear models

specimen_wd18O <- vector(mode="numeric", length=length(samples_final))
specimen_wage <- vector(mode="numeric", length=length(samples_final))

for (i in 1:length(samples_final)){
  # pull out specimen ID
  spec <- samples_final[i] 
  # pull out all calibrated ages and probabilities for that single specimen
  specd18O <- allAges_d18O[which(allAges_d18O$name == spec),]
  # calculate weighted mean d18O
  specimen_wd18O[i] <- weighted.mean(specd18O$d18O, specd18O$probability)
  # calculate weighted mean age
  specimen_wage[i] <- weighted.mean(specd18O$value, specd18O$probability)
}

# match with isotope data
matchedDF_weighted <- cbind(isoDat2_final[,c('UCIAMS_Number', 'Taxon', 'del15N_permil', 'del13C_permil', 'X14C_age_BP')], specimen_wd18O, specimen_wage)

N_model_weighted <- lm(del15N_permil ~ specimen_wd18O, data=matchedDF_weighted)
summary(N_model_weighted)
C_model_weighted <- lm(del13C_permil ~ specimen_wd18O, data=matchedDF_weighted)
summary(C_model_weighted)

par(mfrow=c(1,2))
plot(del15N_permil ~ specimen_wd18O, data=matchedDF_weighted, pch=16)
abline(N_model_weighted)
plot(del13C_permil ~ specimen_wd18O, data=matchedDF_weighted, pch=16)
abline(C_model_weighted)

N_cor.test_weighted <- cor.test(matchedDF_weighted$del15N_permil, matchedDF_weighted$specimen_wd18O)
C_cor.test_weighted <- cor.test(matchedDF_weighted$del13C_permil, matchedDF_weighted$specimen_wd18O)


# Sensitivity analysis  ----
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
  #matchedDF <- cbind(matchedDF, weighted.d18O) ## NATE - No 'weighted.d18O' object found, commented this out.
  
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

par(mfrow=c(1,2))
hist(N_model_res$cor, xlab="N~O correlation")
segments(N_cor.test_weighted$estimate, 0, N_cor.test_weighted$estimate, 17, 
         col="red", lwd=2)
hist(C_model_res$cor, xlab="C~O correlation")
segments(C_cor.test_weighted$estimate, 0, C_cor.test_weighted$estimate, 25, 
         col="red", lwd=2)


# Original analysis, naive mean, no weighting ____
# extract mean d18O value across all ages for per specimen average - code needs to be fixed
mean.d18O<-aggregate( d18O ~ name, allAges_d18O, mean )

#Remove "UCIAMS" from name column 
mean.d18O <- mean.d18O %>%
  mutate_at("name", str_replace, "UCIAMS", "")

#Get rid of extra spaces
mean.d18O[,1:2] <- lapply(mean.d18O[,1:2], trimws)

#Convert d18O back to numeric values 
mean.d18O2<-transform(mean.d18O, d18O = as.numeric(d18O))

# Sort "matchedDF" specimens in assending order like "mean.d18O" specimens
matchedDF2 <- matchedDF[order(matchedDF$UCIAMS_Number),]

# combine isotope and climate files for analysis
d18O.d13C<-cbind(matchedDF2, mean.d18O2)

#Remove extra spaces in taxon levels
d18O.d13C$Taxon <- trimws(d18O.d13C$Taxon, which = c("right"))
d18O.d13C<-transform(d18O.d13C, Taxon = as.factor(Taxon))

####plot data with sample-averaged d18O values
##(these are the data I think we should present in the paper)

##Note - these are Jessica's plots copied from "Iso.Clim.R" with updated data

# Methods potential text:
# For carbon and nitrogen, I fit a linear model that originally included oxygen, taxon, and the interaction between the two variables as independent variables. I then performed stepwise regression to determine a final model. 
# Results potential text
# Stepwise regression indicated that there was no significant interaction between 18O and taxon for either carbon or nitrogen stable isotope values. For carbon, variation in 13C was significantly associated with both 18O and taxon (stats from summary(carbon.lm.final)). For nitrogen, neither 18O nor taxon explained significant variation in 15N, though taxon as a single variable was marginally significant (stats from summary(nitrogen.lm.taxon)).


## models - Carbon

# model with interaction term included
carbon.lm.all.interaction<-lm(del13C_permil~d18O * Taxon, data=d18O.d13C)
summary(carbon.lm.all.interaction)

# model with no interaction term 
# NOTE: this is what the final model is in the end, so THIS  IS WHAT YOU SHOULD REPORT IN THE PAPER.
carbon.lm.all<-lm(del13C_permil~d18O + Taxon, data=d18O.d13C)
summary(carbon.lm.all)

# climate-only model
carbon.lm.clim<-lm(del13C_permil~d18O, d18O.d13C)
summary(carbon.lm.clim)
plot(del13C_permil ~ d18O, data=d18O.d13C, pch=16)
abline(carbon.lm.clim)

# taxon-only model
carbon.lm.taxon<-lm(del13C_permil~Taxon, data=d18O.d13C)
summary(carbon.lm.taxon)

##NOTE - interaction between d13C and climate is weaker, and d13C and taxon
#is stronger, than with IntCal13 data. Variation is now ~equally explained 
#by climate and taxon, rather than dominated by climate - update paper
#discussion accordingly 

#stepwise regression 
carbon.lm.final <- step(lm(del13C_permil~d18O*Taxon, data=d18O.d13C, direction="both"))
summary(carbon.lm.final)
# --> this shows that the d180 + Taxon model (no interaction) is the best final model. 

# t-test of residuals from climate only model
c.t.clim <- t.test(carbon.lm.clim$residuals ~ d18O.d13C$Taxon[which(!is.na(d18O.d13C$d18O))])

# models - Nitrogen

# model with interaction term included
nitrogen.lm.all.interaction<-lm(del15N_permil ~ d18O * Taxon, data=d18O.d13C)
summary(nitrogen.lm.all.interaction)

# model with no interaction term 
nitrogen.lm.all.additive<-lm(del15N_permil ~ d18O + Taxon, data=d18O.d13C)
summary(nitrogen.lm.all.additive)

# climate-only model
nitrogen.lm.clim<-lm(del15N_permil ~ d18O, data=d18O.d13C)
summary(nitrogen.lm.clim)

# taxon-only model
nitrogen.lm.taxon<-lm(del15N_permil~Taxon, data=d18O.d13C)
summary(nitrogen.lm.taxon)

#stepwise regression --> this shows that there is not a good final model
nitrogen.lm.final <- step(lm(del15N_permil~d18O*Taxon, data=d18O.d13C, direction="both"))
# --> there is not a good final model, so I am just going to treat nitrogen the same as carbon for plotting.
nitrogen.lm.final <- nitrogen.lm.all.additive

# t-test of residuals from climate-only model
n.t.clim <- t.test(nitrogen.lm.clim$residuals ~ d18O.d13C$Taxon[which(!is.na(d18O.d13C$d18O))])


# JESSICA's NEW final plot ####

pdf("output/Isoplots/lm_carbon_nitrogen_all.pdf", width=8, height=8)

par(mfcol=c(2,2), mar=c(4,4,1,1), cex.axis=0.8)

# carbon
plot(del13C_permil~d18O, data=d18O.d13C, pch=16, type="n", xlab="", ylab="")
points(del13C_permil~d18O, 
       data=d18O.d13C[which(d18O.d13C$Taxon == "Sylvilagus"),], 
       pch=16, col="darkorange")
points(del13C_permil~d18O, 
       data=d18O.d13C[which(d18O.d13C$Taxon == "Otospermophilus"),], 
       pch=16, col="royalblue2")
abline(carbon.lm.final, lty=2)
abline(lm(del13C_permil~d18O, data=d18O.d13C), lty=1)
mtext(expression({delta}^18*O~'\u2030'), side=1, line=2.25) 
mtext(expression({delta}^13*C~'\u2030'), side=2, line=2) 
legend("bottomleft", legend = c("Otospermophilus", "Sylvilagus"),
       col = c("royalblue2","darkorange"), pch = 16, 
       bty = "n", cex = 0.8)

boxplot(carbon.lm.clim$residuals ~ d18O.d13C$Taxon[which(!is.na(d18O.d13C$d18O))], 
        xlab="", ylab="")
stripchart(carbon.lm.clim$residuals ~ d18O.d13C$Taxon[which(!is.na(d18O.d13C$d18O))], vertical=TRUE, add=TRUE, method="stack", col=c("royalblue2","darkorange"), pch=16)
mtext("Taxon", side=1, line=2.25)
mtext(expression('Residuals ('~{delta}^13*C~'\u2030'~' ~ '~{delta}^18*O~'\u2030'~')'), side=2, line=2.25, cex=0.8)
legend("topright", legend=paste0("t=", round(c.t.clim$statistic,2), "; df=", round(c.t.clim$parameter,2), "; p=", round(c.t.clim$p.value,2)), bty = "n", cex = 0.8)


# nitrogen
plot(del15N_permil~d18O, data=d18O.d13C, pch=16, 
     xlab = "", ylab = "", type="n")
points(del15N_permil~d18O, 
       data=d18O.d13C[which(d18O.d13C$Taxon == "Sylvilagus"),], 
       pch=16, col="darkorange")
points(del15N_permil~d18O, 
       data=d18O.d13C[which(d18O.d13C$Taxon == "Otospermophilus"),], 
       pch=16, col="royalblue2")
abline(nitrogen.lm.final, lty=2)
abline(lm(del15N_permil~d18O, data=d18O.d13C), lty=1)
legend("topleft", legend = c("Otospermophilus", "Sylvilagus"),
       col = c("royalblue2","darkorange"), pch = 16, 
       bty = "n", cex = 0.8)
mtext(expression({delta}^18*O~'\u2030'), side=1, line=2.25)
mtext(expression({delta}^15*N~'\u2030'), side=2, line=2.25)

boxplot(nitrogen.lm.clim$residuals ~ d18O.d13C$Taxon[which(!is.na(d18O.d13C$d18O))],
        xlab="", ylab="")
stripchart(nitrogen.lm.clim$residuals ~ d18O.d13C$Taxon[which(!is.na(d18O.d13C$d18O))], vertical=TRUE, add=TRUE, method="stack", col=c("royalblue2","darkorange"), pch=16)
mtext("Taxon", side=1, line=2.25)
mtext(expression('Residuals ('~{delta}^15*N~'\u2030'~' ~ '~{delta}^18*O~'\u2030'~')'), side=2, line=2.25, cex=0.8)
legend("topright", legend=paste0("t=", round(n.t.clim$statistic,2), "; df=", round(n.t.clim$parameter,2), "; p=", round(n.t.clim$p.value,2)), bty = "n", cex = 0.8)

dev.off()

# I haven't put A-D labels on them yet.
# Figure caption text.
# Figure 4. The relationship between isotope niche and climate. A) & C) show the relationship between 13C or 15N, respectively, and 18O.  In both panels, the dashed line indicates the fitted relationship between 13C or 15N and 18O from the final model which  includes Taxon as an independent variable. The solid line indicates the fitted relationship between 13C or 15N and 18O from a linear model that just includes 18O as the independent variable. B) & D) indicate the residuals from the climate-only model, plotted by taxon.


#plot carbon and climate thru time

ylim.prim <- c(0, 3)   
ylim.sec <- c(-23, -17)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

pdf("output/Isoplots/carbon_climate_time_updated.pdf", width=7, height=4)
C<-ggplot(d18O.d13C, aes(X14C_age_BP, d18O)) +
  geom_line(aes(y = a + del13C_permil*b), color = "brown") +
  geom_line(aes(y = d18O), color = "blue") +
  scale_y_reverse(name = expression({delta}^18*O~'\u2030'), sec.axis = sec_axis(~ (. - a)/b, name = expression({delta}^13*C~'\u2030'))) +
  scale_x_reverse("Radiocarbon Years Before Present")+
  theme_light() 
C + theme(axis.text.y.left = element_text(color="blue"),
          axis.text.y.right = element_text(color="brown"))
dev.off()
#This fig is not very informative with the averaged d18O data, perhaps omit
# from the paper?


