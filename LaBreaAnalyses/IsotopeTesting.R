
# import isotope data
isoDat <- read.delim('data/processed/master_isotopes_file.txt', header=T, sep="\t")
isoDat2 <- read.csv('data/processed/SIBER/SIBER_raw.csv', header=T)

# read in calibrated ages
allAges<- read.csv('output/OxCal/final oxcal models/AllAges_ages_probs.csv', header=T)

# match isotope file to ages file
samples <- paste("UCIAMS", isoDat2$UCIAMS_Number)

# Start of test script - commented out for now through line 76----
# matches <- match(samples, allAges$name)
# 
# # remove the unmatched specimen from the files
# # added ifelse statement to catch instances where there are or are not matches
# if (length(which(is.na(matches)))>0){
#   isoDat2_final <- isoDat2[-which(is.na(matches)),]
#   samples_final <- samples[-which(is.na(matches))]
# }else{
#   isoDat2_final <- isoDat2
#   samples_final <- samples
# }
# 
# N=100
# N_model_res <- as.data.frame(matrix(data=NA, nrow=N, ncol=6))
# C_model_res<- as.data.frame(matrix(data=NA, nrow=N, ncol=6))
# colnames(N_model_res) <- colnames(C_model_res) <- c('Fstat', 'lm_pVal', 'coeff', 'AdjR2', 'cor', 'cor_pVal')
# 
# for (k in 1:N){
# 
# # for each specimen, sample an age
# sampledAges <- vector(mode="numeric", length=length(samples_final))
# 
#   for (i in 1:length(samples_final)){
#     # pull out specimen ID
#     spec <- samples_final[i] 
#     # pull out all calibrated ages and probabilities
#     specAges <- allAges[which(allAges$name == spec),]
#     # sample a single age, with sampling weighted per the probability
#     sampledAges[i] <- sample(specAges$value, 1, prob=specAges$probability)
#   }
#   
#   # match with isotope data
#   matchedDF <- cbind(isoDat2_final[,c('UCIAMS_Number', 'Taxon', 'del15N_permil', 'del13C_permil', 'X14C_age_BP')], sampledAges)
#   
#   N_model <- lm(del15N_permil ~ sampledAges, data=matchedDF)
#   C_model <- lm(del13C_permil ~ sampledAges, data=matchedDF)
#   
#   par(mfrow=c(1,2))
#   plot(del15N_permil ~ sampledAges, data=matchedDF, pch=16)
#   abline(N_model)
#   plot(del13C_permil ~ sampledAges, data=matchedDF, pch=16)
#   abline(C_model)
# 
#   N_cor.test <- cor.test(matchedDF$del15N_permil, matchedDF$sampledAges)
#   C_cor.test <- cor.test(matchedDF$del13C_permil, matchedDF$sampledAges)
#   
#   N_model_res[k, 1]<- summary(N_model)$fstatistic[1]
#   N_model_res[k, 2]<- summary(N_model)$coefficients[8]
#   N_model_res[k, 3]<- summary(N_model)$coefficients[2]
#   N_model_res[k, 4]<- summary(N_model)$adj.r.squared
#   N_model_res[k, 5]<- N_cor.test$estimate
#   N_model_res[k, 6]<- N_cor.test$p.value
#   
#   C_model_res[k, 1]<- summary(C_model)$fstatistic[1]
#   C_model_res[k, 2]<- summary(C_model)$coefficients[8]
#   C_model_res[k, 3]<- summary(C_model)$coefficients[2]
#   C_model_res[k, 4]<- summary(C_model)$adj.r.squared
#   C_model_res[k, 5]<- C_cor.test$estimate
#   C_model_res[k, 6]<- C_cor.test$p.value
#   
# }
# 
# hist(N_model_res$cor, xlim=c(-0.1,0.1))
# hist(C_model_res$cor, xlim=c(0.3,0.45))


# Start of Nate's script. This is the final one that is being modified ----
# this script does the sampling and then calculates the correlation with the ages
# but...you probably want to add a step where you match the sampledAges with the d18O climate curve, and then calculate the correlation

dat2<- read.delim("data/raw/climate/hendy2002data.txt")

xout <- abs(allAges$value)
Hendy_extracted <- approx(x=dat2$HendyAge, y=dat2$pach.d18O, method="linear", xout=xout)
d18O<-Hendy_extracted$y
allAges_d18O<-cbind(allAges, d18O)

matches <- match(samples, allAges_d18O$name)

# remove the unmatched specimen from the files
# added ifelse statement to catch instances where there are or are not matches
if (length(which(is.na(matches)))>0){
  isoDat2_final <- isoDat2[-which(is.na(matches)),]
  samples_final <- samples[-which(is.na(matches))]
}else{
  isoDat2_final <- isoDat2
  samples_final <- samples
}

N=100
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
  matchedDF <- cbind(isoDat2_final[,c('UCIAMS_Number', 'Taxon', 'del15N_permil', 'del13C_permil', 'X14C_age_BP')], sampledd18O)
  #matchedDF <- cbind(matchedDF, weighted.d18O) ## NATE - No 'weighted.d18O' object found, commented this out.
  
  N_model <- lm(del15N_permil ~ sampledd18O, data=matchedDF)
  C_model <- lm(del13C_permil ~ sampledd18O, data=matchedDF)
  
  par(mfrow=c(1,2))
  plot(del15N_permil ~ sampledd18O, data=matchedDF, pch=16)
  abline(N_model)
  plot(del13C_permil ~ sampledd18O, data=matchedDF, pch=16)
  abline(C_model)
  
  N_cor.test <- cor.test(matchedDF$del15N_permil, matchedDF$sampledd18O)
  C_cor.test <- cor.test(matchedDF$del13C_permil, matchedDF$sampledd18O)
  
  N_model_res[k, 1]<- summary(N_model)$fstatistic[1]
  N_model_res[k, 2]<- summary(N_model)$coefficients[8]
  N_model_res[k, 3]<- summary(N_model)$coefficients[2]
  N_model_res[k, 4]<- summary(N_model)$adj.r.squared
  N_model_res[k, 5]<- N_cor.test$estimate
  N_model_res[k, 6]<- N_cor.test$p.value
  
  C_model_res[k, 1]<- summary(C_model)$fstatistic[1]
  C_model_res[k, 2]<- summary(C_model)$coefficients[8]
  C_model_res[k, 3]<- summary(C_model)$coefficients[2]
  C_model_res[k, 4]<- summary(C_model)$adj.r.squared
  C_model_res[k, 5]<- C_cor.test$estimate
  C_model_res[k, 6]<- C_cor.test$p.value
  
}

hist(N_model_res$cor)
hist(C_model_res$cor)

min(C_model_res$AdjR2)
max(C_model_res$AdjR2)

# these plots are going back to just one of the sampled times. Commented out.
# d18O_N<-lm(del15N_permil ~ sampledd18O, data=matchedDF)
# summary(d18O_N)
# plot(del15N_permil ~ sampledd18O, data=matchedDF, pch=16)
# abline(d18O_N)
# 
# d18O_C.lm<-lm(del13C_permil ~ sampledd18O, data=matchedDF)
# summary(d18O_C.lm)
# plot(del13C_permil ~ sampledd18O, data=matchedDF, pch=16)
# abline(d18O_C.lm)


#extract mean d18O value across all ages for per specimen average - code needs to be fixed
mean.d18O<-aggregate(allAges_d18O[,8], list(allAges_d18O$name), mean) 
#
names<-paste("UCIAMS", matchedDF$UCIAMS_Number)
match(names, weighted.d18O$Group.1) # these numbers don't make sense to me

