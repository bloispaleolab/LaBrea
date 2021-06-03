
# import isotope data
isoDat <- read.delim('data/processed/master_isotopes_file.txt', header=T, sep="\t")
isoDat2 <- read.csv('data/processed/SIBER/SIBER_raw.csv', header=T)

# read in calibrated ages
allAges<- read.csv('output/OxCal/final oxcal models/AllAges_ages_probs.csv', header=T)

# match isotope file to ages file
samples <- paste("UCIAMS", isoDat2$UCIAMS_Number)
matches <- match(samples, allAges$name)

# remove the unmatched specimen from the files
isoDat2_final <- isoDat2[-which(is.na(matches)),]
samples_final <- samples[-which(is.na(matches))]

N=100
N_model_res <- as.data.frame(matrix(data=NA, nrow=N, ncol=6))
C_model_res<- as.data.frame(matrix(data=NA, nrow=N, ncol=6))
colnames(N_model_res) <- colnames(C_model_res) <- c('Fstat', 'lm_pVal', 'coeff', 'AdjR2', 'cor', 'cor_pVal')

for (k in 1:N){

# for each specimen, sample an age
sampledAges <- vector(mode="numeric", length=length(samples_final))

  for (i in 1:length(samples_final)){
    # pull out specimen ID
    spec <- samples_final[i] 
    # pull out all calibrated ages and probabilities
    specAges <- allAges[which(allAges$name == spec),]
    # sample a single age, with sampling weighted per the probability
    sampledAges[i] <- sample(specAges$value, 1, prob=specAges$probability)
  }
  
  # match with isotope data
  matchedDF <- cbind(isoDat2_final[,c('UCIAMS_Number', 'Taxon', 'del15N_permil', 'del13C_permil', 'X14C_age_BP')], sampledAges)
  
  N_model <- lm(del15N_permil ~ sampledAges, data=matchedDF)
  C_model <- lm(del13C_permil ~ sampledAges, data=matchedDF)
  
  par(mfrow=c(1,2))
  plot(del15N_permil ~ sampledAges, data=matchedDF, pch=16)
  abline(N_model)
  plot(del13C_permil ~ sampledAges, data=matchedDF, pch=16)
  abline(C_model)

  N_cor.test <- cor.test(matchedDF$del15N_permil, matchedDF$sampledAges)
  C_cor.test <- cor.test(matchedDF$del13C_permil, matchedDF$sampledAges)
  
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

hist(N_model_res$cor, xlim=c(-0.1,0.1))
hist(C_model_res$cor, xlim=c(0.3,0.45))


# Note for Nate
# this script does the sampling and then calculates the correlation with the ages
# but...you probably want to add a step where you match the sampledAges with the d18O climate curve, and then calculate the correlation