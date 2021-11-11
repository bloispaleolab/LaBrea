# compare Hendy and ngrip records

hendy<- read.delim("data/raw/climate/hendy2002data.txt")
ngrip<- read.delim("data/raw/climate/ngrip.txt")

# clean data

#ngrip dates are "before 2000", so need to subtract 50 years
ngrip$Age <- ngrip$Age-50
ngrip <- ngrip[-which(ngrip$Age>65000),]
ngrip <- ngrip[-which(ngrip$Age<0),]

# some NAs in hendy d18O
hendy <- hendy[-which(is.na(hendy$pach.d18O)),]

ageOlder <- 65000
ageYounger <- 0

# calculate anomalies for plotting
ngrip_z <- (ngrip$d18O - mean(ngrip$d18O)) / sd(ngrip$d18O)
hendy_z <- (hendy$pach.d18O - mean(hendy$pach.d18O)) / sd(hendy$pach.d18O)
hendy_z <- -hendy_z #get it on same direction as ngrip


plot(ngrip$Age, ngrip_z, type="l", col="black", lty=1, xlim=c(65000,0), ylim=c(-2.2, 3.6), xlab="Age (years BP)", ylab="d18O z-score")
lines(hendy$HendyAge, hendy_z, type="l", col="red")

#calculate anomaly between two records
xout <- seq(65000, 0, by=-100)
Hendy_extracted <- approx(x=hendy$HendyAge, y=hendy_z, method="linear", xout=xout)
ngrip_extracted <- approx(x=ngrip$Age, y=ngrip_z, method="linear", xout=xout)
climDiff <- Hendy_extracted$y - ngrip_extracted$y

plot(climDiff ~ xout, type="l", col="black", lty=1, xlim=c(65000,0))
segments(65000, 0, 0,0)

plot(Hendy_extracted$y~ngrip_extracted$y, pch=16)
