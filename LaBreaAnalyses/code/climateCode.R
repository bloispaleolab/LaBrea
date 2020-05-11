
# NOTE that this plot uses the 2003 data, whereas the next one (final) uses the 2002 data
dat<- read.delim("data/raw/climate/hendy2003data.txt")
plot(-smoothed~Y2Kage, data=dat, type="l", 
     xlab="Years before 2000", ylab="d18O benthic forams, smoothed",
     bty="n", xlim=c(75000, 10000),
     lab=c(12, 8, 7), xaxs="i", yaxs="r")

### This is the final plot, based on 2002 data
pdf(file="output/ClimateCurve-updated.pdf", height=6, width=8)
  dat2<- read.delim("data/raw/climate/hendy2002data.txt")
  plot(pach.d18O~HendyAge, dat=dat2, type="l", 
       xlab="Years before 2000", ylab="d18O Neogloboquadrina pachyderma",
       bty="n", xlim=c(60000, 0), ylim=c(3.1, 0.25), # modifiy to time bin relevant to P23. Original set from 75000 to 0 BP
       lab=c(12, 8, 7), xaxs="i", yaxs="r", col="gray")
  x1.05<- loess(dat2$pach.d18O~dat2$HendyAge, span=0.05)
  lines(x1.05$fitted~x1.05$x, type="l", col="red")
dev.off()


# Interpolate climate to calibrated radiocarbon ages of small mammals 

dat2<- read.delim("data/raw/climate/hendy2002data.txt")
#HendyAge_corr<-dat2$HendyAge -50 #correct for difference between RC age (before 1950) and Hendy Age (before 2000)
#dat2_corrected<-cbind(dat2, HendyAge_corr)

ages <- read.delim("data/processed/Interpolation/Calibrated_dates.txt") # mean, error, 

xout <- ages$Calibrated.age
Hendy_extracted <- approx(x=dat2$HendyAge, y=dat2$pach.d18O, method="linear", xout=xout)
Hendy_extracted$y # interpolated d18O at each radiocarbon date

write.csv(Hendy_extracted$y, file = "data/processed/Interpolation/interpolated_d18O.csv")


# weighted mean test
dat2<- read.delim("data/raw/climate/hendy2002data.txt")
ages <- read.delim("data/processed/Interpolation/P23_calibrated.txt") 
xout <- seq(29770, 48810, by=10)

seq_extracted <- approx(x=dat2$HendyAge, y=dat2$pach.d18O, method="linear", xout=xout)
seq_extracted$y
interpolated_seq<-cbind(seq_extracted$y, xout)
write.csv(interpolated_seq, file = "data/processed/Interpolation/interpolated_seq.csv")

