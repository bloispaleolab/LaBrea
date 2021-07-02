## Final climate code. This will do the following:

# Section 1: Create Figure 2 of the community paper
# Section 2: Calculate climate values based on deposit age ranges

### This is the final plot, based on 2002 data
dat2<- read.delim("data/raw/climate/hendy2002data.txt")
ages <- read.csv("output/OxCal/final oxcal models/Table Y - DateSummary.csv", header=T)
# add 50 years to all ages to make it years before 2000
ages[,4:12] <- ages[,4:12] + 50

grDevices::cairo_pdf(file="output/Figure2_ClimateCurve-updated.pdf", height=6, width=8)
  plot(pach.d18O~HendyAge, dat=dat2, type="l", 
       xlab="Years before 2000", ylab = expression({delta}^18*O~'\u2030'),
       bty="n", 
       xlim=c(60000, 0), ylim=c(3, 0),
       lab=c(12, 8, 7), xaxs="i", yaxs="r", col="gray")
  x1.05<- loess(dat2$pach.d18O~dat2$HendyAge, span=0.05)
  lines(x1.05$fitted~x1.05$x, type="l", col="black")

  # Add age ranges
  segments(x0=ages$olderbound_median_age, x1=ages$youngerbound_median_age, y0=c(0, 0.5, 0.25, 0.75), y1=c(0, 0.5, 0.25, 0.75), col="darkgray")
  text(x=ages$youngerbound_median_age-1000, y=c(0, 0.5, 0.25, 0.75), labels=paste("Box", ages$P23_box[c(1,3,2,4)]), adj=0, cex=0.75)
  dev.off()


## Final climate code ----
# Interpolate climate to calibrated radiocarbon ages of small mammals 

dat2<- read.delim("data/raw/climate/hendy2002data.txt")

# from everything I can tell, they do not necessarily use 'before 2000' in the 2002 paper, so commented out below!
#HendyAge_corr<-dat2$HendyAge -50 #correct for difference between RC age (before 1950) and Hendy Age (before 2000)
#dat2_corrected<-cbind(dat2, HendyAge_corr)

ages <- read.delim("data/processed/Interpolation/P23_calibrated.txt") # mean, error, 

xout <- ages$Median
Hendy_extracted <- approx(x=dat2$HendyAge, y=dat2$pach.d18O, method="linear", xout=xout)
Hendy_extracted$y # interpolated d18O at each radiocarbon date

matched_ages<-cbind(ages, Hendy_extracted$y)
colnames(matched_ages)[ncol(matched_ages)]<-"P23_d18O"
matched_ages$Deposit<-as.factor(matched_ages$Deposit)
d18O_mean<-aggregate(matched_ages$P23_d18O ~ matched_ages$Deposit, FUN=mean)
d18O_med<-aggregate(matched_ages$P23_d18O ~ matched_ages$Deposit, FUN=median)
Box_d18O_summary<-cbind(d18O_mean, d18O_med[,2])
colnames(Box_d18O_summary)<-c("Deposit", "P23_d18O_mean", "P23_d18O_med")
write.csv(matched_ages, file = "data/processed/Interpolation/matched_ages.csv", row.names = FALSE)
write.csv(Box_d18O_summary, file = "data/processed/Interpolation/Box_summary.csv", row.names = FALSE)

# weighted mean test
dat2<- read.delim("data/raw/climate/hendy2002data.txt")
ages <- read.delim("data/processed/Interpolation/P23_calibrated.txt") 
xout <- seq(29770, 48810, by=10)

seq_extracted <- approx(x=dat2$HendyAge, y=dat2$pach.d18O, method="linear", xout=xout)
seq_extracted$y
interpolated_seq<-cbind(seq_extracted$y, xout)
write.csv(interpolated_seq, file = "data/processed/Interpolation/interpolated_seq.csv")


# Add calibrated date boundaries to the plot ----
# boundary ages for each deposit calculated in OxCal
# Based on IntCal20

boundaries <- read.csv("output/OxCal/Boundaries.csv", header=T)
boundaries_no14 <- boundaries %>% filter(Deposit != 14)
OlderBound_mid <- (boundaries_no14$OlderBound_Old + boundaries_no14$OlderBound_Young)/2 
YoungerBound_mid <- (boundaries_no14$YoungerBound_Old + boundaries_no14$YoungerBound_Young)/2 

# dat2<- read.delim("data/raw/climate/hendy2002data.txt")
y0=y1=c(0.5, 0.75, 1)

pdf(file="output/ClimateCurve-updated.pdf", height=6, width=8)
plot(pach.d18O~HendyAge, dat=dat2, type="l", 
     xlab="Years before 2000", ylab = expression({delta}^18*O~'\u2030'),
     bty="n", xlim=c(65000, 0), ylim=c(3.1, 0.25), # modifiy to time bin relevant to P23. Original set from 75000 to 0 BP
     lab=c(12, 8, 7), xaxs="i", yaxs="r", col="gray")
x1.05<- loess(dat2$pach.d18O~dat2$HendyAge, span=0.05)
lines(x1.05$fitted~x1.05$x, type="l", col="red")


# add older boundary 95% range
#arrows(x0=boundaries_no14$OlderBound_Old, y0, x1=boundaries_no14$OlderBound_Young, y1, code=3, angle=90, length = 0.05, lwd=1.5, col="gray")

# add younger boundary 95% range
#arrows(x0=boundaries_no14$YoungerBound_Old, y0, x1=boundaries_no14$YoungerBound_Young, y1, code=3, angle=90, length = 0.05, lwd=1.5, col="gray")

# Add connector
#segments(x0=OlderBound_mid, y0, x1=YoungerBound_mid, y1, col=c(rgb(228/255,26/255,28/255), rgb(55/255,126/255,184/255), rgb(77/255,175/255,74/255)), lwd=3)

# This is the oldest bound of the oldest age to youngest bound of youngest age
segments(x0=boundaries_no14$OldestAge_Old, y0+.1, x1=boundaries_no14$YoungestAge_Young, y1+.1, col=c(rgb(228/255,26/255,28/255), rgb(55/255,126/255,184/255), rgb(77/255,175/255,74/255)), lwd=3)

dev.off()

cairo_pdf(file="output/ClimateCurve-updated-SPD.pdf", height=8, width=8)

layout(mat=as.matrix(c(1,2,3,4)), heights=c(2,2,2,10))
#layout.show(n=4)
par(mar=c(0,4.5,0,0))
plot(Box7b.spd, yaxt="n", xaxt="n", bty="n")#, col=rgb(77/255,175/255,74/255)) 
segments(x0=boundaries_no14$OldestAge_Old[3], 0, x1=boundaries_no14$YoungestAge_Young[3], 0, col=c(rgb(77/255,175/255,74/255)), lwd=6)

plot(Box13.spd, yaxt="n", xaxt="n", bty="n")#,col=rgb(55/255,126/255,184/255)) 
segments(x0=boundaries_no14$OldestAge_Old[2], 0, x1=boundaries_no14$YoungestAge_Young[2], 0, col=rgb(55/255,126/255,184/255), lwd=6)

plot(Box1.spd, yaxt="n", xaxt="n", bty="n")#, col=rgb(228/255,26/255,28/255)) 
segments(x0=boundaries_no14$OldestAge_Old[1], 0, x1=boundaries_no14$YoungestAge_Young[1], 0, col=rgb(228/255,26/255,28/255), lwd=6)

par(mar=c(4,4.5,0,0))
plot(pach.d18O~HendyAge, dat=dat2, type="l", 
     xlab="Years before 2000", ylab = expression({delta}^18*O~'\u2030'),
     bty="n", xlim=c(55000, 30000), ylim=c(3.1, 0.25), # modifiy to time bin relevant to P23. Original set from 75000 to 0 BP
     lab=c(12, 8, 7), xaxs="i", yaxs="r", col="gray")
x1.05<- loess(dat2$pach.d18O~dat2$HendyAge, span=0.05)
lines(x1.05$fitted~x1.05$x, type="l", col="red")

dev.off()

cairo_pdf(file="output/ClimateCurve-55-0.pdf", height=5, width=8)
par(mar=c(4,4.5,0,0))
plot(pach.d18O~HendyAge, dat=dat2, type="l", 
     xlab="Years before 2000", ylab = expression({delta}^18*O~'\u2030'),
     bty="n", xlim=c(55000, 0), ylim=c(3.1, 0.25), # modifiy to time bin relevant to P23. Original set from 75000 to 0 BP
     lab=c(12, 8, 7), xaxs="i", yaxs="r", col="gray")
x1.05<- loess(dat2$pach.d18O~dat2$HendyAge, span=0.05)
lines(x1.05$fitted~x1.05$x, type="l", col="red")
dev.off()
