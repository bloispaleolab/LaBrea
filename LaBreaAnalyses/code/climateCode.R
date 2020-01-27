
# NOTE that this plot uses the 2003 data, whereas the next one (final) uses the 2002 data
dat<- read.delim("data/original_data/climate/hendy2003data.txt")
plot(-smoothed~Y2Kage, data=dat, type="l", 
     xlab="Years before 2000", ylab="d18O benthic forams, smoothed",
     bty="n", xlim=c(75000, 10000),
     lab=c(12, 8, 7), xaxs="i", yaxs="r")

### This is the final plot, based on 2002 data
pdf(file="output/ClimateCurve-updated.pdf", height=6, width=8)
  dat2<- read.delim("data/original_data/climate/hendy2002data.txt")
  plot(pach.d18O~HendyAge, dat=dat2, type="l", 
       xlab="Years before 2000", ylab="d18O Neogloboquadrina pachyderma",
       bty="n", xlim=c(75000, 0), ylim=c(3.1, 0.25),
       lab=c(12, 8, 7), xaxs="i", yaxs="r", col="gray")
  x1.05<- loess(dat2$pach.d18O~dat2$HendyAge, span=0.05)
  lines(x1.05$fitted~x1.05$x, type="l", col="red")
dev.off()




