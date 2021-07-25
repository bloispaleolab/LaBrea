## Final climate code. This will do the following:

# Section 1: Create Figure 2 of the community paper
# Section 2: Calculate climate values based on deposit age ranges

### Read in common data files
## things to do still - determine whether the -50 offset is appropriate

dat2<- read.delim("data/raw/climate/hendy2002data.txt")
ages <- read.csv("output/OxCal/final oxcal models/Table Y - DateSummary.csv", header=T)
#ages[,4:12] <- ages[,4:12] - 50 # subtract 50 years to all ages to make it years before 2000


# OLD NOTE: from everything I can tell, they do not necessarily use 'before 2000' in the 2002 paper, so commented out below!
#HendyAge_corr<-dat2$HendyAge -50 #correct for difference between RC age (before 1950) and Hendy Age (before 2000)
#dat2_corrected<-cbind(dat2, HendyAge_corr)

  ### SECTION 1 CODE ----
  
### This is the final plot, based on 2002 data

grDevices::cairo_pdf(file="output/Figure2_ClimateCurve-updated.pdf", height=6, width=8)
  plot(pach.d18O~HendyAge, dat=dat2, type="l", 
       xlab="Years before present", ylab = expression({delta}^18*O~'\u2030'),
       bty="n", 
       xlim=c(60000, 0), ylim=c(3, 0),
       lab=c(12, 8, 7), xaxs="i", yaxs="r", col="gray")
  x1.05<- loess(dat2$pach.d18O~dat2$HendyAge, span=0.05)
  lines(x1.05$fitted~x1.05$x, type="l", col="black")

  # Add age ranges
  segments(x0=ages$olderbound_median_age, x1=ages$youngerbound_median_age, y0=c(0, 0.5, 0.25, 0.75), y1=c(0, 0.5, 0.25, 0.75), col="darkgray")
  text(x=ages$youngerbound_median_age-1000, y=c(0, 0.5, 0.25, 0.75), labels=paste("Box", ages$P23_box[c(1,3,2,4)]), adj=0, cex=0.75)
  dev.off()

  ### SECTION 2 CODE ----
  
  # find max and min ages for overall ages across all 4 deposits
  max_age_95 <- max(ages$olderbound_median_age)
  min_age_95 <- min(ages$youngerbound_median_age)
  max_age_event <- max(ages$event_older_95_age)
  min_age_event <- min(ages$event_younger_95_age)
  overallMax <- max(max_age_95, max_age_event)
  overallMin <- min(min_age_95, min_age_event)
  
    xout <- seq(overallMax, overallMin, by=-10)
  Hendy_extracted <- approx(x=dat2$HendyAge, y=dat2$pach.d18O, method="linear", xout=xout)  # interpolated d18O at every 10 years

# calculate climate stats for each deposit - full 95% CI
overallOlder_indices <- match(ages$olderbound_median_age, xout)
overallYounger_indices  <- match(ages$youngerbound_median_age, xout)
overallAgeRange_d18O <- as.data.frame(matrix(data=NA, nrow=length(ages$P23_box), ncol=5), row.names=ages$P23_box)
colnames(overallAgeRange_d18O) <- c('mean', 'median', 'max', 'min', 'sd')
  
for (i in 1:nrow(ages)){
  overallAgeRange_d18O[i,'mean'] <- mean(Hendy_extracted$y[overallOlder_indices[i]:overallYounger_indices[i]])
  overallAgeRange_d18O[i,'median'] <- median(Hendy_extracted$y[overallOlder_indices[i]:overallYounger_indices[i]])
  overallAgeRange_d18O[i,'max'] <- max(Hendy_extracted$y[overallOlder_indices[i]:overallYounger_indices[i]])
  overallAgeRange_d18O[i,'min'] <- min(Hendy_extracted$y[overallOlder_indices[i]:overallYounger_indices[i]])
  overallAgeRange_d18O[i,'sd'] <- sd(Hendy_extracted$y[overallOlder_indices[i]:overallYounger_indices[i]])
}
  
overallAgeRange_d18O <- cbind(ages$olderbound_median_age, ages$youngerbound_median_age, overallAgeRange_d18O)
colnames(overallAgeRange_d18O)[1:2] <- c("OlderMedianAge", "YoungerMedianAge")

# calculate climate stats for each deposit - event 95% CI
eventOlder_indices <- match(ages$event_older_95_age, xout)
eventYounger_indices  <- match(ages$event_younger_95_age, xout)
eventAgeRange_d18O <- as.data.frame(matrix(data=NA, nrow=length(ages$P23_box), ncol=5), row.names=ages$P23_box)
colnames(eventAgeRange_d18O) <- c('mean', 'median', 'max', 'min', 'sd')

for (i in 1:nrow(ages)){
  eventAgeRange_d18O[i,'mean'] <- mean(Hendy_extracted$y[eventOlder_indices[i]:eventYounger_indices[i]])
  eventAgeRange_d18O[i,'median'] <- median(Hendy_extracted$y[eventOlder_indices[i]:eventYounger_indices[i]])
  eventAgeRange_d18O[i,'max'] <- max(Hendy_extracted$y[eventOlder_indices[i]:eventYounger_indices[i]])
  eventAgeRange_d18O[i,'min'] <- min(Hendy_extracted$y[eventOlder_indices[i]:eventYounger_indices[i]])
  eventAgeRange_d18O[i,'sd'] <- sd(Hendy_extracted$y[eventOlder_indices[i]:eventYounger_indices[i]])
}

eventAgeRange_d18O <- cbind(ages$event_older_95_age, ages$event_younger_95_age, eventAgeRange_d18O)
colnames(eventAgeRange_d18O)[1:2] <- c("OlderEventAge", "YoungerEventAge")

write.csv(overallAgeRange_d18O, file = "output/d18O_overallAgeRange.csv", row.names = FALSE)
write.csv(eventAgeRange_d18O, file = "output/d18O_eventAgeRange.csv", row.names = FALSE)

