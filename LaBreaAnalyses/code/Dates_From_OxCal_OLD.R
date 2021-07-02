# load the final calibrated dates and deposit chronologies and get them all in order for input 


boxes <- c('1', '7b', '13', '14')
all.dates <- matrix(data=NA, nrow=0, ncol=8)
all.dates <- as.data.frame(all.dates)

for (i in 1:length(boxes)){
  file <- paste0('output/OxCal/final oxcal models/box ', boxes[i], '/Box', boxes[i], '_final_event.csv')

  headers = read.csv(file, skip = 3, header = F, nrows = 1, as.is = T)
  dat = read.csv(file, skip = 6, header = F)
  colnames(dat)= headers
  colnames(dat)[1] <- 'Name'
  dat <- dat[-grep('Phase', dat$Name),]
  head(dat)
  tail(dat)
  
  # find rows with dates
  boundary.rows <- grep("Boundary", dat$Name)
  event.rows <- grep("Event", dat$Name)
  ignore.rows <- c(boundary.rows,event.rows)
  ignore.rows
  
  dates <- dat[-ignore.rows,c(1, 14:19)]
  dates$box <- boxes[i]
  colnames(all.dates) <- colnames(dates)
  all.dates <- rbind(all.dates, dates)
}

write.table(all.dates, file = 'output/OxCal/final oxcal models/calibrated.dates.txt', sep="\t", row.names=F)
