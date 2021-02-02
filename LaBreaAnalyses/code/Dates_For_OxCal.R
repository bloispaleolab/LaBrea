library('tidyverse')

# work with dates in radiocarbon space

# read in dates
dates_orig <- read.delim("data/processed/master_dates_file.txt",
                         sep="\t", stringsAsFactors = F)

# add 'processed' names ----
# read in mammal taxonomy file
mammalTaxa <- read.delim("data/raw/TaxonomyMatchingFile.txt", sep="\t")
matchedTaxa<- mammalTaxa[match(dates_orig$prelim_taxon_name, mammalTaxa$OriginalName),
                         c('RevisedName', 'Genus')]
dates <- cbind(dates_orig, matchedTaxa)

# Add column indicating a "too old" flag
GreaterThan <- rep(0, nrow(dates)) # set all to 0
GreaterThan[grep(">", dates$X14C_age_BP)] <- 1 # find the greater thans and re-set those to 1
dates <- cbind(dates, GreaterThan)
dates$X14C_age_BP <- gsub(">", "", dates$X14C_age_BP) #remove the > from the old dates and convert to numeric from factor
dates$X14C_age_BP <- as.numeric(dates$X14C_age_BP)

# Rename 14C
colnames(dates)[which(colnames(dates)=="X14C_age_BP")] <- "C14_age_BP"
colnames(dates)[which(colnames(dates)=="X14C_age_error")] <- "C14_age_error"

# re-order factor levels for plotting
## may not need to do this for OxCal script, but kept it in
dates$box <- factor(dates$box, levels = rev(c("14", "7b", "13", "1", "HC", "999", "4", "10")))
dates$RevisedName <- factor(dates$RevisedName, levels = rev(c("Sylvilagus sp", "Sylvilagus bachmani", "Sylvilagus audubonii", "Lepus sp", "Otospermophilus beecheyi", "Neotoma sp", "Thomomys sp", "Mustela frenata", "Canis latrans")))

# select only the columns I need
dates <- select(dates, "UCIAMS_Number", "Museum_Number", "C14_age_BP", "C14_age_error", "box", "Canister", "RevisedName", "Genus", "GreaterThan")


box <- c(1, 13, 14) # do these separate from 7b
for (i in 1:length(box)){
  box_dates <- dates[which(dates$box == box[i]),] #find all dates from a box

  # remove any "greater than" dates
  if (any(box_dates$GreaterThan == 1)){
    box_dates <- box_dates[-which(box_dates$GreaterThan == 1),]
  }

  # create the dataframe to store the dates
  dates_for_input <- as.data.frame(matrix(data=NA, nrow=nrow(box_dates), ncol=2))
  colnames(dates_for_input) <- c("type", "code")
  
  # first, create the basic dataframe for all dates
  dates_for_input$type <- rep(paste0("R_date"), nrow(box_dates))
  dates_for_input$code <- paste0("R_Date('UCIAMS ", box_dates$UCIAMS_Number, "',", box_dates$C14_age_BP, ",", box_dates$C14_age_error, ")")
  
  # then, if there are any duplicated dates, replace "R_date" with "R_combine"
  # find any dates duplicated across individuals. These should be entered into OxCal using the R_Combine function
  
  if (any(duplicated(box_dates$Museum_Number))){
    dup_rows <- which(duplicated(box_dates$Museum_Number))
    dup_MuseumNum <- unique(box_dates[dup_rows,'Museum_Number']) # can't just work with row numbers. need to get down to museum number in case multiple different specimens are re-dated, or there are three dates for an individual, etc.
    
    for (m in 1:length(dup_MuseumNum)){
      # which rows contain the first set of duplicated museum numbers?
      dups <-  which(box_dates$Museum_Number == dup_MuseumNum[m])
      dates_for_input$type[dups] <- rep(paste0("R_combine_",m), length(dups))
    }
  }
  
  write.csv(dates_for_input, file=paste0('output/OxCal/dates_for_input_box', box[i], '.csv'), row.names=F)
}

# now write script for 7b
box <- '7b'
for (i in 1:length(box)){
  box_dates <- dates[which(dates$box == box[i]),] #find all dates from a box
  
  # remove any "greater than" dates
  if (any(box_dates$GreaterThan == 1)){
    box_dates <- box_dates[-which(box_dates$GreaterThan == 1),]
  }
  
  # create the dataframe to store the dates
  dates_for_input <- as.data.frame(matrix(data=NA, nrow=nrow(box_dates), ncol=2))
  colnames(dates_for_input) <- c("type", "code")
  
  # first, create the basic dataframe for all dates
  dates_for_input$type <- rep(paste0("R_date"), nrow(box_dates))
  dates_for_input$code <- paste0("R_Date('UCIAMS ", box_dates$UCIAMS_Number, "',", box_dates$C14_age_BP, ",", box_dates$C14_age_error, ")")
  
  # then, if there are any duplicated dates, replace "R_date" with "R_combine"
  # find any dates duplicated across individuals. These should be entered into OxCal using the R_Combine function
  
  if (any(duplicated(box_dates$Museum_Number))){
    dup_rows <- which(duplicated(box_dates$Museum_Number))
    dup_MuseumNum <- unique(box_dates[dup_rows,'Museum_Number']) # can't just work with row numbers. need to get down to museum number in case multiple different specimens are re-dated, or there are three dates for an individual, etc.
    
    for (m in 1:length(dup_MuseumNum)){
      # which rows contain the first set of duplicated museum numbers?
      dups <-  which(box_dates$Museum_Number == dup_MuseumNum[m])
      dates_for_input$type[dups] <- rep(paste0("R_combine_",m), length(dups))
    }
  }
  
  write.csv(dates_for_input, file=paste0('output/OxCal/dates_for_input_box', box[i], '.csv'), row.names=F)
}
