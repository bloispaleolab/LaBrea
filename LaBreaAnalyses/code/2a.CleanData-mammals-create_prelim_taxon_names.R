library(readr)

# Read in data exported from Google Drive ----
deposits <- c("1","7b","13","14", "misc_1", "misc_7b", "misc_13", "misc_14", "HC")
files <- list.files(
  "data/original_google_data/GoogleDriveExports-mammals", 
  full=T) #Remove "Old_versions" when error is resolved 
files <- files[-1] #when the "old versions" issue is dealt with, delete that folder and this line

# create a master spreadsheet with standardized taxonomic names ----
master <- NULL
for (i in 1:length(files)){
  # read in data ----
  original<- read_tsv(files[i], trim_ws=T) 
 
  # keep the relevant columns and make sure in same order
  colsToKeep <- match(c("Museum_Number", "UCM_Number", "Canister", "Class", "Order", "Family", "Subfamily", "Genus", "Species"), colnames(original))
  
  data <- original[,colsToKeep]

  # data cleaning ----
  # if sp. has a period, remove it!
  if (length(which(data$Species=="sp.")) > 0){
    data$Species[which(data$Species == "sp.")] <- "sp"
  }
  
  # replace cf. with cf
  if (length(grep("cf. ", data$Species)>0)){
    data$Species[grep("cf. ", data$Species)] <- gsub("cf. ", "cf ", data$Species[grep("cf. ", data$Species)])
  }
  if (length(grep("cf. ", data$Genus)>0)){
    data$Genus[grep("cf. ", data$Genus)] <- gsub("cf. ", "cf ", data$Genus[grep("cf. ", data$Genus)])
  }
  
  # figure out prelim_taxon_name ----
  allRows <- seq(1, nrow(data))
  rowsToSpecies <- intersect(which(data$Species != "sp"), which(data$Species != "")) #which rows have been identified to species?
  rowsToGenus <- intersect(which(data$Species == "sp"), which(data$Genus != "")) #which rows have been identified only to genus?
  rowsToSubfamily <- intersect(which(data$Subfamily != ""), which(data$Genus == ""))
  rowsToFamily <- intersect(which(data$Family != ""), which(data$Genus == ""))
  rowsToFamily <- rowsToFamily[-match(intersect(rowsToFamily, rowsToSubfamily), rowsToFamily)]
  otherRows <- allRows[-c(rowsToSpecies, rowsToGenus, rowsToSubfamily, rowsToFamily)]
  
  if (length(allRows) == length(otherRows) + length(rowsToSpecies) + length(rowsToGenus) + length(rowsToSubfamily) + length(rowsToFamily)){
    print(paste0(i, ": PROCEED: All rows accounted for"))
  }else{
    print(paste0(i, ": STOP: Not all rows accounted for"))
  }
  
  # assign preliminary taxon name
  data$prelim_taxon_name <- vector(length=nrow(data))
  data$prelim_taxon_name[rowsToSpecies] <- paste(data$Genus[rowsToSpecies], data$Species[rowsToSpecies], sep=" ")
  data$prelim_taxon_name[rowsToGenus] <- paste(data$Genus[rowsToGenus], data$Species[rowsToGenus], sep=" ")
  data$prelim_taxon_name[rowsToSubfamily] <- paste0(data$Family[rowsToSubfamily], " (", data$Subfamily[rowsToSubfamily], ")")
  data$prelim_taxon_name[rowsToFamily] <- as.character(data$Family[rowsToFamily])
  data$prelim_taxon_name[otherRows] <- paste(data$Class[otherRows], data$Order[otherRows], sep="-")
  
  # Add Box number to dataframe
  box <- sub('.*Deposit ', '', files[i])
  box <- sub(".tsv", '', box)
  data$box <- box 

  # Add Misc bones indicator dataframe
  if (length(grep("Misc", files[i]))>0) { 
    misc <- "y"
  }else{
    misc <- "n"}
  data$misc <- misc 

    
  # Add onto the master spreadsheet
  if (i ==1){
    master <- rbind(master, data)
    print(paste0(i, ": Rows added"))
  }else{
    if (all(colnames(data) == colnames(master))){
      master <- rbind(master, data)
      print(paste0(i, ": Rows added"))
    }else{
      print(paste0(i, ": CHECK COLNAMES"))
    }
  }
  
} 

#ERROR HERE WITH UPDATED DATAFILES - nothing unusual about updated datafile structure, just minor edits to "Misc Dep 1" and "Hancock"
# Jessica discussion: I don't see any issues with running this with the updated files. Is the error here, or in a different (downstream) script?

# deal with specimens with repeated catalog numbers ----

# first, clean up so Museum_Number so it matches. Most of the time, repeats separated with a semi-colon.
if (any(grep(":", master$Museum_Number))){
  master$Museum_Number[grep(":", master$Museum_Number)] <- gsub(":", ";", master$Museum_Number[grep(":", master$Museum_Number)])
}

# then find all rows with repeats
rowsWithRepeats <- grep(";", master$Museum_Number)

# scroll through each, copy row to end and separate catalog numbers
newRowsMaster <- NULL
for (i in 1:length(rowsWithRepeats)){
  # find row with repeats and split out catalog number
  oldRow <- master[rowsWithRepeats[i],]
  splitNumbers <-strsplit(as.character(oldRow$Museum_Number), "; ")[[1]] 
  
  # create new rows and assign individual catalog numbers
  l <- length(splitNumbers)
  newRows <- oldRow[1,]
  for (j in 2:l){
    newRows <- rbind(newRows, oldRow)
  }
  newRows$Museum_Number <- splitNumbers
  
  # add to master newRows dataframe
  newRowsMaster <- rbind(newRowsMaster, newRows)
  rm(oldRow, newRows)
}

# merge with master - delete old row, add new rows
# have to do this outside the loop, otherwise rownumbers get thrown off
master <- master[-rowsWithRepeats,]
master <- rbind(master, newRowsMaster)

# This block of code can be used for the taxonomy matching file
unique_names <- unique(master$prelim_taxon_name)
unique_names <- sort(unique_names)
unique_names

# replace "7B" with "7b"
master$box[which(master$box == "7B")] <- "7b"


# export master file ----
write.table(master, file="data/processed/master_mammal_file.txt", sep="\t")

