library(readr)

# Read in data exported from Google Drive ----
deposits <- c("1","14","7b","13")
files <- list.files(
  "data/original_data/GoogleDriveExports-plants", 
  full=T)

# create a master spreadsheet with standardized taxonomic names ----
master <- NULL
for (i in 1:length(files)){
  # read in data ----
  original<- read_tsv(files[i], trim_ws=T) 
  
  # change column names
  if (length(grep("Box 1 Loan 1", files[i])) == 1){
    colnames(original)[match(c('Gill temporary #', "identification", 'how many'), colnames(original))] <- c("Number", "Identification", "NISP")}else{
      colnames(original)[match(c('how many'), colnames(original))] <- c("NISP")
    }
  
  # keep the relevant columns and make sure in same order
  colsToKeep <- match(c("Number", "Identification", "NISP"), colnames(original))
  
  data <- original[,colsToKeep]
  
  # change names
  colnames(data) <- c("SpecimenNumber", "Species", "NISP")
  
  # Add Box number to dataframe
  box <- sub('.*Box ', '', files[i])
  if (length(grep("Loan", box)) == 1){
    box <- gsub(' Loan.*', '', box)}else{
    box <- gsub(".tsv", '', box)
  }

  data$box <- box
  
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

# data cleaning ----
allSpecies <- unique(master$Species)
allSpecies <- allSpecies[order(allSpecies)]



# merge with master - delete old row, add new rows
# have to do this outside the loop, otherwise rownumbers get thrown off
master <- master[-rowsWithRepeats,]
master <- rbind(master, newRowsMaster)


# This block of code can be used for the taxonomy matching file
unique_names <- unique(master$prelim_taxon_name)
unique_names <- sort(unique_names)
unique_names

# export master file ----
write.table(master, file="data/processed/master_mammal_file.txt", sep="\t")



### OLD CODE



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

