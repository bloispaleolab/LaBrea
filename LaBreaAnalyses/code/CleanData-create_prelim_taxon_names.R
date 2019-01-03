# Read in data exported from Google Drive ----
deposits <- c("1","17","7b","13")
files <- list.files(
  "/Users/jessicablois/Documents/GitHub/LaBrea/original_data/GoogleDriveExports-mammals", 
  full=T)

# create standardized taxonomic names ----
master <- NULL
for (i in 1:length(files)){
  original<- read.delim(files[i], sep="\t")
  
  original <- original[,-which(colnames(original)== "element")] # remove element column for now
  
  # Data cleaning
  if (i==1){
    original[414,'Family'] <- "" # specimen number entered into family field
  }
  
  if (i == 2 || i == 4){
    colnames(original)[1:2]<- c("UCM_Number", "Canister") # rename first two columns
  }
  
  if (i == 3){
    original <- original[,-c(13:15)] # remove unnecessary columns
  }
    
  # figure out prelim_taxon_name
  allRows <- seq(1, nrow(original))
  rowsToSpecies <- intersect(which(original$Species != "sp."),which(original$Species != ""))
  rowsToGenus <- intersect(which(original$Species == "sp."), which(original$Genus != ""))
  rowsToSubfamily <- intersect(which(original$Subfamily != ""), which(original$Genus == ""))
  rowsToFamily <- intersect(which(original$Family != ""), which(original$Genus == ""))
  rowsToFamily <- rowsToFamily[-match(intersect(rowsToFamily, rowsToSubfamily), rowsToFamily)]
  otherRows <- allRows[-c(rowsToSpecies, rowsToGenus, rowsToSubfamily, rowsToFamily)]
  
  if (length(allRows) == length(otherRows) + length(rowsToSpecies) + length(rowsToGenus) + length(rowsToSubfamily) + length(rowsToFamily)){
    print(paste0(i, ": PROCEED: All rows accounted for"))
  }else{
    print(paste0(i, ": STOP: Not all rows accounted for"))
  }
  
  # assign preliminary taxon name
  original$prelim_taxon_name[rowsToSpecies] <- paste(original$Genus[rowsToSpecies], original$Species[rowsToSpecies], sep=" ")
  original$prelim_taxon_name[rowsToGenus] <- paste(original$Genus[rowsToGenus], original$Species[rowsToGenus], sep=" ")
  original$prelim_taxon_name[rowsToSubfamily] <- paste0(original$Family[rowsToSubfamily], " (", original$Subfamily[rowsToSubfamily], ")")
  original$prelim_taxon_name[rowsToFamily] <- as.character(original$Family[rowsToFamily])
  original$prelim_taxon_name[otherRows] <- paste(original$Class[otherRows], original$Order[otherRows], sep="-")
  
  # Add Box number to dataframe
  box <- sub('.*Deposit ', '', files[i])
  box <- sub(".tsv", '', box)
  original$box <- box 
  
  # Add onto the master spreadsheet
  if (i ==1){
    master <- rbind(master, original)
    print(paste0(i, ": Rows added"))
  }else{
    if (all(colnames(original) == colnames(master))){
      master <- rbind(master, original)
      print(paste0(i, ": Rows added"))
    }else{
      print(paste0(i, ": CHECK COLNAMES"))
    }
  }
  
}

# deal with specimens with repeated catalog numbers ----

# first, clean up so Museum_Number so it matches
master$Museum_Number[grep(":", master$Museum_Number)] <- gsub(":", ";", master$Museum_Number[grep(":", master$Museum_Number)])

# then find all rows with repeats
rowsWithRepeats <- grep(";", master$Museum_Number)
newRowsMaster <- NULL

# scroll through each, copy row to end and separate catalog numbers
for (i in 1:length(rowsWithRepeats)){
  rm(oldRow, newRows)
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
}

# merge with master - delete old row, add new rows
# have to do this outside the loop, otherwise rownumbers get thrown off
master <- master[-rowsWithRepeats,]
master <- rbind(master, newRowsMaster)

# export master file ----
write.table(master, file="data/processed/master_mammal_file.txt", sep="\t")

