# Convert specimen data from La Brea into a format suitable for entry into Neotoma database
# We need two files for each deposit:
# 1) Main data file, which should have the following columns:
# 2) specimen-level spreadsheet, with the following columns  

library(readr)
library(tidyr)

# This is Daniela's test edit 
# This is Daniela's second test 
# STEP 1: Read in data exported from Google Drive and clean in same manner as La Brea project----
deposits <- c("1","7b","13","14", "misc_1", "misc_7b", "misc_13", "misc_14", "HC")
files <- list.files(
  "data/GoogleDriveExports-mammals", 
  full=T) 

# create a master spreadsheet with standardized taxonomic names ----
master <- NULL
taxonomy_file <- read.delim("data/TaxonomyMatchingFile-forNeotomaDB.txt", sep="/t", header=T)

for (i in 1:length(files)){
  
    # read in data ----
  original<- read_tsv(files[i], trim_ws=T) 
  
  # keep the relevant columns and make sure in same order
  colsToKeep <- match(c("Museum_Number", "UCM_Number", "Canister", "Class", "Order", "Family", "Subfamily", "Genus", "Species", "Element_original"), 
                      colnames(original))
  
  data <- original[,colsToKeep]
  
  # data cleaning ----
  # if sp. has a period, remove it!
  if (length(which(data$Species=="sp.")) > 0){
    data$Species[which(data$Species == "sp.")] <- "sp"
  }
  
  # fix terminology for Peromyscus cf. californicus
  if (length(which(data$Species=="cf P. californicus")) > 0){
    data$Species[which(data$Species == "cf P. californicus")] <- "cf californicus"
  }
  
  # replace cf with cf.
  if (length(grep("cf ", data$Species)>0)){
    data$Species[grep("cf ", data$Species)] <- gsub("cf ", "cf. ", data$Species[grep("cf ", data$Species)])
  }
  if (length(grep("cf ", data$Genus)>0)){
    data$Genus[grep("cf ", data$Genus)] <- gsub("cf ", "cf. ", data$Genus[grep("cf ", data$Genus)])
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
  
  # 2nd step taxon name cleaning ### JESSICA MAKE SURE THIS WORKS
  # replace the prelim_taxon_name with the correct name in taxonomy_file$OriginalName_cleaned
  data[,'prelim_taxon_name'] <- taxonomy_file[match(taxonomy_file$OriginalName_cleaned, data$prelim_taxon_name),'OriginalName_cleaned']
  
  # delete the rows with NA taxa # JESSICA CHECK THIS TOO!!
  if (length(which(data$prelim_taxon_name=="NA-NA")) > 0){
    data <- data[-which(data$prelim_taxon_name=="NA-NA"),]
  }
  
  # Add Box number to dataframe
  box <- sub('.*Deposit ', '', files[i])
  box <- sub(".tsv", '', box)
  if (length(grep("Hancock", files[i]))>0) {
    box <- "HC"}
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


# clean up taxonomy for Tilia



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

# replace "7B" with "7b"
master$box[which(master$box == "7B")] <- "7b"

# export master file ----
# all specimens
write.table(master, file="data/output/master_mammal_file - all specimens combined.txt", sep="\t", row.names = F)

# STEP 2: Break cleaned data back into spreadsheets for each box ----
master <- read.delim(file="data/output/master_mammal_file - all specimens combined.txt", sep="\t")

deposits <- unique(master$box)

for (k in 1:length(deposits)){
  
  # First, create the Tilia main 'Data' table ----
  # 1) Main data file, which should have the following columns:
  # Name --> Taxon ID
  # Element	 --> "bone/tooth"
  # Units --> "present/absent"
  # Taphonomy	> blank
  # Group	> blank
  # Analysis units --> add one column per analysis unit, which corresponds to canister
  dat <- master[which(master$box == deposits[k]),]
  
  # find unique analysis units
  AnUnits <- unique(dat$Canister)
  AnUnits <- na.omit(AnUnits)
  
  taxa <- unique(dat$prelim_taxon_name)
  
  master_data <- matrix(nrow = length(taxa), ncol=length(AnUnits))
  master_data <- as.data.frame(master_data)
  colnames(master_data) <- AnUnits
  
  # assign presences to the main Analysis Units/Canisters
  for (i in 1:length(AnUnits)){
    temp <- dat[which(dat$Canister == AnUnits[i]),]
    master_data[match(unique(temp$prelim_taxon_name), taxa), which(colnames(master_data) == AnUnits[i])] <- 1
    rm(temp)
  }
  
  # add on misc bones
  # these should only be bones marked 'y' for misc, but without an assigned canister
  if (length(intersect(which(dat$misc=="y"), which(is.na(dat$Canister))))>0){
    master_data$misc <- NA
    temp <- dat[intersect(which(dat$misc=="y"), which(is.na(dat$Canister))),]
    master_data[match(unique(temp$prelim_taxon_name), taxa), 'misc'] <- 1
    rm(temp)
  }
  
  # check - all Analysis Units have data?
  cat("Deposit", as.character(deposits[k]))
  print(any(colSums(master_data, na.rm=T)==0)) # should be FALSE
  
  # delete NAs 
  master_data[is.na(master_data)] <- ""
  
  # add on other info for Tilia
  Name <- taxa
  Element <- rep("bone/tooth", length(taxa))
  Units <- rep("present/absent", length(taxa))  
  
  # reorder columns:
  master_data <- cbind(Name, Element, Units, master_data) 
  
  # all specimens
  write.table(master_data, file=paste0("data/output/Tilia_data-Box", deposits[k], ".txt"), sep="\t", row.names = F)
 
  ## JESSICA NOTES: ----
  ## 1) Do you want to add another set of rows with the NISP data too?
  ## 2) We still need to do some taxonomy cleanup here - match to lookup file?
  
  # Second, arrange the specimen data for easy input into the specimens table of Tilia ----
  # 2) specimen-level spreadsheet, with the following columns  
  # Spec ID --> Museum_Number
  # Depth
  # Anal Unit --> Canister
  # Taxon
  # element --> Element_original
  # symmetry
  # portion
  # maturity
  # sex
  # domestic status
  # taphonomy
  # preservative
  # nisp --> 1
  # repository --> "Los Angeles County Museum of Natural History" (or "LACM")
  # Spec Nr --> Museum_Number 
  # Field Nr
  # Arctos Nr
  # GenBank Nr
  # Notes
  
  # First, replace the "NA" with "misc" in the "Canister"
  master_specimens <- dat
  master_specimens[which(is.na(master_specimens$Canister)), 'Canister'] <- "misc"
  
  
  # start to deal with elements
  #symmetry:
  symmetry <- NULL
  symmetry[grep("rt", master_specimens$Element_original)] <- "right"
  symmetry[grep("lt", master_specimens$Element_original)] <- "left"
}
