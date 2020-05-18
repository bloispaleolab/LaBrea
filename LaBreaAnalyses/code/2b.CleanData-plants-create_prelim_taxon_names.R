library(readr)
library(stringr)

# Read in data exported from Google Drive ----
deposits <- c("1","14","7b","13")
files <- list.files(
  "data/original_google_data/GoogleDriveExports-plants", 
  full=T)

# create a master spreadsheet with standardized taxonomic names ----
master <- NULL
for (i in 1:length(files)){
  # read in data ----
  original<- read_tsv(files[i], trim_ws=T, skip=1) 
  
  # standardize column names
  if (length(grep("Box 1 Loan 1", files[i])) == 1){
    colnames(original)[match(c('Gill temporary #', "identification"), colnames(original))] <- c("Lab_Number", "Identification")}else{
      colnames(original)[match(c('Number'), colnames(original))] <- c("Lab_Number")
    }
  
  colnames(original)[match(c('how many', "New Museum catalog # to be used"), colnames(original))] <-  
    c("NISP", 'Catalog_Number')
  
  
  # split out Box, grid, level and attach as new columns
  original$'Box, grid, level' <- sub('Box ', 'Split ', original$'Box, grid, level')
  original$'Box, grid, level' <- sub('grid ', 'Split ', original$'Box, grid, level')
  original$'Box, grid, level' <- sub('level ', 'Split ', original$'Box, grid, level')
  out <- strsplit(original$'Box, grid, level', "Split")
  out <- lapply(out, function(x){x[!x ==""]})
  out <- lapply(out, function(x){str_trim(x, side = "both")})
  out <- t(data.frame(out))
  rownames(out) <- seq(1:nrow(out))
  colnames(out) <- c("Box", "Grid", "Level")
  original <-cbind(original, out) 
  
  # keep the relevant columns and make sure in same order
  colsToKeep <- match(c("Lab_Number", "Catalog_Number", "Box", "Grid", "Level", "Identification", "NISP"), colnames(original))
  
  data <- original[,colsToKeep]
  
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
allSpecies <- unique(master$Identification)
allSpecies <- allSpecies[order(allSpecies)]

# export master file ----
write.table(master, file="data/processed/plant macros/master_plant_file.txt", sep="\t", row.names=F)
write.table(allSpecies, file="data/processed/plant macros/master_plant_taxa.txt", sep="\t", row.names=F, col.names=c("Original_Identification"))