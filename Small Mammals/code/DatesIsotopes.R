# Match specimen numbers to radiocarbon dates and isotopes

# read in all data ----

# master mammal data 
mammals <- read.delim("data/processed/master_mammal_file.txt", sep="\t")

# radiocarbon dates and isotopes
files <- list.files(
  "/Users/jessicablois/Documents/GitHub/LaBrea/original_data/GoogleDriveExports-dates_isotopes", 
  full=T)
dates <- read.delim(file=files[grep("Dates_Master", files)], sep="\t")
isotopes <- read.delim(file=files[grep("Isotopes_Master", files)], sep="\t")

# remove extra Museum_Number column
dates <- dates[, -match('Museum_Number.1', colnames(dates))]

# match specimens with dates ----

# dates with catalog number matches - add box and taxon to the date dataframe
tempD <- dates[which(!is.na(match(dates$Museum_Number, mammals$Museum_Number))),]
tempM <- mammals[na.omit(match(dates$Museum_Number, mammals$Museum_Number)), ]
if (all(as.character(tempD$Museum_Number) == as.character(tempM$Museum_Number))){
  tempDates <- cbind(tempD, tempM[,c('prelim_taxon_name', 'box')])
}else{
  print("STOP: specimens not matching!")
}

# dates without catalog number matches - add box and taxon to the date dataframe manually
### NEED TO FIX THIS CODE ONCE THE CATALOG NUMBER ISSUE IS RESOLVED! ####
tempD <- dates[which(is.na(match(dates$Museum_Number, mammals$Museum_Number))),]
prelim_taxon_name <- c("Sylvilagus sp.", "Sylvilagus sp.", "Canis latrans", "Canis latrans", "Otospermophilus beecheyi", "Otospermophilus beecheyi", "Otospermophilus beecheyi", "Otospermophilus beecheyi", "Sciuridae", "Sciuridae")
box <- c(1,1,1,1,1,1,1,1,14,14)

tempD <- cbind(tempD, prelim_taxon_name, box)
tempDates <- rbind(tempDates, tempD)

# write dates to processed files
write.table(tempDates, file="data/processed/master_dates_file.txt", sep="\t")

# match specimens with isotopes ----

# isotopes with catalog number matches - add box and taxon to the date dataframe
tempI <- isotopes[which(!is.na(match(isotopes$Museum_Number, mammals$Museum_Number))),]
tempM <- mammals[na.omit(match(isotopes$Museum_Number, mammals$Museum_Number)), ]
if (all(as.character(tempI$Museum_Number) == as.character(tempM$Museum_Number))){
  tempIsotopes <- cbind(tempI, tempM[,c('prelim_taxon_name', 'box')])
}else{
  print("STOP: specimens not matching!")
}

# dates without catalog number matches - add box and taxon to the date dataframe manually
### NEED TO FIX THIS CODE ONCE THE CATALOG NUMBER ISSUE IS RESOLVED! ####
tempI <- isotopes[which(is.na(match(isotopes$Museum_Number, mammals$Museum_Number))),]
prelim_taxon_name <- c("Sylvilagus sp.", "Sylvilagus sp.", "Canis latrans", "Canis latrans", "Otospermophilus beecheyi", "Otospermophilus beecheyi", "Otospermophilus beecheyi", "Otospermophilus beecheyi", "Sciuridae", "Sciuridae", "Sciuridae")
box <- c(1,1,1,1,1,1,1,1,14,14,14)

tempI <- cbind(tempI, prelim_taxon_name, box)
tempIsotopes <- rbind(tempIsotopes, tempI)

# write dates to processed files
write.table(tempIsotopes, file="data/processed/master_isotopes_file.txt", sep="\t")

