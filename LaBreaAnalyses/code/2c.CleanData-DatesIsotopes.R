# Match specimen numbers to radiocarbon dates and isotopes

# read in all data ----

# master mammal data 
mammals <- read.delim("data/processed/master_mammal_file.txt", sep="\t")
levels(mammals$Museum_Number)[grep("LAMMP23", levels(mammals$Museum_Number))] <- "LACMP23-30238"
mammals[grep("LAMMP23", mammals$Museum_Number),'Museum_Number'] <- "LACMP23-30238"

# radiocarbon dates and isotopes
files <- list.files(
  "data/original_google_data/GoogleDriveExports-dates_isotopes", 
  full=T)

dates_all <- read.delim(file=files[grep("Dates_Master", files)], sep="\t")
isotopes_all <- read.delim(file=files[grep("Isotopes_Master", files)], sep="\t")

mammal_dates <-dates_all[which(dates_all$SampleType=="mammal"),] 
mammal_isotopes <-isotopes_all[which(isotopes_all$SampleType=="mammal"),] 

plant_dates <-dates_all[which(dates_all$SampleType=="plant"),] 
plant_isotopes <-isotopes_all[which(isotopes_all$SampleType=="plant"),] 

# match specimens with dates - mammals ----


# dates with catalog number matches - add box and taxon to the date dataframe
tempD <- mammal_dates[which(!is.na(match(mammal_dates$Museum_Number, mammals$Museum_Number))),]
tempM <- mammals[na.omit(match(mammal_dates$Museum_Number, mammals$Museum_Number)), ]
if (all(as.character(tempD$Museum_Number) == as.character(tempM$Museum_Number))){
  tempDates <- cbind(tempD, tempM[,c('prelim_taxon_name', 'box', 'Canister')])
}else{
  print("STOP: specimens not matching!")
}

# dates without catalog number matches - add box and taxon to the date dataframe manually
### NEED TO FIX THIS CODE ONCE THE CATALOG NUMBER ISSUE IS RESOLVED! ####
tempD <- mammal_dates[which(is.na(match(mammal_dates$Museum_Number, mammals$Museum_Number))),]

prelim_taxon_name <- c("Canis latrans", "Canis latrans", "Otospermophilus beecheyi", "Otospermophilus beecheyi", "Otospermophilus beecheyi", "Otospermophilus beecheyi", "Sylvilagus sp", "Sylvilagus sp", "Sylvilagus sp", "Sylvilagus sp", "Otospermophilus beecheyi", "Otospermophilus beecheyi") # assigning the Spermpophilus to "Otospermophilus beecheyi" for now.

box <- c(1,1,1,1,1,1,4,4, 999, 999, 10, 999)
Canister <- c("B1/L3", "B1/L7", "B2/L8", "B1/L8", "B1/L4", "B2/L4", rep("NA", 6))

tempD <- cbind(tempD, prelim_taxon_name, box, Canister)
levels(tempDates$box) <- c(levels(tempDates$box),4,10,999)
tempDates <- rbind(tempDates, tempD)

# write dates to processed files
write.table(tempDates, file="data/processed/master_dates_file.txt", sep="\t")

# match specimens with isotopes - mammals ----

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
tempI <- mammal_isotopes[which(is.na(match(mammal_isotopes$Museum_Number, mammals$Museum_Number))),]
prelim_taxon_name <- c("Canis latrans", "Canis latrans", "Otospermophilus beecheyi", "Otospermophilus beecheyi", "Otospermophilus beecheyi", "Otospermophilus beecheyi", "Sylvilagus sp", "Sylvilagus sp", "Sylvilagus sp", "Sylvilagus sp", "Otospermophilus beecheyi", "Otospermophilus beecheyi") # assigning the Spermpophilus to "Otospermophilus beecheyi" for now.

box <- c(1,1,1,1,1,1,4,4, 999, 999, 10, 999)
Canister <- c("B1/L3", "B1/L7", "B2/L8", "B1/L8", "B1/L4", "B2/L4", rep("NA", 6))

tempI <- cbind(tempI, prelim_taxon_name, box)
levels(tempIsotopes$box) <- c(levels(tempIsotopes$box),4,10,999)
tempIsotopes <- rbind(tempIsotopes, tempI)

# write dates to processed files
write.table(tempIsotopes, file="data/processed/master_isotopes_file.txt", sep="\t")

