# Match specimen numbers to radiocarbon dates and isotopes

# read in master mammal data 
mammals <- read.delim("data/processed/master_mammal_file.txt", sep="\t")

# read in radiocarbon dates and isotopes

files <- list.files(
  "/Users/jessicablois/Documents/GitHub/LaBrea/original_data/GoogleDriveExports-dates_isotopes", 
  full=T)
dates <- read.delim(file=files[grep("Dates_Master", files)], sep="\t")
isotopes <- read.delim(file=files[grep("Isotopes_Master", files)], sep="\t")

# find specimens
match(dates$Museum_Number, mammals$Museum_Number)
