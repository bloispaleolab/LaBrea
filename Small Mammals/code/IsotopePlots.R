# make isotope plots ----

# read in files
iso <- read.delim("data/processed/master_isotopes_file.txt", sep="\t")
dates <- read.delim("data/processed/master_dates_file.txt", sep="\t")

# merge dates and isotopes
iso$Museum_Number
dates$Museum_Number

match(isotopes$Museum_Number, mammals$Museum_Number)
