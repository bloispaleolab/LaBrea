library(tidyverse)
library(stringr)

# function to arrange dates for OxCal ####
# after running this function, 
# paste into the BBEdit file 'oxcal_input_file.txt', 
# then copy that into oxcal
format_for_oxcal <- function(file){
  text <- NULL
  for (i in 1:nrow(file)){
    text <- c(text, cat(
      'R_Date("', 
      file$UCIAMS_Number[i], 
      '", ',
      file$C14_age_BP[i], 
      ", ", 
      file$C14_age_error[i],
      ");", sep=""))
  }
  return(text)
}


# filter dates in various ways ----
dates <- read.delim(file="data/processed/master_dates_file.txt", sep="\t", stringsAsFactors = F)

# filter out "no" dates and Box 999 dates
dates <- dates[-which(dates$UseSample=="N"),]
dates <- dates[-which(dates$box=="999"),]

# create separate files for each box
box1_all <- dates %>%
  filter(box == 1)
box14_all <- dates %>%
  filter(box == 14)
box7b_all <- dates %>%
  filter(box == "7b")
box13_all <- dates %>%
  filter(box == 13)

#box1_squirrels <- dates %>%
  #filter(box == 1, prelim_taxon_name == "Otospermophilus beecheyi")

#box1_rabbits <- dates %>%
  #filter(box == 1) %>%
  #filter(str_detect(tolower(prelim_taxon_name), pattern = "sylvilagus"))

#box1_S.audubonii <- dates %>%
  #filter(box == 1, prelim_taxon_name == "Sylvilagus audubonii")

format_for_oxcal(box1_all)


