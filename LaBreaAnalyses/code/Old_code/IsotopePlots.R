# make isotope plots ----
library(tidyverse)

# read in files
iso <- read.delim("data/processed/master_isotopes_file.txt", sep="\t", stringsAsFactors = F)
dates <- read.delim("data/processed/master_dates_file.txt", sep="\t", stringsAsFactors = F)

# merge dates and isotopes ----

# first, fix suspect UCIAMS number in isotope data
iso[which(iso$UCIAMS_Number == 191244),'UCIAMS_Number'] <- 191242

# find matches based on UCIAMS number
matches <- match(iso$UCIAMS_Number, dates$UCIAMS_Number)

# isotopes without date matches - remove these rows from iso dataframe
iso <- iso[-which(is.na(matches)),]

# isotopes with date  matches - add dates to the iso dataframe
iso <- cbind(iso, dates[na.omit(matches),c('X14C_age_BP', 'X14C_age_error')])

# replace > dates with numbers
for (i in 1:nrow(iso)){
  if (length(grep(">", iso[i,'X14C_age_BP']))>0){
    iso[i,'X14C_age_BP'] <- 50000
  }
}

iso$X14C_age_BP <- as.numeric(iso$X14C_age_BP)
iso$box <- as.factor(iso$box)

#rename audubonii so it's standardized  # Delete this and add script below once plot is fixed
iso$prelim_taxon_name[c(grep("Sylvilagus audubonii", iso$prelim_taxon_name),
  grep("Sylvilagus cf. auduboni", iso$prelim_taxon_name),
  grep("Sylvilagus cf audubonii", iso$prelim_taxon_name))] <- "Sylvilagus audubonii"


##change lines 32-35 to this once output plot is fixed so all Sylvilagus are included, total iso_filtered specimens should = 82 but I can't get LACMP23-31423 to "grep"
#iso$prelim_taxon_name[c(grep("Sylvilagus audubonii", iso$prelim_taxon_name),
#grep("Sylvilagus cf audubonii", iso$prelim_taxon_name),
#grep("Sylvilagus ?audubonii", iso$prelim_taxon_name), #not grepping this specimen, LACMP23-31423, for some reason 
#grep("Sylvilagus bachmani", iso$prelim_taxon_name),
#grep("Sylvilagus sp", iso$prelim_taxon_name),
#grep("Sylvilagus cf bachmani", iso$prelim_taxon_name))] <- "Sylvilagus"

#Include all squirrels
iso$prelim_taxon_name[c(grep("Otospermophilus beecheyi", iso$prelim_taxon_name),
grep("Sciuridae", iso$prelim_taxon_name))] <- "Otospermophilus beecheyi"

iso$prelim_taxon_name <- as.factor(iso$prelim_taxon_name)

select(iso, Species = prelim_taxon_name, N15 = del15N_permil, C13 = del13C_permil, box, Age_14C= X14C_age_BP, Error_14C = X14C_age_error)

iso_filtered <- select(iso, Species = prelim_taxon_name, N15 = del15N_permil, C13 = del13C_permil, Box=box, Age_14C= X14C_age_BP, Error_14C = X14C_age_error) 

iso_filtered <- iso_filtered %>%
  filter(Species == "Sylvilagus audubonii" | Species == "Otospermophilus beecheyi")

long_iso <- iso_filtered %>%
  gather(Isotope, Value, N15:C13)

facet_names <- list(
  'C13'=expression({delta}^13*C),
  'N15'=expression({delta}^14*N)
  )
facet_labeller <- function(variable,value){
  return(facet_names[value])
}

matchingColors <-
  # Sylvilagus #00BE67
  # Otospermophilus #00BDD1
par <- ggplot(long_iso, aes(x=Age_14C,
                                y=Value,
                                color=Species,
                                shape=Box)) +
  geom_point() +
  geom_smooth(method="loess", span=0.75) + #Doesn't seem appropriate to force a linear fit on the plot since climate isnt linear
  scale_x_reverse() +
  scale_color_manual(values= c("Otospermophilus beecheyi"='#00BDD1', "Sylvilagus audubonii"='#00BE67')) +
  xlab("Age (14C)") +
  ylab(expression(~'\u2030')) +
  #theme_bw()
  theme_light()

pdf(file="output/isotopes.pdf", width=15, height=7.5, encoding="MacRoman") # 15, 7.5 #
par + 
  facet_grid(Isotope ~ ., scales="free", labeller=facet_labeller) +
  theme(strip.text.y = element_text(size = 12, colour = "black", angle = -90))
dev.off()
# current plot doesnt accomidate all 81-82 rabbit and squirrel specimens 
# Also, I'm confused why Otospermophilus gets two trend lines. Is each "greped" taxon treated as a unique group? 
# If so, would it be easier to just filter out the non-Sylvilagus/Otospermophilus taxa from "iso"?

  