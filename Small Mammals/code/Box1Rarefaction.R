
# Read in file
dat <- read.csv("input/raw/P23_Dep1_SmallMammals.csv", header=T)

# create new species name that combines genus and species
dat$name <- paste(dat$Genus, dat$Species)

# for purposes of this analysis, combine different species into groups

all.species<- c(
  "Canis dirus", 
  "Canis latrans", 
  "Canis sp", 
  "Dipodomys sp", 
  "Mephitis mephitis", 
  "Microtus sp",
  "Mustela frenata",
  "Neotamias sp",
  "Neotamias merriami",
  "Neotoma sp",
  "Onychomys torridus",
  "Otospermophilus sp",
  "Otospermophilus beecheyi",
  "Perognathus sp.",
  "Peromyscus sp.",
  "Reithrodontomys megalotis",
  "Reithrodontomys sp",
  "Spilogale sp",
  "Sylvilagus sp",
  "Sylvilagus audubonii",
  "Sylvilagus bachmani",
  "Taxidea taxus",
  "Thomomys bottae",
  "Thomomys sp")

raw.names <- levels(as.factor(dat$name))
raw.names

# Replace orginal names with processed names
dat$name[which((dat$name == "Canis cf. dirus") | (dat$name == "Canis dirus"))] <- "Canis dirus"
dat$name[which((dat$name == "Canis cf. latrans") | (dat$name == "Canis latrans"))] <- "Canis latrans"
dat$name[which((dat$name == "Canis sp."))] <- "Canis sp"
dat$name[which((dat$name == "cf. Dipodomys sp.") | (dat$name == "Dipodomys sp."))] <- "Dipodomys sp"
dat$name[which((dat$name == "cf. Microtus ") | (dat$name == "cf. Microtus sp.") | (dat$name == "Microtus sp. ") | (dat$name == "Microtus ") | (dat$name == "Microtus sp."))] <- "Microtus sp" 
dat$name[which((dat$name == "Neotoma ") | (dat$name == "Neotoma sp."))] <- "Neotoma sp"
dat$name[which((dat$name == "Mustela cf. frenata"))] <- "Mustela frenata"
dat$name[which((dat$name == "Neotamias cf. merriami") | ("cf. Neotamias"))] <- "Neotamias sp"
dat$name[which((dat$name == "cf. Neotamias"))] <- "Otospermophilus sp"




## NATE: Finish this section

# Finish doing the name match, then write the file as a clean file
write.csv(dat, file="input/processed/P23_Dep1_SmallMammals-clean.csv")
