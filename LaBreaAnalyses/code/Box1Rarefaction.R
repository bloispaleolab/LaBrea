### Read in initial data ####

# Read in original file
dat <- read.csv("input/raw/P23_Dep1_SmallMammals.csv", header=T)

# Read in taxon matching file
taxonomy <- read.delim("input/raw/TaxonomyMatchingFile.txt", sep="\t", header=T)

### Match original taxa with revised names ####
# create new species name that combines genus and species
dat$name <- paste(dat$Genus, dat$Species)
unique(dat$name)

dat$name <- taxonomy[match(dat$name, taxonomy$OriginalName), 'RevisedName']

# Examine the revised data ####
head(dat)
tail(dat)
dat$name
unique(dat$name)

### Write the revised file to the processed folder ####
# Finish doing the name match, then write the file as a clean file
write.csv(dat, file="input/processed/P23_Dep1_SmallMammals-clean.csv")






### OLD CODE ####
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
  "Perognathus sp",
  "Peromyscus sp",
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
dat$name[which((dat$name == "Neotamias cf. merriami") | (dat$name == "cf. Neotamias"))] <- "Neotamias sp"
dat$name[which((dat$name == "cf. Otospermophilus"))] <- "Otospermophilus sp"
dat$name[which((dat$name == "Otospermophilus cf. beecheyi")| (dat$name == "cf. Otospermophilus beecheyi"))] <- "Otospermophilus beecheyi"
dat$name[which((dat$name == "Perognathus sp."))] <- "Perognathus sp"
dat$name[which((dat$name == "cf. Peromyscus") | (dat$name == "Peromyscus sp.") | (dat$name == "Peromyscus "))] <- "Peromyscus sp"
dat$name[which((dat$name == "Reithrodontomys cf. megalotis"))] <- "Reithrodontomys megalotis"
dat$name[which((dat$name == "cf. Reithrodontomys ") | (dat$name == "Reithrodontomys sp."))] <- "Reithrodontomys sp"
dat$name[which((dat$name == "cf. Spilogale "))] <- "Spilogale sp"
dat$name[which((dat$name == "cf. Sylvilagus ") | (dat$name == "cf. Sylvilagus sp.") | (dat$name == "Sylvilagus sp.") | (dat$name == "Sylvilagus "))] <- "Sylvilagus sp"
dat$name[which((dat$name == "Sylvilagus cf. audubonii"))] <- "Sylvilagus audubonii"
dat$name[which((dat$name == "Sylvilagus cf. bachmani") | (dat$name == "Sylvilagus cf. bachmanii"))] <- "Sylvilagus bachmani"
dat$name[which((dat$name == "cf. Taxidea taxus"))] <- "Taxidea taxus"
dat$name[which((dat$name == "Thomomys cf. bottae"))] <- "Thomomys bottae"
dat$name[which((dat$name == "cf. Thomomys ") | (dat$name == "Thomomys sp."))] <- "Thomomys sp"
