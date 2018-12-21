# Read in data exported from Google Drive

deposits <- c("1","17","7b","13")
files <- list.files("input/raw/GoogleDriveExports", full=T)

i <- 1
original<- read.delim(files[i], sep="\t")

# figure out prelim_taxon_name
rowsToSpecies <- which(original$Species != "")
rowsToGenus <- which(original$Genus != "")
rowsToSubfamily <- which(original$Subfamily == "")
