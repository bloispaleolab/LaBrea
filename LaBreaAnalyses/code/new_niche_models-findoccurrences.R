# species occurrence data ----
# retrieve and clean

library(spocc)

species_list <- "Peromyscus maniculatus"

# query selected database for occurrence records
i=1
results <- spocc::occ(query = species_list[i], from = "gbif", limit = 5000, has_coords = TRUE, gbifopts = list(basisOfRecord =))
# retrieve data table from spocc object
results.data <- results[["gbif"]]$data[[formatSpName("Peromyscus maniculatus")]]
# remove rows with duplicate coordinates
occs.dups <- duplicated(results.data[c('longitude', 'latitude')])
occs <- results.data[!occs.dups,]
# make sure latitude and longitude are numeric (sometimes they are characters)
occs$latitude <- as.numeric(occs$latitude)
occs$longitude <- as.numeric(occs$longitude)
# give all records a unique ID
occs$occID <- row.names(occs)

