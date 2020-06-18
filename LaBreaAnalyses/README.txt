# README

# Note, this is set up as an R project, so I always first double click on the 'LaBreaAnalyses.Rproj' file to open it in R. 
# Only after I've done that do I look at scripts, etc. This ensures that all the paths are correct and the scripts run properly


Our master workflow for mammals looks like this:
- Original data are of specimens.
- So each row in the google sheet is a single specimen.
	- *This isn’t entirely true because the data aren’t 100% clean, so we deal with this in the data cleanup step in R
- Each row in the dates or isotopes sheet is a measurement (date) or set of measurements (isotopes) from a single individual

Step 1: Make sure you have the most up-to-date version of the original data downloaded from Google Drive

- Github: LaBreaAnalyses/data/original_google_data
- I have downloaded different kinds of data as tsv files from Google Drive (mammals, plants), or for the dates/isotopes I have exported tab delimited text files of the excel sheets.


Step 2: Gather any auxiliary data you may need

- This is stored in the folder “raw” on Github: LaBreaAnalyses/data/raw
- This has files related to climate, traits, or our taxonomy matching files


Step 3: Code is stored on Github: LaBreaAnalyses/code

The main three scripts outline:
1. Downloading the data
2. (a, b, c) Cleaning the data 
            *For the specimens, this mainly deals with taxonomy, plus a few other things
3. NISP: Aggregating specimen-level data into summarized form
    1. We have created a separate “taxonomy matching file”.
    2. This takes the original taxonomic ID we assigned in the data cleaning step and matches it with a final ID for analysis
        1. Mainly, this pools some uncertain taxa together for purposes of analysis
        2. This taxonomy matching file is done offline -> we’ve taken the output from lines 129-132 in script 2A, copied that over into a column in excel/txt, then created a new column that indicates the “revised name”. We have also stored info about the Order, Family, and Genus in case we want to do our pooling at different levels (eg, pool all data to genus)


**In the end, there are two primary sets of files for the mammals:

1. The cleaned specimen level data file. This contains the master information about each specimen -
    1. Museum number
    2. Original taxonomic assignment
    3. Box and Canister info
    4. Main taxonomic ID (used for taxonomy matching later on)
2. Then, in the script “nisp.R”, we load the cleaned specimen level data, then aggregate it into a table with the following columns:
    1. RevisedName, Family, Genus, Box_1, Box_13, Box_14, Box_7b
    2. Each row now represents a unique taxon
    3. Values in the box columns represent the number of specimens of that taxon found in each box


