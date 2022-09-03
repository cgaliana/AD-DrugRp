# CGaliana
# parse STITCH annotation files 
# Version: 1.0


# loading packages 

suppressPackageStartupMessages({
library(readr)
library(dplyr)})


# Set working directory
currentDataDir = getwd()
 

## loading files
chemicals_links <- read.csv("9606.protein_chemicals.links.v5.0.tsv", sep = "\t")
chemical_inchi <-read_delim("chemicals.inchikeys.v5.0.tsv", "\t", col_names = TRUE)
drugs_inchi <- read.csv("drugs_inchikey.csv")

## protein anotation from Biomart
protein_annot <-read.table("mart_export.txt", sep = "\t", header = TRUE) %>%rename(protein = Protein.stable.ID)

## save inchi_key
search_inchi <- drugs_inchi %>% .$inchi_key

## parse by inchi_key

drugs_inchi_CID <- drugs_inchi %>% left_join(chemical_inchi, by = "inchi_key")

table(is.na(drugs_inchi_CID$flat_chemical_id))

drugs_inchi_CID_links <- drugs_inchi_CID %>% 
	rename(chemical = flat_chemical_id) %>% 
  	left_join(chemicals_links, by = "chemical") %>%
  	mutate(across(everything(), gsub, pattern = "9606.", replacement = "")) %>%
  	left_join(protein_annot, by = "protein")

table(is.na(drugs_inchi_CID_links$Gene.name))

stitch_annot <- drugs_inchi_CID_links  %>% filter(combined_score >= 700 & !is.na(Gene.name)) %>% select(cmap_name, Gene.name)

table(is.na(stitch_annot$Gene.name))

write.csv(stitch_annot, "stitch_annot.csv")
