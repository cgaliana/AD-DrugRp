# CGaliana
# drugs-target genes annotation
# Version: 1.0


# loading packages 

suppressPackageStartupMessages({
library(dplyr)
library(tidyr)})


# Set working directory
currentDataDir = getwd()
 

# loading files
## LINCS metadata 
compound_beta <- read.csv("compoundinfo_beta.txt", sep = "\t") %>% select(cmap_name, target, inchi_key)

## drugs from our library
drugs_include <- read.csv("/metadata_DrugsDB.csv") %>% select(drug) %>% rename(cmap_name = drug) %>% distinct()

## include target from Drug Repurposing Hub
rephub <- read.csv("repurposing_drugs_20200324.txt", sep = "\t") %>% separate_rows(target) %>% rename(cmap_name = pert_iname) %>% select(cmap_name, target) %>% distinct()

## include inchi_key/target/moa to our drug database from LINCS metadata files

drugs_include <- drugs_include %>% left_join(compound_beta, by = "cmap_name") %>% select(cmap_name, target) %>% distinct()
drugs_inchi <- drugs_include %>% select(cmap_name, inchi_key) %>% distinct()

write.csv(drugs_inchi, "/home/cristina/Escritorio/TFM/SCRIPTS/stitch/drugs_inchikey.csv")

drugs_targets <- drugs_include %>% select(cmap_name, target) %>% distinct()

drugs_targets <- rbind(drugs_targets, rephub) %>% distinct() 

## get target from STITCH

stitch <- read.csv("stitch_annot.csv") %>% rename(target = Gene.name) %>% select(cmap_name, target) %>% distinct()


drugs_targets <- rbind(drugs_targets, stitch) %>% 
  distinct() %>% 
  mutate_all(na_if,"") %>% 
  filter(!is.na(target)) %>% 
  group_by(cmap_name) %>% 
  mutate(targets = paste0(target, collapse = ";")) %>% 
  select(cmap_name, targets) %>% 
  distinct() 


table(is.na(drugs_targets$targets))

write.table(drugs_targets, file = "annot_targets_tostitch.txt", sep = "\t")

