# CGaliana
# parse XML DrugBank Database (The latest release of DrugBank Online (version 5.1.9, released 2022-01-03)
# Version: 1.0


#----------------------------install.packages---------------------------------------------------#
#install.packages("XML", dependencies = TRUE)

#-------------------------------loading packages-------------------------------------------------# 

suppressPackageStartupMessages({
library(XML)
library(purrr)
library(tidyverse)})


# Set working directory
currentDataDir = getwd()
 

#loading & reading the xml
url <- "full database.xml"
xmldoc <-xmlParse(url)
rootnode <-xmlRoot(xmldoc)

# parent node: <drug>
# getting root's children
all_drugs <- xmlChildren(rootnode)  


## -------------------------------functions to parse------------------------------#

## getting the basic info: drugbank-id and other identifiers that may be associated with the drug.
getDinfo <- function(drug) {    
    tibble(
      
        ## Each drug has 3 keys: one primary (that always exists) and two optional
        DrugBank_ID = xmlValue(drug["drugbank-id"][[1]]),
        secondary_key = ifelse(length(drug["drugbank-id"]) > 1, xmlValue(drug["drugbank-id"][[2]]), NA),
        third_key = ifelse(length(drug["drugbank-id"]) > 2, xmlValue(drug["drugbank-id"][[3]]), NA),
        
        ## drug attribues
        type = xmlGetAttr(node = drug, name = "type"),
        
        ## drug name
        name = xmlValue(drug[["name"]]),
        
        ## groups
        FDA_label = xmlValue(drug[["fda-label"]])
    )
}
           

#Function to obtain children info 
getChildInfo <- function(base_node,                 
                        child_node,                
                        sub_child_node = NULL,     
                        id = "drugbank-id") {     
    
    
    if (!is.null(base_node[[child_node]])) {
        if (is.null(sub_child_node)) {
            df <- xmlToDataFrame(base_node[[child_node]], 
                                 stringsAsFactors = FALSE)
        } else {
            df <- xmlToDataFrame(base_node[[child_node]][[sub_child_node]], 
                                 stringsAsFactors = FALSE)
        }
    } else {
        df <- NULL
    }
    
    if (!is.null(df) && nrow(df) > 0) {
        parent_key <- NULL
      
        if (!is.null(id)) {
            parent_key <- xmlValue(base_node[id][[1]])
        }
       
        if (!is.null(parent_key)) {
            df$parent_key <- parent_key
        }
    }
    
    return(df)
}

#Function to obtain ATC codes 
get_atc_codes_rec <- function(r, drug_key) {
    tibble(
        atc_code = xmlGetAttr(r, name = "code"),
        
        level_1 = xmlValue(r[[4]]),
        code_1 = xmlGetAttr(r[[4]], name = "code"),
        
        level_2 = xmlValue(r[[3]]),
        code_2 = xmlGetAttr(r[[3]], name = "code"),
        
        level_3 = xmlValue(r[[2]]),
        code_3 = xmlGetAttr(r[[2]], name = "code"),
        
        level_4 = xmlValue(r[[1]]),
        code_4 = xmlGetAttr(r[[1]], name = "code"),
        
        parent_key = drug_key
    )
}


get_atc_codes_df <- function(drug) {
    return (map_df(xmlChildren(drug[["atc-codes"]]),
                   ~ get_atc_codes_rec(.x, xmlValue(drug["drugbank-id"][[1]]))))
}


drugs_IDs <- map_df(all_drugs, ~getDinfo(.x)) 
drugs_IDs <- drugs_IDs %>% dplyr::rename(parent_key = DrugBank_ID)
saveRDS(drugs_IDs, "basic_info.rds")

groups <- map_df(all_drugs, ~getChildInfo(.x, child_node = "groups"))
combine_group <- groups %>% group_by(parent_key) %>% summarise(status = paste0(text,collapse = ", "))
saveRDS(combine_group, "status.rds")

synonyms <-  map_df(all_drugs, ~getChildInfo(.x, child_node = "synonyms"))
combine_synonyms <- synonyms %>% group_by(parent_key) %>% summarise(synonyms = paste0(text,collapse = ", "))

saveRDS(combine_synonyms, "synonyms.rds")
categories <- map_df(all_drugs, ~getChildInfo(.x, child_node = "categories"))
saveRDS(categories, "categories.rds")


atc_codes <- map_df(all_drugs, ~get_atc_codes_df(.x))
saveRDS(atc_codes, "atc_codes.rds")

calculated_properties <- map_df(all_drugs, ~getChildInfo(.x, child_node = "calculated-properties"))

smiles <- calculated_properties %>% filter(kind == "SMILES") %>% select(parent_key, value) %>% dplyr::rename(SMILES = value)
InChI <- calculated_properties %>% filter(kind == "InChI") %>% select(parent_key, value) %>% dplyr::rename(InChI = value)
InChIKey <- calculated_properties %>% filter(kind == "InChIKey") %>% select(parent_key, value) %>% dplyr::rename(InChIKey = value)


#put all data frames into list
df_list <- list(smiles,InChI, InChIKey)

chem_ids <- df_list %>% reduce(full_join, by='parent_key')

#put all data frames into list
df_list2 <- list(drugs_IDs, combine_group, combine_synonyms, atc_codes, chem_ids)

DrugBankDB <- df_list2 %>% reduce(full_join, by='parent_key')

saveRDS(DrugBankDB, "DrugBankDB.rds")

# filter by drug_type
DrugBankDB_mol <- DrugBankDB %>% filter(type == "small molecule")
saveRDS(DrugBankDB_mol, "DrugBankDB_mol.rds")
DrugBankDB_biotec <- DrugBankDB %>% filter(type == "biotech")
saveRDS(DrugBankDB_biotec, "DrugBankDB_biotech.rds")

# Pharmacology
  
indications <- map_df(all_drugs, ~getChildInfo(.x, child_node = "indication"))

  
  

