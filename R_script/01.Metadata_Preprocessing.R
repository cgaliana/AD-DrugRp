# CGaliana
# Metadata Preprocessing
# Version: 1.0.1


#----------------------------install.packages---------------------------------------------------#
#install.packages("cmapR", dependencies = TRUE)

#-------------------------------loading packages-------------------------------------------------# 

library(cmapR)
library (dplyr, warn.conflicts = FALSE)


# Set working directory
currentDataDir = getwd()

#-------------------------------loading metadata files------------------------------------------- # 

# instinfo_beta.txt <- LEVEL 3 LINCS metadata  (dataset CMap_LINCS2020)
# Repurposing_drug_Hub file 
# gene_info
# path to GCTx files

LINCS_meta <- read.delim(paste0(currentDataDir,"/MetaFiles/instinfo_beta.txt"), header = T, sep = "\t", na.strings = c("","NA"))

rephub <- read.delim(paste0(currentDataDir,"/MetaFiles/repurposing_drugs_20200324.txt"), header = T, sep = "\t", na.strings = c("","NA"))

gene_info <- read.delim(paste0(currentDataDir,"/MetaFiles/geneinfo_beta.txt"), header = T, sep = "\t", na.strings = c("","NA")) %>%
  filter(feature_space == "landmark") %>%
  select(everything())

GCTx_cp_path <- "/home/cgaliana/lincs_trt_cp/level3_beta_trt_cp_n1805898x12328.gctx"
GCTx_ct_path <- "/home/cgaliana/lincs_ctl/level3_beta_ctl_n188708x12328.gctx"


#-------------------------------standard workflow functions--------------------------------------- #

# getDrugs()
## Description:  Drug selection function based on clinical phase.
##               Repurposing Drug Hub (RDH) as reference annotation database. 
##               Intersection by cmap_name (curated name of drug by the CMap team).   
## Parameters:
###   development_phase (atomic vector of characters) options: (RDH designation)
###                                  “Launched”, “Phase 1”,  “Phase 2”,  “Phase 3”, “Phase 1/Phase 2”,  “Phase 2/Phase 3”,  “Preclinical”,  “Withdrawn”.  
###   database (data.frame)  from the RDH 
## Values:
###   names (atomic vector of characters) drugs names


getDrugs <- function(development_phase, database) {
  
  cphase <- database %>% 
    distinct(clinical_phase)%>%
    .$clinical_phase
  
  drugs <- database %>% 
    dplyr::rename(cmap_name = pert_iname) %>%
    inner_join(LINCS_meta, by = "cmap_name") %>%
    select(cmap_name, clinical_phase) 
  
  names <- drugs %>% filter(clinical_phase %in% development_phase) %>% distinct(cmap_name) %>% .$cmap_name
   if (length(names) == 0) {
     stop("Unknown input!")
  }
  return(names)
}

# getDrugsout()
## Description:  Drug selection function for unannotated drugs in the RDH.   
## Parameters:  
###   database (data.frame)  from the RDH 
## Values:
###   names (atomic vector of characters) drugs names


getDrugsout <- function(database) {

  database <- database %>% dplyr::rename(cmap_name = pert_iname) 
  
  names<- LINCS_meta  %>%
    full_join(database, by =  "cmap_name") %>% 
    filter(is.na(clinical_phase)) %>%
    distinct(cmap_name) %>% .$cmap_name

return(names)
  }


# getGeneids()
## Description ::: Function to obtain the specific feature_space gene_ids in Lexicographical Order
## Parameters
###   feature: (character) (default: "landmark") 
###   gene_metadata (data.frame)  Gene annotation
## Values:
###   gene_ids (character)

getGeneids <- function(feature = "landmark", gene_metadata) {
  #get the Entrez_IDs (rows of the GCTx file)  
  rid <- read_gctx_meta(GCTx_cp_path, dim="row")
  
  #get the genes for the selected feature_space. 
  gene_ids <- merge(rid, gene_metadata, 
                    by.x = "id", by.y = "gene_id", 
                    all.x = TRUE, all.y = TRUE, 
                    sort = TRUE,  
                    incomparables = "NA", 
  ) %>%
    filter(feature_space == feature)  %>%
    .$id
  
  return(gene_ids)
}

# cleanmeta()
## Description: filter function by: 1st) cell line, 2nd) drugs, 3rd) dose and time, and 4th) technical replicates >=3.
## Parameters
###   cell (character) cell line
###   drugs (character) drugs previously selected with getDrugs() or getDrugsout() functions. 
###   dose (character)  doses 
###   time (character)  times 
## Values:
###   data (data.frame) shortest metadata file with all instance/profiles (cell line-drug-dose-time combination). 
###                     it include new covariate (Batch_id)  

cleanmeta <- function(cell, drugs, dose, time) {
  data <-LINCS_meta %>%
    
    filter(cell_iname == cell & pert_type == "trt_cp") %>%
    filter(cmap_name %in% dg) %>%
    filter(pert_idose %in% ds) %>%
    filter(pert_itime %in% t) %>%
    add_count(cmap_name, pert_idose, pert_itime, wt = NULL, name = "technical_replicates", sort = TRUE) %>%
    filter(technical_replicates >= 3) %>%
    rowwise() %>%
    mutate(batch_id = paste(strsplit(rna_plate, "_")[[1]][1:3], collapse="_")) %>%
    select(sample_id, pert_type, cmap_name, pert_id, pert_idose, pert_itime, project_code, cell_mfc_name, batch_id, rna_plate, det_plate, rna_well, det_well)
  
  return(data)
}

# getcomb()
## Description: To get the all instances into a object (without metadata information)
## Parameters
###   filtered_metadata: shortest metadata file with all instance/profiles (cell line-drug-dose-time combination).
## Values:
###   data (data.frame) 

getcomb <- function(filtered_metadata) {
  data <-filtered_metadata %>%
    select(cmap_name,pert_idose, pert_itime) %>%
    distinct() %>%
    
    return(data)
}

# setting the Conditions----------------------------------------------------------------------------------------·# 

cell <- "NPC"
#drug_phase <- c("Phase 1",  "Phase 2",  "Phase 3", "Phase 1/Phase 2",  "Phase 2/Phase 3",  "Preclinical",  "Withdrawn")
# only for the first cell line analysis, after that load the .RData
ds <- c("0.04 uM", "0.12 uM", "0.37 uM", "1.11 uM", "3.33 uM", "10 uM","30 uM", "60 uM", "125 uM")
t <- c("6 h", "24 h") 


# main ------------------------------------------------------------------------------------------------------------#

#dg <- getDrugs(drug_phase, rephub)
dg <- getDrugsout(rephub)
gene_id <- getGeneids("landmark", gene_info)

save(dg, file = paste0(currentDataDir,"/Data/drugnames.RData"))
save(gene_id, file = paste0(currentDataDir,"/Data/genes_id.RData"))
save(ds, file = paste0(currentDataDir,"/Data/doses.RData"))
save(t, file = paste0(currentDataDir,"/Data/times.RD"))
f_meta <- cleanmeta(cell, dg, ds, t)
save(f_meta, file = paste0(currentDataDir,"/Data/meta",cell,".RData"))


profiles <- getcomb(f_meta)
save(profiles, file = paste0(currentDataDir,"/Data/profiles",cell,".RData"))

cat(paste("drugs screened in", cell, ":\n"), file=stdout())
cat(paste("Total drugs: ", length(profiles %>% distinct(cmap_name)%>% .$cmap_name), "\n"), file=stdout()) 
cat(paste("Total profiles: ", nrow(profiles), "\n"), file=stdout()) 

