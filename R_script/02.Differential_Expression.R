# CGaliana
# Differential Expression 
# Version: 1.0.5


## Standard pipeline for each drug, dose and time combination: 
###       1. Filter treatment metadata (cell_linemeta.RData) for a specified profile. Dataframe for the SummarizedExperiment object. Function: getmetacomb() 
###       2. Save the unique identifiers of the treatment (sample_ids).
###       3. Identify the corresponding control profiles (pert_type = "ctl_vehicle", DMSO) coming from the same plates as the treatment.  
###       4. Save the metadata for the controls. Dataframe for the SummarizedExperiment object. Function: cleanmeta_ctl()
###       5. Save the unique identifiers of the controls (sample_ids).
###       6. Build a GCT object. Function: buildGCTx
###       7. Assemble the expression data and metadata into a SummarizedExperiment object and return it.
###       8. Perform differential expression.

#-------------------------------loading packages-------------------------------------------------# 

suppressPackageStartupMessages({
library(cmapR)
library (dplyr)
library(SummarizedExperiment)
library(limma)})


# Set working directory
setwd("/home/cgaliana/NPC_screening")
currentDataDir = getwd()

#-------------------------------loading metadata files--------------------------------------------# 

# instinfo_beta.txt <- LEVEL 3 LINCS metadata  (dataset CMap_LINCS2020)
LINCS_meta <- read.delim(paste0(currentDataDir,"/MetaFiles/instinfo_beta.txt"), header = T, sep = "\t", na.strings = c("","NA"))

# path to GCTx files
GCTx_cp_path <- "/home/cgaliana/lincs_trt_cp/level3_beta_trt_cp_n1805898x12328.gctx"
GCTx_ct_path <- "/home/cgaliana/lincs_ctl/level3_beta_ctl_n188708x12328.gctx"

#-------------------------------standard workflow functions---------------------------------------#


getmetacomb <-function(drug, dose, time){
  ids <- f_meta %>%
    filter(cmap_name == drug, pert_idose == dose, pert_itime == time)
  }


cleanmeta_ctl <- function(cell, rna_plates) {
  
  data <-LINCS_meta %>%
    filter(cell_iname == cell & pert_type == "ctl_vehicle") %>%
    filter(rna_plate %in% rna_plates) %>%
    add_count(cmap_name, pert_idose, pert_itime, wt = NULL, name = "technical_replicates", sort = TRUE) %>%
    filter(technical_replicates >= 3) %>%
    rowwise() %>%
    mutate(batch_id = paste(strsplit(rna_plate, "_")[[1]][1:3], collapse="_")) %>%
    select(sample_id, pert_type, cmap_name, pert_id, pert_idose, pert_itime, project_code, cell_mfc_name, batch_id, rna_plate, det_plate, rna_well, det_well)
  
  return(data)
  }


buildGCTx <- function(cp_path, ctl_path, id_genes, ids_cp, ids_ctl, metadata_compounds, metadata_controls){
  # build separately the corresponding GCT objects
  GCT_cp <- parse_gctx(cp_path, rid = id_genes, cid = ids_cp, matrix_only = FALSE)
  GCT_clt <- parse_gctx(ctl_path, rid = id_genes, cid = ids_ctl, matrix_only = FALSE)
  
  #Merge the GCT objects
  GCT_merge <- merge_gct(GCT_cp, GCT_clt, dim = "column", matrix_only = FALSE)
  
  # jointly metadata 
  metadata <- rbind(metadata_compounds, metadata_controls)
  GCT_merge <- annotate_gct(GCT_merge, metadata, dim = "column", keyfield = "sample_id")
  
  return(GCT_merge)
  }


diffexprs <- function(se){
  # defining factor groups
  treatment <- as.factor(colData(se)[, "pert_type"])
  batch <- as.factor(colData(se)[, "batch_id"])
  
  form <- "~ treatment";
  if(length(levels(batch))>1) 
    form <- paste(form,"batch",sep="+")
  
  cat("\ndesign model\n",form)
  
  design <- model.matrix(as.formula(form))
  
  # fit the linear models-gene-wide limma test
  fit <- lmFit(assay(se), design)
  # adjust the hierarchical model
  fit2 <- eBayes(fit)
  DEG_limma <- topTable(fit2, coef = 2, number = Inf, adjust = "BH", sort.by = "none", genelist = gene_id)
  
  return(DEG_limma)
  
}

# setting the conditions----------------------------------------------------------------------------------------·#  

cell <- "NPC"
load(paste0(currentDataDir,"/Data/profiles",cell,".RData")) 
load(paste0(currentDataDir,"/Data/meta",cell,".RData"))
load(paste0(currentDataDir,"/Data/genes_id.RData"))


t <- proc.time() 

#df_container
info_samples <- data.frame(matrix(ncol =8, nrow = nrow(profiles)))
colnames(info_samples) <- c("profile", "drug", "dose", "time", "nsamples_trt", "nsamples_ctl", "sample_size", "n_batch")

for (row in 1:nrow(profiles)) {     
  drug <- profiles[row, "cmap_name"]
  dose  <- profiles[row, "pert_idose"]
  time <-  profiles[row, "pert_itime"]
  profile <- gsub(" ", "", paste(drug, dose, time, sep = "_")) 
  
  cat(paste("\n Analyzing profile: ",profile,"\n"), file=stdout()) 
  
  # get metadata df of treatments  
  meta_cp <- getmetacomb(drug, dose, time)
  
  # get treatment sample_ids
  ids_cp <- meta_cp[["sample_id"]]
  
  # get rna_plate IDS
  rna_plates <- meta_cp[["rna_plate"]]
  
  # metada df of controls ERROR HANDLING
  meta_ctl <- tryCatch(
      cleanmeta_ctl(cell, rna_plates),
      error=function(e) e
      )
      if(!inherits(meta_ctl, "error")){
      #PASS

      # save number of batchs	
       n_batch <- nrow(meta_ctl %>% select(batch_id) %>% distinct())
  
      # get controls sample_ids
      ids_ctl <- meta_ctl[["sample_id"]]
  
      # Build the GCTx object from treatments & Controls
      gct <- buildGCTx(GCTx_cp_path, GCTx_ct_path, gene_id, ids_cp, ids_ctl, meta_cp, meta_ctl)
  
      # Convert GCT to se
      se <- as(gct, "SummarizedExperiment")
      saveRDS(se,paste0(currentDataDir,"/se/",desig_class, "/",profile,".rds"))	
      # Container for save info
      nsamples_trt <- length(colData(se)[colData(se)$pert_type == "trt_cp", "pert_type"])
      nsamples_ctl <- length(colData(se)[colData(se)$pert_type == "ctl_vehicle", "pert_type" ])
      sample_size <- length(colData(se)[, "pert_type"])
    
      info_samples[row, ] <- data.frame(profile, drug, dose, time, nsamples_trt, nsamples_ctl, sample_size, n_batch)

      #identify Differentially Expressed Genes (DEG)
      deg <- diffexprs(se)
      save(deg, file = paste0(currentDataDir, paste0("/DEG/",profile,".RData")))

      cat(paste("\nsuccessfully completed! ",profile,"\n"), file=stdout())
      } 
}

write.csv(info_samples, file = paste0(currentDataDir, "/resfiles/info_samples.csv"))

proc.time()-t  

