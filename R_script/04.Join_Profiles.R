# CGaliana
# 03.Join_Profiles
# Version: 1.0.0


#-------------------------------loading packages-------------------------------------------------#

library(stringr)

# Set working directory
setwd("/home/cgaliana/NPC_screening")
currentDataDir = getwd()

#-------------------------------loading metadata files--------------------------------------------# 

load(paste0(currentDataDir,"/Data/genes_id.RData"))


# setting the conditions--------------------------------------------------------------------------# 

clinical_phase <- "rest"

# ------------------------------------------------------------------------------------------------# 

f_list <-list.files(paste0(currentDataDir,"/DEG/",clinical_phase), pattern="*.RData", full.names="TRUE")

#create containers for: 
#save Fold Change (FC) values
df_FC <- data.frame(matrix(ncol = length(f_list), nrow = 978))
rownames(df_FC) <- gene_id
#save average expression (Ave) values
df_Ave <- data.frame(matrix(ncol = length(f_list), nrow = 978))
rownames(df_Ave) <- gene_id

names_col <- vector("list", length(f_list))

#save the number of differentially-expressed genes (DEGs) called significant at an adjusted p-value cutoff (FDR) threshold of 0.05
df_info <- data.frame(matrix(ncol = 8, nrow = length(f_list)))
colnames(df_info) <- c("name_profile", "cmap_name", "clinical_phase", "pert_idose", "pert_itime", "ss_genes", "UP", "DOWN")


for (i in seq_along(f_list)) {
  
  load(f_list[i])
  
  name_profile <- str_match(f_list[i], "DEG/rest/\\s*(.*?)\\s*.RData")[,2]
  comb <- strsplit(name_profile, "_")
  cmap_name <- comb[[1]][1]
  pert_idose <- comb[[1]][2]
  pert_itime <- comb[[1]][3]
  ss_genes <- length(deg[deg$adj.P.Val < 0.05, "adj.P.Val"])
  UP <- length(deg[deg$logFC > 0 & deg$adj.P.Val < 0.05, "logFC"])
  DOWN <- length(deg[deg$logFC < 0 &deg$adj.P.Val < 0.05, "logFC"])
  
  names_col[[i]] <- name_profile
  
  df_info[i, ] <- data.frame(name_profile, cmap_name, clinical_phase, pert_idose, pert_itime, ss_genes, UP, DOWN)
  
  df_FC[, i] <- deg$logFC 
  df_Ave[, i] <- deg$AveExpr
  colnames(df_FC) <- names_col
  colnames(df_Ave) <- names_col
}

write.csv(df_info, file = paste0(currentDataDir, "/resfiles/",clinical_phase,"_ss_genes.csv"))
write.csv(df_FC, file = paste0(currentDataDir, "/resfiles/",clinical_phase,"_FC.csv"))
write.csv(df_Ave, file = paste0(currentDataDir, "/resfiles/",clinical_phase,"_Ave.csv"))

