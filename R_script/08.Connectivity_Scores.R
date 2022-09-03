# CGaliana

################### Connectivity Methods for Drug Repurposing (based on L1000 data) ###################
# Version: 1.0.5                                                                                      #
# Adapted GSEA methods: Weighted Signed Kolmogorov-Smirnov Statistic                                  #
## CMap 1.0 : Connectivity Score (CS) (Lamb et al. 2006)                                              #                                                                                
## CMap 2.0 : Weighted Connectivity Score (WCS) (Subramanian et al.2017)                              #                                                                                                     
## ZangScore : Connection Strength Score (CSS)                                                        #
## Extreme Pairwise Metrics: XSum                                                                     #
## Pairwise metrics: Cor, Spe, Cos                                                                    #
#                                                                                                     #
## CMap 1.0, CSS, XSum implemented in RCSM package.                                                   #
## CMap 2.0 adapted from RCSM package.                                                                #
#######################################################################################################  


suppressPackageStartupMessages({
  library(RCSM)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(doParallel)
  library(lsa)})


# Set working directory

currentDataDir = "/home/cgaliana/AD/"


# settings: gene set size, number of cores, number of permutations 
ncore <- 20
npermut <- 1000
set.seed(1546)
cutoff_FDR <- 0.05

# format database 
# FC <- readRDS("/home/cgaliana/logFC_alllines.rds")
# geneinfo <- read.delim("/home/cristina/Escritorio/TFM/SCRIPTS/GENES/geneinfo_beta.txt", sep = "\t")
# entrez_ID <- as.data.frame(rownames(FC)); names(entrez_ID) <- "gene_id"
# gene_symbol <- geneinfo %>% dplyr::select(gene_id, gene_symbol)
# symbol_FC <- merge(entrez_ID, gene_symbol, by ="gene_id")
# rownames(FC) <- symbol_FC$gene_symbol
#saveRDS(FC, "/home/cgaliana/drugs_DB.rds")


#loading the Drug Profiles DataBase: logFC from analysis routines of DEGs
Drugs_DB <- readRDS("/home/cgaliana/drugs_DB.rds")
#loading the expression signatures of disease
Disease_profiles <- readRDS(paste0(currentDataDir,"AD_profiles.rds"))
# loading metadata for profiles 
meta <- readRDS("/home/cgaliana/meta_alllines.rds") %>% dplyr::rename(pert = name_profile) 


# CMap 2.0 adapted from RCSM package------------------------------------------------------------------------------#

weightedKScore <- function(refMatrix, queryUp, queryDown, permuteNum = npermut, pAdjMethod = "BH", mcCore = ncore) {

  if (is.data.frame(refMatrix)) {refMatrix <- as.matrix(refMatrix)}

  ## Convert the gene expression matrix to ranked list
  matrixToRankedList <- function(refMatrix) {
    refList <- vector("list", ncol(refMatrix))
    for(i in 1:ncol(refMatrix)) {
      ## Sorting genes based on logFC
      refList[[i]] <- refMatrix[order(refMatrix[, i], decreasing=TRUE), i]
    }
    return(refList)
  }
  ## Compute the EScore
  EScore <- function(refList, query) {
    tagIndicator <- sign(match(names(refList), query, nomatch=0))
    noTagIndicator <- 1 - tagIndicator
    N <- length(refList)
    Nh <- length(query)
    Nm <-  N - Nh
    correlVector <- abs(refList)
    sumCorrelTag <- sum(correlVector[tagIndicator == 1])
    normTag <- 1.0/sumCorrelTag
    normNoTag <- 1.0/Nm
    RES <- cumsum(tagIndicator * correlVector * normTag -
                    noTagIndicator * normNoTag)
    maxES <- max(RES)
    minES <- min(RES)
    maxES <- ifelse(is.na(maxES), 0, maxES)
    minES <- ifelse(is.na(minES), 0, minES)
    ifelse(maxES > - minES, maxES, minES)
  }
  # Computing WTCS as in Subramanian et al. 2017
  weight1 <- function(refList, queryUp, queryDown) {
    scoreUp <- EScore(refList, queryUp)
    scoreDown <- EScore(refList, queryDown)
    ifelse(scoreUp * scoreDown <= 0, (scoreUp - scoreDown)/2, 0)
  }

  ## Prepare the ranked reference lists
  refList <- matrixToRankedList(refMatrix)
  ## Prepare the up and down signatures
  queryUp <- intersect(queryUp, rownames(refMatrix))
  queryDown <- intersect(queryDown, rownames(refMatrix))
  WCSscore <- mclapply(refList, weight1, queryUp = queryUp,
                    queryDown = queryDown, mc.cores = mcCore)
  WCSscore <- as.vector(do.call(rbind, WCSscore))
  
  #null distribution of ES
  permuteScore <- matrix(0, ncol = permuteNum, nrow = ncol(refMatrix))
  for(n in 1:permuteNum) {
    bootUp <- sample(rownames(refMatrix), size = length(queryUp))
    bootDown <- sample(rownames(refMatrix), size = length(queryDown))
    bootScore <- mclapply(refList, weight1, queryUp = bootUp,
                          queryDown = bootDown, mc.cores = mcCore)
    permuteScore[, n] <- as.vector(do.call(rbind, bootScore))
  }
  permuteScore[is.na(permuteScore)] <- 0
  ## Compute the p-value based on bootstrap method
  pValue <- rowSums(abs(permuteScore) >= abs(WCSscore)) / permuteNum
  pAdjust <- p.adjust(pValue, method = pAdjMethod)

  scoreResult <- data.frame(Score = WCSscore, pValue = pValue, pAdjValue = pAdjust)
  rownames(scoreResult) <- colnames(refMatrix)
  
  scoreResult  <- scoreResult %>%
                     arrange(Score) %>%
                     rownames_to_column() %>%
                     rowwise() %>% 
                     mutate(pert = strsplit(rowname, split="__")[[1]][1]) %>%
                     inner_join(meta, by = "pert") %>% 
                     dplyr::rename(name_profile = pert) %>%
                     select(name_profile, Score, pValue, pAdjValue, cmap_name, pert_idose, pert_itime, cell_iname)
  # compute NCS
  cell_type <- paste(scoreResult$cell_iname, as.character(ifelse(scoreResult$Score > 0, "UP", "DOWN")), sep = "_")
  ES.values <-as.numeric(scoreResult$Score)
  ES.values[ES.values == 0] <- NA
  mean_by_celltype <- tapply(ES.values, cell_type, mean, na.rm=TRUE)
  ncs <- as.numeric(scoreResult$Score) / abs(mean_by_celltype[cell_type]) 
  scoreResult <- cbind(scoreResult, ncs) %>% filter(Score < 0 & pAdjValue < cutoff_FDR) 
  return(scoreResult)
}


 #-----------------------------------------------------------------------------------------------------------------------------------------#
#Cor
Corr <-function(Disease_profiles, drugs_db, corrmethod, ncore){
   
    query0 <- Disease_profiles[, c("gene_symbol", "logFC")] %>% 
                   filter(gene_symbol %in% shared) %>% 
                   column_to_rownames("gene_symbol") %>% 
                   select(logFC)
    
    query0_subset <- as.numeric(query0[c(shared),])
    drug_profile <- as.matrix(drugs_db[c(shared),])
  
 
  XCor <-function(drug_profile, query_subset, corrmethod) {
	stopifnot(is.vector(query_subset), dim(drug_profile)[1] == length(query_subset))
	cor(drug_profile, query_subset, use = "complete.obs", method = corrmethod)
}

cl <- makeCluster(ncore) 
registerDoParallel(cl)

cor_values <- foreach (i = 1:ncol(drug_profile)) %dopar% {
  XCor(drug_profile[,i], query0_subset,corrmethod)
}

stopCluster(cl)
score <- as.vector(do.call(rbind, cor_values))

# container
permut.values <- matrix(0, ncol = npermut, nrow = ncol(drug_profile))
  for(n in 1:npermut) {
    ## Prepare the random query signatures
    bootquery.names <- sample(rownames(drug_profile), size = length(query0_subset ))
    bootquery.values <- list(drug_profile[bootquery.names,])
    ## Compute the random scores 
    boot.scores <- mclapply(bootquery.values, XCor, query = query0_subset, corrmethod = corrmethod, mc.cores = ncore)
    permut.values[, n] <- as.vector(do.call(rbind, boot.scores))
  }
  permut.values[is.na(permut.values)] <- 0
  
  ## Compute the p-values
  pValue <- rowSums(abs( permut.values) >= abs(score)) / npermut
  ## Compute the adjusted p-values
  pAdjust <- p.adjust(pValue, method = "BH")

  scoreResult <- data.frame(pert= colnames(drugs_db), Score = score, pValue = pValue, pAdjValue = pAdjust) %>% 
    filter(pAdjValue < cutoff_FDR & Score < 0) %>% 
    arrange(Score) %>% 
    rowwise() %>% 
    mutate(pert = strsplit(pert, split="__")[[1]][1]) %>% 
    inner_join(meta, by = "pert") %>% 
    dplyr::rename(name_profile = pert) %>%
    select(name_profile, Score, pValue, pAdjValue, cmap_name, pert_idose, pert_itime, cell_iname)
  
}
#---------------------------------------------------------------------------------------------------------------------------------------#
# Cos
Cosine <-function(Disease_profiles, drugs_db, ncore){
   
  # Building the query: get the logFC values from UP & DOWN queries
  query0 <- Disease_profiles[, c("gene_symbol", "logFC")] %>% 
                    filter(gene_symbol %in% shared) %>% 
                    column_to_rownames("gene_symbol") %>% 
                    select(logFC)
    
  query0_subset <- as.numeric(query0[c(shared),])
  drug_profile <- as.matrix(drugs_db[c(shared),])

 
  XCosi <-function(drug_profile, query_subset) {
	stopifnot(is.vector(query_subset), dim(drug_profile)[1] == length(query_subset))
	lsa::cosine(as.numeric(drug_profile), query_subset)
}

  cl <- makeCluster(ncore) 
registerDoParallel(cl)

cor_values <- foreach (j = 1:ncol(drug_profile)) %dopar% {
  XCosi(drug_profile[,j], query0_subset)
}

stopCluster(cl)
score <- as.vector(do.call(rbind, cor_values))

# container
permut.values <- matrix(0, ncol = npermut, nrow = ncol(drug_profile))
  for(n in 1:npermut) {
    ## Prepare the random query signatures
    bootquery.names <- sample(rownames(drug_profile), size = length(query0_subset))
    bootquery.values <- drug_profile[bootquery.names,]
    bootquery.values <- split(bootquery.values, rep(1:ncol(bootquery.values), each = nrow(bootquery.values)))
    ## Compute the random scores 
    boot.scores <- mclapply(bootquery.values, XCosi, query = query0_subset, mc.cores = ncore)
    permut.values[, n] <- as.vector(do.call(rbind, boot.scores))
  }
  permut.values[is.na(permut.values)] <- 0
  
  ## Compute the p-values
  pValue <- rowSums(abs( permut.values) >= abs(score)) / npermut
  ## Compute the adjusted p-values
  pAdjust <- p.adjust(pValue, method = "BH")

  scoreResult <- data.frame(pert= colnames(drugs_db), Score = score, pValue = pValue, pAdjValue = pAdjust) %>% 
    filter(pAdjValue < cutoff_FDR & Score < 0) %>% 
    arrange(Score) %>% 
    rowwise() %>% 
    mutate(pert = strsplit(pert, split="__")[[1]][1]) %>% 
    inner_join(meta, by = "pert") %>% 
    dplyr::rename(name_profile = pert) %>%
    select(name_profile, Score, pValue, pAdjValue, cmap_name, pert_idose, pert_itime, cell_iname)
  
} 
#---------------------------------------------------------------------------------------------------------------------#
t <- proc.time() 

for (i in 1:length(Disease_profiles)) {
  p_name <- names(Disease_profiles[i])
  cat("\nAnalyzing:", p_name) 
  
  # getting the shared genes 
  genes_drugs <- Drugs_DB %>% rownames_to_column("gene_symbol") %>% .$gene_symbol 
  genes_EM <- Disease_profiles[[i]]$gene_symbol
  shared <- intersect(genes_drugs, genes_EM)
  cat("\nshared genes: ", length(shared))
  
  Drugs_DB <- Drugs_DB[shared,]
  
  # Building the query input 
  ## keep the logFC and gene symbol and filter by shared genes 
  query0 <- Disease_profiles[[i]][, c("gene_symbol", "logFC","p.adjust.fdr")] %>% filter(gene_symbol %in% shared) 
  query_UP <- query0 %>% filter(p.adjust.fdr < 0.05 & logFC > 0) %>% arrange(desc(logFC)) %>% .$logFC
  names(query_UP) <- query0 %>% filter(p.adjust.fdr < 0.05 & logFC > 0) %>% arrange(desc(logFC)) %>% .$gene_symbol
  UP <- names(query_UP)
  
  query_DOWN <- query0 %>% filter(p.adjust.fdr < 0.05 & logFC < 0) %>% arrange(desc(logFC)) %>% .$logFC
  names(query_DOWN) <- query0 %>% filter(p.adjust.fdr < 0.05 & logFC < 0) %>% arrange(desc(logFC)) %>% .$gene_symbol
  DOWN <- names(query_DOWN)
  
  cat("\ngene set size:\n-UP:", length(UP), "\n-DOWN:", length(DOWN))
  
  # CMap 1.0-RCSM
  cat("\nCMap 1.0 method:\n")
  cmap1 <- KSScore(Drugs_DB, UP, DOWN, permuteNum = npermut, pAdjMethod = "BH", mcCore = ncore)
  ## Computing scaled score [1,-1] dividing positive and negative scores with the maximum positive score and the absolute value of the minimum        negative score, respectively. 
  max.pos <- max(cmap1 %>% filter(Score > 0) %>% .$Score)
  max.neg <- max(abs(cmap1 %>% filter(Score < 0) %>% .$Score))
  cmap1 <- cmap1 %>% mutate(scaled_score = ifelse(Score > 0, Score/max.pos, Score/max.neg)) %>% 
                     filter(Score < 0 & pAdjValue < cutoff_FDR) %>%
                     arrange(Score) %>%
                     rownames_to_column() %>%
                     rowwise() %>% 
                     mutate(pert = strsplit(rowname, split="__")[[1]][1]) %>%
                     inner_join(meta, by = "pert") %>% 
                     dplyr::rename(name_profile = pert) %>%
                     select(name_profile, Score, scaled_score, pValue, pAdjValue, cmap_name, pert_idose, pert_itime, cell_iname) 
    
  saveRDS(cmap1,paste0(currentDataDir, p_name,"_cmap1.rds" ))        
  write.csv(cmap1, paste0(currentDataDir, p_name,"_cmap1.csv"))
  
  
  
  # CMap 2.0-Adapted from RCSM
  cat("\nCMap 2.0 method:\n")
  cmap2 <- weightedKScore(Drugs_DB, UP, DOWN, permuteNum = npermut, pAdjMethod = "BH", mcCore = ncore)
    
   saveRDS(cmap2,paste0(currentDataDir, p_name, "_cmap2.rds" ))   
   write.csv(cmap2, paste0(currentDataDir, p_name, "_cmap2.csv"))
  
  
  # ZangScore- RCSM
  cat("\nZangScore method:\n")
  res_ZangScore <- ZhangScore(Drugs_DB, UP, DOWN, permuteNum = npermut, pAdjMethod = "BH", mcCore = ncore)
  ZangS <- res_ZangScore %>% 
           filter(pAdjValue < cutoff_FDR & Score < 0) %>% arrange(Score) %>% 
           rownames_to_column() %>% 
           rowwise() %>% 
           mutate(pert = strsplit(rowname, split="__")[[1]][1]) %>% 
           inner_join(meta, by = "pert") %>% 
           dplyr::rename(name_profile = pert) %>%
           select(name_profile, Score, pValue, pAdjValue, cmap_name, pert_idose, pert_itime, cell_iname)
  
  saveRDS(ZangS,paste0(currentDataDir, p_name, "_ZangS.rds" )) 
  write.csv(ZangS, paste0(currentDataDir, p_name, "_ZangS.csv"))
  
  #XSum 
  cat("\nXSum method-size 200:\n")
  res_Xsum <- XSumScore(Drugs_DB, UP, DOWN, topN = 200, permuteNum = npermut, pAdjMethod = "BH", mcCore = ncore)
  Xsum200 <- res_Xsum %>% 
           filter(pAdjValue < cutoff_FDR & Score < 0) %>% arrange(Score) %>% 
           rownames_to_column() %>% 
           rowwise() %>% 
           mutate(pert = strsplit(rowname, split="__")[[1]][1]) %>% 
           inner_join(meta, by = "pert") %>% 
           dplyr::rename(name_profile = pert) %>%
           select(name_profile, Score, pValue, pAdjValue, cmap_name, pert_idose, pert_itime, cell_iname)
  
  write.csv(Xsum200, paste0(currentDataDir, p_name, "_200_XsumS.csv"))
  saveRDS(Xsum200,paste0(currentDataDir, p_name, "_200_XsumS.rds" )) 
  
  
  #XSum 
  cat("\nXSum method-size 300:\n")
  res_Xsum <- XSumScore(Drugs_DB, UP, DOWN, topN = 300, permuteNum = npermut, pAdjMethod = "BH", mcCore = ncore)
  Xsum300 <- res_Xsum %>% 
           filter(pAdjValue < cutoff_FDR & Score < 0) %>% arrange(Score) %>% 
           rownames_to_column() %>% 
           rowwise() %>% 
           mutate(pert = strsplit(rowname, split="__")[[1]][1]) %>% 
           inner_join(meta, by = "pert") %>% 
           dplyr::rename(name_profile = pert) %>%
           select(name_profile, Score, pValue, pAdjValue, cmap_name, pert_idose, pert_itime, cell_iname)
  
  write.csv(Xsum300, paste0(currentDataDir, p_name, "_300_XsumS.csv"))
  saveRDS(Xsum300,paste0(currentDataDir, p_name, "_300_XsumS.rds" ))


  #XSum 
  cat("\nXSum method-size 400:\n")
  res_Xsum <- XSumScore(Drugs_DB, UP, DOWN, topN = 400, permuteNum = npermut, pAdjMethod = "BH", mcCore = ncore)
  Xsum400 <- res_Xsum %>% 
           filter(pAdjValue < cutoff_FDR & Score < 0) %>% arrange(Score) %>% 
           rownames_to_column() %>% 
           rowwise() %>% 
           mutate(pert = strsplit(rowname, split="__")[[1]][1]) %>% 
           inner_join(meta, by = "pert") %>% 
           dplyr::rename(name_profile = pert) %>%
           select(name_profile, Score, pValue, pAdjValue, cmap_name, pert_idose, pert_itime, cell_iname)
  
  write.csv(Xsum400, paste0(currentDataDir, p_name, "_400_XsumS.csv"))
  saveRDS(Xsum400,paste0(currentDataDir, p_name, "_400_XsumS.rds" ))


  #Cor method
  cat("\nCor method:\n")
  Cor <- Corr(Disease_profiles[[i]], Drugs_DB, "pearson", ncore)
  write.csv(Cor, paste0(currentDataDir, p_name, "_Cor.csv"))
  saveRDS(Cor,paste0(currentDataDir, p_name, "_Cor.rds" )) 

  #Spe method
  cat("\nSpe method:\n")
  Spe <- Corr(Disease_profiles[[i]], Drugs_DB, "spearman", ncore)
  write.csv(Spe, paste0(currentDataDir, p_name,  "_Spe.csv"))
  saveRDS(Spe,paste0(currentDataDir, p_name, "_Spe.rds" )) 
   
  #Cos method
  cat("\nCos method:\n")
  Coss <- Cosine(Disease_profiles[[i]], Drugs_DB, ncore)
  write.csv(Coss, paste0(currentDataDir, p_name, "_Coss.csv"))
  saveRDS(Coss,paste0(currentDataDir, p_name,  "_Coss.rds" )) 


}

proc.time()-t  






