
## package
{
  library(data.table)
  library(dplyr)
  library(stringr)
  library(openxlsx)
  library(clusterProfiler)
  library(tidyr)
  library(org.Hs.eg.db)
  library(survival)
  library(survminer)
  library(phenoTest)
  library(GSVA)
  library(Seurat)
  library(harmony)
  library(DoubletFinder)
  library(CellChat)
  library(Startrac)
  library(AUCell)
  library(BayesPrism)
  library(DESeq2)
  library(alakazam)
  library(scRepertoire)
  library(monocle)
  library(CytoTRACE) 
  library(ggplot2)
  library(ggthemes)
  library(scales)
  library(viridis)
  library(patchwork)
  library(ComplexHeatmap)
  library(ggsignif)
  library(ggalluvial)
  library(ggpie)
  library(forestploter)
}

## Functions
fucs <- list.files(function_path,full.names=T)
lapply(fucs,function(x) {source(x,encoding="utf-8");return("Yes")})

## local fucntion
{
  groupMeans <- function(mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
      if (sparse) {
        Matrix::rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
      }
      else {
        rowMeans(mat[, which(groups == x), drop = F], na.rm = na.rm)
      }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
  }
}

## color
{
  # sample
  {
    color_sample <- c(
      "MA1_pre"="#1f77b4",
      "MA1_post"="#aec7e8",
      "MA2_pre"="#ff7f0e",
      "MA2_post"="#ffbb78",
      "MA3_pre"="#2ca02c",
      "MA3_post"="#98df8a",
      "MA4_pre"="#e377c2",
      "MA4_post"="#f7b6d2",
      "MA5_pre"="#9467bd",
      "MA5_post"="#c5b0d5",
      "MA6_pre"="#8c564b",
      "MA6_post"="#c49c94",
      "MA7_pre"="#d62728",
      "MA7_post"="#ff9896"
    )
    color_patient <- c(
      "MA1"="#1f77b4",
      "MA2"="#ff7f0e",
      "MA3"="#2ca02c",
      "MA4"="#e377c2",
      "MA5"="#9467bd",
      "MA6"="#8c564b",
      "MA7"="#d62728"
    )
  }
  
  # label
  {
    color_combine <- c('pre_R'='#71afe5',
                       'post_R'='#106ebe',
                       'pre_NR'='#ffa2a4',
                       'post_NR'='#ff0700')
    color_combine_fill <- c('pre_R'='#a3c6e6',
                            'post_R'='#88ABCB',
                            'pre_NR'='#ffd4d2',
                            'post_NR'='#ffa58e')
    color_ISG <- c('ISG_high'='#F9C556',
                   'ISG_low'='#376FB3')
    color_treat <- c('Pre'='#7AB774',
                     'Post'='#B0BBDE')
    color_response_fill <- c('R'='#a3c6e6','NR'='#ffd3d3')
    color_response <- c('R'='#88ABCB','NR'='#ffa3a5') 
  }
  
  # all cells
  {
    color_all_celltype1 <- c('CD4_Naive'='#c8e6c9',
                             'CD4_Mem'='#9EBEE2',
                             'CD4_GZMK'= '#BCBD22',
                             'CD4_EFF'='#ffa58d',
                             'CD4_Tfh'='#9E76C1',
                             'CD4_Treg'='#FF9332',
                             'CD8_Naive'="#0a7b37",
                             'CD8_Tcm'="#3BBEB2",
                             'CD8_GZMK'="#5D9FCF",
                             'CD8_CTL'='#cc101c',
                             'GDT'='#C8969A',
                             'MAIT'='#FB9697',
                             'NKT'='#5976ba',
                             'Prolif.T'='#547689',
                             'NK_CD16'="#1F78B4",
                             'NK_CD56'="#98D083",
                             "Mono_CD14"='#e67f26',
                             "Mono_FCGR3A"='#87b7d8',
                             "Macro_C1QC"='#e078a8',
                             "Macro_NFKB1"='#e1ab8c',
                             "Macro_CXCL10"='#f7cc76',
                             "Macro_LGALS3"='#86d4d8',
                             "Macro_MT1H"='#7bed9f',
                             "Neutrophils"='#8b7042',
                             "cDC1"='#f7ce9c',
                             "cDC2"='#3478b3',
                             'cDC3_LAMP3'='#b83570',
                             "pDC"='#a7cfe7',
                             'Naive B'="#1f77b4",
                             'Memory B'="#aec7e8",
                             'Plasma'="#2ca02c",
                             'Plasmablast'="#98df8a",
                             "Epithelial"='#9668BC'
    )
    color_umap1 <- c('CD4_Naive'='#c8e6c9',
                     'CD4_Mem'='#9EBEE2',
                     'CD4_GZMK'= '#BCBD22',
                     'CD4_EFF'='#ffa58d',
                     'CD4_Tfh'='#9E76C1',
                     'CD4_Treg'='#FF9332',
                     'CD8_Naive'="#0a7b37",
                     'CD8_Tcm'="#3BBEB2",
                     'CD8_GZMK'="#5D9FCF",
                     'CD8_CTL'='#cc101c',
                     'NKT'='#96CC86',
                     'GDT'='#C8969A',
                     'MAIT'='#FB9697',
                     'Prolif.T'='#547689',
                     'NK'='#B1C1D7',
                     "Monocyte"='#e67f26',
                     "Macrophage"='#e1ab8c',
                     "Neutrophil"='#8b7042',
                     "DC"='#6F9C5A',
                     "pDC"='#a7cfe7',
                     "Plasma"='#9d9d82',
                     "B cells"='#d5c8a0',
                     "Epithelial"='#9668BC')
  }
  
  #T
  { 
    color_T <- c('CD4_Naive'='#c8e6c9',
                 'CD4_Mem'='#9EBEE2',
                 'CD4_GZMK'= '#BCBD22',
                 'CD4_EFF'='#ffa58d',
                 'CD4_Tfh'='#9E76C1',
                 'CD4_Treg'='#FF9332',
                 'CD8_Naive'="#0a7b37",
                 'CD8_Tcm'="#3BBEB2",
                 'CD8_GZMK'="#5D9FCF",
                 'CD8_CTL'='#cc101c',
                 'NKT'='#96CC86',
                 'GDT'='#C8969A',
                 'MAIT'='#FB9697',
                 'Prolif.T'='#547689'
    )
  }
  
  # Myeloid
  {
    color_Myeloid <- c("Mono_CD14"='#e67f26',
                       "Mono_FCGR3A"='#87b7d8',
                       "Macro_C1QC"='#e078a8',
                       "Macro_NFKB1"='#e1ab8c',
                       "Macro_CXCL10"='#f7cc76',
                       "Macro_LGALS3"='#86d4d8',
                       "Macro_MT1H"='#9467BD',
                       "Neutrophils"='#8b7042',
                       "cDC1"='#f7ce9c',
                       "cDC2"='#3478b3',
                       'cDC3_LAMP3'='#b83570',
                       "pDC"='#a7cfe7')
    color_macro <- c(
      "Macro_C1QC"='#e078a8',
      "Macro_NFKB1"='#e1ab8c',
      "Macro_CXCL10"='#f7cc76',
      "Macro_LGALS3"='#86d4d8',
      "Macro_MT1H"='#9467BD')
  }
  
  color_B <- c('Naive B'="#1f77b4",'Memory B'="#aec7e8",
               'Plasma'="#2ca02c",'Plasmablast'="#98df8a")
  color_NK <- c('NK_CD16'="#1F78B4",'NK_CD56'="#98D083")
  color_DC <- c("cDC1"='#f7ce9c',"cDC2"='#3478b3',
                'cDC3_LAMP3'='#b83570',"pDC"='#a7cfe7')
  
}

## data load
{
  # object
  object <- readRDS('object.rds')
  # object <- readRDS(file.path(afterN_path,'Figure\\write\\script\\Data','object.rds'))
  
  # object_metadata
  all.metadata <- data.table(cellid=colnames(object),object@meta.data)
  
  #AUCell score
  load(file.path(afterN_path,'AUCell_score_all.rda'))
  
  # T
  object_T <- subset(object,subset= celltype_main=='T')
  object_T.meta <- all.metadata[match(colnames(object_T),cellid),]
  
  #CytoTRACE
  GeneCounts <- as.matrix(object_T@assays$RNA@counts)
  iOrd <- rowSums(GeneCounts>0)
  GeneCounts <- GeneCounts[iOrd>10,]
  object.CytoTRACE <- CytoTRACE(GeneCounts, enableFast = TRUE, ncores = 1, subsamplesize = 1000)
  
  # TCR
  load(file.path(fig_tcr,"productive filtered TCR.contigs and clonotypes for all samples.rda"))
  TCR_meta <- left_join(object_T.meta,cell2clone[,-c('sample','Patient')],by=c('cellid'='barcode'))
  TCR_meta <- TCR_meta[!is.na(sample_cloneID),]
  
  # CD8
  object_cd8 <- subset(object,subset= CellType_2=='CD8+ cells')
  object_cd8.meta <- all.metadata[match(colnames(object_cd8),cellid),]
  
  # cd4
  object_cd4 <- subset(object,subset= CellType_2=='CD4+ cells')
  object_cd4.meta <- all.metadata[match(colnames(object_cd4),cellid),]
  
  # NK
  object_NK <- subset(object,subset= CellType_2=='NK')
  object_NK_meta <- all.metadata[match(colnames(object_NK),cellid),]
  
  # B
  object_B <- subset(object,subset= CellType_2=='B Lineage')
  object_B_meta <- all.metadata[match(colnames(object_B),cellid),]
  
  # Myeloid
  object_Myeloid <- subset(object,subset= CellType_2%in%c('Myeloid','pDC'))
  object_Myeloid_meta <- all.metadata[match(colnames(object_Myeloid),cellid),]
  
  # Macrophage
  object_macro <- subset(object,subset= CellType_umap1%in%c('Macrophage'))
  object_macro.meta <- all.metadata[match(colnames(object_macro),cellid),]
  
  #DC
  object_DC <- subset(object,subset= CellType_umap1%in%c('DC','pDC'))
  object_DC_meta <- all.metadata[match(colnames(object_DC),cellid),]
  
}

## AUcell
{
  #ranking
  use_mat <- GetAssayData(object)
  cells_rankings <- AUCell_buildRankings(use_mat, nCores=8, plotStats=TRUE)
  
  #scores
  NeoTCR_CD8 <- AUCell_calcAUC(NeoTCR_CD8, cells_rankings)
  cytotoxicity_score <- AUCell_calcAUC(sig_zhangNK_cytotoxicity, cells_rankings)
  ISG_score <- AUCell_calcAUC(sig_ISG, cells_rankings)
  
  msigdb_go_bp_b_activate <- msigdb_go_bp[gs_name=='GOBP_B_CELL_ACTIVATION',gene_symbol]
  b_activate <- AUCell_calcAUC(msigdb_go_bp_b_activate, cells_rankings)
  
  msigdb_go_bp_antigen <- msigdb_go_bp[gs_name=='GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION',gene_symbol]
  msigdb_go_bp_antigen.score <- AUCell_calcAUC(msigdb_go_bp_antigen, cells_rankings)
  
  AUC_res <- data.table(cellid=colnames(object),
                        NeoTCR_CD8=assay(NeoTCR_CD8)[1,],
                        ISG_score=assay(ISG_score)[1,],
                        b_activate=assay(b_activate)[1,],
                        msigdb_go_bp_antigen=assay(msigdb_go_bp_antigen.score)[1,])
}





