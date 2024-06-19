
##Fig4A.
DimPlot(object_Myeloid,cols=color_Myeloid,group.by = 'CellType_n',label =T,label.size = 3,pt.size = 2,raster=T,shuffle=T)+ggtitle("")+NoAxes()

##Fig4B.
{
  DEG_table_cell <- data.table()
  for (i in unique(object_Myeloid_meta$CellType_n)) {
    tmp_cell <- object_Myeloid_meta[CellType_n==i,cellid]
    load(file.path(fig_all,'celltyp_n_diff',paste0(i,' R vs NR diff.rda')))
    diff_gene <- use_diff[p_val_adj<0.05,gene]
    tmp_mat <- GetAssayData(object_Myeloid)[diff_gene,tmp_cell]
    tmp_index <- colSums(tmp_mat>0) 
    DEG_table_cell <- rbind(DEG_table_cell,data.table(cellid=names(tmp_index),DEG=tmp_index))
  }
  object_Myeloid$DEG <- DEG_table_cell[match(colnames(object_Myeloid),cellid)]$DEG
  FeaturePlot(object_Myeloid,features = 'DEG',pt.size = 1.5,raster = T)+
    scale_color_gradientn(colours = c(colorRampPalette(c("#fdf2b5","#F67B51"))(10),
                                      colorRampPalette(c("#F67B51","#A30023"))(8)))+NoAxes()
}

##Fig4C.
{
  use_diff <- FindMarkers(object_macro,ident.1 = 'R',ident.2 = 'NR',group.by = 'response',test.use = "wilcox",
                          min.pct = 0.05,logfc.threshold = 0)
  use_diff$gene <- rownames(use_diff)
  use_diff <- data.table(use_diff)
  setorder(use_diff,p_val_adj)
  
  bp1 <- enrichGO(use_diff[avg_log2FC>0.25&p_val_adj<0.05,gene],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
  bp2 <- enrichGO(use_diff[avg_log2FC<(-0.25)&p_val_adj<0.05,gene],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
  up_select <- c('cytokine-mediated signaling pathway',
                 'pattern recognition receptor signaling pathway',
                 'leukocyte cell-cell adhesion',
                 'toll-like receptor signaling pathway',
                 'regulation of adaptive immune response',
                 'response to interferon-gamma')
  down_select <- c(
    'oxidative phosphorylation',
    'cellular respiration',
    'ATP metabolic process',
    'ATP biosynthetic process')
  plot_data <- rbind(data.table(data.table(bp1@result)[Description%in%up_select,],label='R'),
                     data.table(data.table(bp2@result)[Description%in%down_select,],label='NR'))
  plot_data$q <- -log10(plot_data$qvalue)
  plot_data$Description <- factor(plot_data$Description,levels = rev(plot_data$Description))
  plot_data[label=='NR','q'] <- -(plot_data[label=='NR',q])
  up <- plot_data[label=='R',]
  down <- plot_data[label=='NR',]
  
  ggplot(plot_data,aes(x =q, y = Description, fill = label)) + 
    geom_col(alpha=0.6) +
    theme_bw()+
    scale_fill_manual(values = color_pdsd)+
    geom_text(data = up,aes(x = -0.2, y = Description, label = Description),size = 2,hjust = 1)+ 
    geom_text(data = down,aes(x = 0.2, y = Description, label = Description),size = 2,hjust = 0)+
    scale_x_continuous(breaks = c(-20,-10,0,10,20))+
    labs(x = '-log10(qvalue)', y = ' ', title = 'Macrophage R vs NR Enriched Pathway') +
    theme(axis.text.x=element_text(size = 6,colour = 'black'),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(color = 'black',linewidth  = 0.235),
          axis.title.y = element_text(size = 6),
          axis.title.x = element_text(size = 6),
          plot.title = element_text(hjust=0.5,size = 8),
          legend.position = 'right',
          legend.text = element_text(size = 6),
          legend.title = element_blank(),
          legend.key.size = unit(0.1, "inches")
    )
}

##Fig4D.
{
  ## diff R vs NR
  {
    use_diff <- FindMarkers(object_macro,ident.1 = 'R',ident.2 = 'NR',group.by = 'response',test.use = "wilcox",
                            min.pct = 0.05,logfc.threshold = 0)
    use_diff$gene <- rownames(use_diff)
    use_diff <- data.table(use_diff)
    setorder(use_diff,p_val_adj)
    use_diff_R_NR <- use_diff
  }
  
  ## diff M1 vs M2
  {
    m1_gene <- intersect(sig_M1,rownames(object_macro))
    m2_gene <- intersect(sig_M2,rownames(object_macro))
    use_gene <- intersect(c(sig_M1_zhang,sig_M2_zhang),rownames(object_macro))
    mat_mac <- GetAssayData(object_Myeloid)[use_gene,]
    res1 <- apply(mat_mac, 2, function(x){t.test(x[m1_gene],x[m2_gene])$statistic})
    res_t <- data.table(cellid=names(res1),res1)
    
    use_data <- object_macro_meta
    use_data$M1vsM2_new <- res_t[match(use_data$cellid,cellid)]$res1
    use_data$M1M2_label <- 'M1/M2 high'
    use_data[M1vsM2_new<0]$M1M2_label <- 'M1/M2 low'
    object_macro$M1M2_label <- use_data$M1M2_label
    
    use_diff <- FindMarkers(object_macro,ident.1 = 'M1/M2 high',ident.2 = 'M1/M2 low',group.by = 'M1M2_label',
                            test.use = "wilcox",min.pct = 0.05,logfc.threshold = 0)
    use_diff$gene <- rownames(use_diff)
    use_diff_M1_M2 <- use_diff %>% data.table()
  }
  
  ## diff ISG_hi vs ISG_low
  {
    use_data <- object_macro_meta
    use_data <- left_join(use_data,AUC_res)
    use_data$ISG_score_scale <- scale(use_data$ISG_score)
    use_data$ISG_label <- 'ISG high'
    use_data[ISG_score_scale<0]$ISG_label <- 'ISG low'
    
    object_macro$ISG_label <- use_data$ISG_label
    
    use_diff <- FindMarkers(object_macro,ident.1 = 'ISG high',ident.2 = 'ISG low',group.by = 'ISG_label',
                            test.use = "wilcox",min.pct = 0.05,return.thresh = 1,logfc.threshold = 0)
    use_diff$gene <- rownames(use_diff)
    use_diff_ISG <- use_diff %>% data.table()
  }
  
  colnames(use_diff_R_NR) <- paste0('RvsNR_',colnames(use_diff_R_NR))
  colnames(use_diff_R_NR)[6] <- 'gene'
  colnames(use_diff_M1_M2) <- paste0('M1vsM2_',colnames(use_diff_M1_M2))
  colnames(use_diff_M1_M2)[6] <- 'gene'
  colnames(use_diff_ISG) <- paste0('ISG_',colnames(use_diff_ISG))
  colnames(use_diff_ISG)[6] <- 'gene'
  
  plot_data <- inner_join(use_diff_R_NR,use_diff_M1_M2)
  plot_data <- inner_join(plot_data,use_diff_ISG)
  plot_data$color_label <- plot_data$ISG_avg_log2FC
  
  colour_bk <- c(colorRampPalette(c("#2166ac","#d1e5f0"))(83),
                 colorRampPalette(c("#d1e5f0","#f7f7f7"))(15),
                 colorRampPalette(c("#f7f7f7","#fddbc7"))(15),
                 colorRampPalette(c("#fddbc7","#b2182b"))(84))
  
  ggplot(plot_data,aes(x=RvsNR_avg_log2FC,y=M1vsM2_avg_log2FC))+
    geom_point(aes(color=color_label),size=0.5)+
    scale_color_gradientn(name="Fold change (log2)",colours = colour_bk,limits=c(-1,1))+
    geom_hline(yintercept = 0, linetype = "dashed", color = "#999999",linewidth  = 0.235)+
    geom_vline(xintercept = 0, linetype = "dashed", color = "#999999",linewidth  = 0.235)+
    ggrepel::geom_text_repel(aes(label = label),data = plot_data,color="black",
                             size=2,segment.size=0.25,
                             min.segment.length = 0, 
                             segment.color='black',max.overlaps = Inf)+
    theme_classic()+ggtitle('')+ylab('Expression fold change\n(log2, M1/M2 hi vs. M1/M2 low)')+
    xlab('Expression fold change\n(log2, R vs. NR)')+labs(fill = "")+
    theme(axis.text.x=element_text(size = 6),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.y = element_text(size = 6),
          plot.title = element_text(hjust=0.5,size = 8),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.key.size = unit(0.1, "inches"),
          legend.position=c(0.8,0.1),
          axis.line = element_line(linewidth  = 0.235), 
          axis.ticks = element_line(linewidth  = 0.235)
    )
}

##Fig4E.
{
  object_Myeloid$M1vsM2_test <- res_t[match(colnames(object_Myeloid),cellid)]$res1
  
  FeaturePlot(object_Myeloid,cells=use_data[response=='R',cellid],features = 'M1vsM2_test',pt.size = 3,raster = T)+
    scale_color_gradientn(colours = c(colorRampPalette(c("#313695","#abd9e9"))(10),
                                      colorRampPalette(c("#abd9e9","#f0f0f0"))(3),
                                      colorRampPalette(c("#f0f0f0","#fee090"))(3),
                                      colorRampPalette(c("#fee090","#a50026"))(10)))+
    NoAxes()+ggtitle('M1/M2 score in Responders')
  
  FeaturePlot(object_Myeloid,cells=use_data[response=='NR',cellid],features = 'M1vsM2_test',pt.size = 2.5,raster = T)+
    scale_color_gradientn(colours = c(colorRampPalette(c("#313695","#abd9e9"))(10),
                                      colorRampPalette(c("#abd9e9","#f0f0f0"))(3),
                                      colorRampPalette(c("#f0f0f0","#fee090"))(5),
                                      colorRampPalette(c("#fee090","#a50026"))(10)))+
    NoAxes()+ggtitle('M1/M2 score in Nonresponders')
  
  object_Myeloid$ISG_score <- use_data$ISG_score
  object_Myeloid$ISG_score <- use_data$ISG_score
  FeaturePlot(object_Myeloid,cells=use_data[response=='R',cellid],features = 'ISG_score',pt.size = 3,raster = T)+
    scale_color_gradientn(colours = c(
      colorRampPalette(c("#363d99","#a0cae1"))(10),
      colorRampPalette(c("#a0cae1","#fedf90"))(5),
      colorRampPalette(c("#fedf90","#a8092a"))(10)))+
    NoAxes()+ggtitle('ISG score in Responders')
  
  FeaturePlot(object_Myeloid,cells=use_data[response=='NR',cellid],features = 'ISG_score',pt.size = 2.5,raster = T)+
    scale_color_gradientn(colours = c(
      colorRampPalette(c("#363d99","#a0cae1"))(10),
      colorRampPalette(c("#a0cae1","#fedf90"))(5),
      colorRampPalette(c("#fedf90","#a8092a"))(10)))+
    NoAxes()+ggtitle('ISG score in Nonresponders')
}

##Fig4G.
{
  tmp_cell <- c(object.meta[CellType_2%in%c('CD4+ cells','CD8+ cells'),cellid],object.meta[CellType_umap1%in%c('Macrophage'),cellid])
  use_obj <- subset(object,cells=tmp_cell)
  cellchat <- createCellChat(object = use_obj@assays$RNA@data,meta = use_obj@meta.data,group.by = 'CellType_n')
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  netVisual_aggregate(cellchat, signaling = 'CXCL', layout = "chord",color.use=color_tmp)
}

##Fig4H.
{
  load(file = 'Data/GO30140_IMbrave150/expr_process.rda')
  use_cli <- EGA_cli[drug!='Sorafenib'&treat=='Pre-treatment',]
  use_mat_tpm <- EGA_tpm[,use_cli$sampleid]
  
  ligand_gene <- c('CXCL9','CXCL10','CXCL11',"CXCL13",'CCL2','CCL3','CCL4','CCL5','CCL19')
  recept_gene <- c("CCRL2","CXCR3","CCR2","CCR3","CCR4","CCR5","CCR8")
  LR <- read.xlsx(file.path('combined Ligand_receptor pairs.xlsx')) %>% data.table()
  lr_pairs <- data.table()
  for (i in ligand_gene) {
    tmp <- LR[Ligand.ApprovedSymbol==i&Receptor.ApprovedSymbol%in%recept_gene,]
    lr_pairs <- rbind(lr_pairs,tmp[,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')])
  }
  colnames(lr_pairs) <- c('ligand','receptor')
  
  LR_mean_sample <- c()
  for(s in colnames(use_mat_tpm)){
    LR_mean <- c()
    for(lr in 1:nrow(lr_pairs)){
      lr_mean <-ifelse(lr_pairs[lr,ligand] %in% row.names(use_mat_tpm) & lr_pairs[lr,receptor] %in% row.names(use_mat_tpm),
                       mean(use_mat_tpm[lr_pairs[lr,ligand],s],use_mat_tpm[lr_pairs[lr,receptor],s],na.rm=T),NA)
      LR_mean <- c(LR_mean,lr_mean)
    }
    LR_mean_sample  <- c(LR_mean_sample,mean(LR_mean,na.rm=T))
  }
  use_cli$use_score <- LR_mean_sample
  plot_data <- use_cli[response%in%c('CR','PR','SD','PD')]
  plot_data$response <- factor(plot_data$response,levels = c('CR','PR','SD','PD'))
  ggplot(plot_data,aes(response,use_score))+
    geom_boxplot(aes(fill=response),width=0.5,color='black',outlier.shape = NA,linewidth  = 0.235)+
    geom_jitter(fill='white',color='black',shape=21,width =0.2,size=0.8,stroke = 0.1)+
    scale_fill_manual(values=c("#E7F0F9","#ABD0E6",'#3D91C8','#9267A9'))+
    stat_compare_means(comparisons=list(c('CR','PD'),c('PR','PD'),c('SD','PD')) ,method ="wilcox.test",size=2)+
    theme_classic()+ggtitle('GO30140&IMbrave150 cohort')+ylab('Ligand-receptor interaction intensity')+xlab('')+
    theme(
      axis.text.x=element_text(size = 6,colour = 'black'),
      axis.title.y = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      axis.text.y = element_text(size = 6,colour = 'black'),
      axis.line = element_line(linewidth  = 0.235), 
      axis.ticks = element_line(linewidth  = 0.235),
      plot.title = element_text(hjust=0.5,size = 8),
      legend.position = 'none'
    )
}

##FigS4A.
{
  markers <- list(combine=c(c('FCN1','VCAN'),
                            c('CD14','FCGR3A','CD68','CD163'),
                            c("APOE",'SPP1',"C1QC",'TREM2'),
                            c('NFKB1','IL15','CD44'),
                            c('CXCL10','CXCL11'),
                            c('LGALS3','CXCL3','ACP5'),
                            c('MT1H','MT1G'),
                            c('S100A8','S100A9','IL1R2','CSF3R'),
                            c('CADM1','CLEC9A','CLNK'),
                            c('CD1C','CLEC10A','FCGR2B'),
                            c('LAMP3','FSCN1','CD40'),
                            c('IL3RA',"GZMB",'JCHAIN')))
  DotPlot(object_Myeloid, features = markers,group.by = "CellType_n",dot.scale=2)+
    scale_y_discrete(limits=rev(names(color_Myeloid)))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5,size=6,color = 'black'),
          axis.text.y=element_text(angle=0,hjust = 1,vjust=0.5,size=6,color = 'black'),
          legend.title=element_text(size = 6),
          legend.text=element_text(size = 6),
          legend.key.size = unit(0.1, "inches"),
          panel.border=element_rect(size=0.5),
          strip.text=element_text(size=6),
          strip.background=element_blank())+
    scale_color_gradientn(colours = c("#2166ac",'#abd1e5',"#f0f3f4",'#f9bfa3',"#b31b2c"))+
    labs(x=NULL,y=NULL)+
    guides(size=guide_legend(order=3))
}

##FigS4B. 
{
  markers <- list(
    M1=c("IL23A","TNF","CXCL9","CXCL10","CXCL11","CD86","IL1A","IL1B","IL6","CCL5",
         "IRF5","IRF1","CD40","IDO1","KYNU","CCR7"),
    M2=c("IL4R","CCL4","CCL13","CCL20","CCL17","CCL18","CCL22","CCL24","LYVE1",
         "VEGFA","VEGFB","CTSA","CTSB","CTSC","CTSD",
         "TGFB1","MMP14","MMP19","MMP9","CLEC7A","WNT7B",
         "FASLG","TNFSF12","TNFSF8","CD276","MSR1","FN1","IRF4"))
  DotPlot(object_macro, features = markers,group.by = "CellType_n",dot.scale=2)+
    theme_bw()+
    scale_y_discrete(limits=rev(names(color_macro)))+
    theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5,size=6,colour = 'black'),
          axis.text.y=element_text(angle=0,hjust = 1,vjust=0.5,size=6,colour = 'black'),
          legend.title=element_text(size = 6),
          legend.text=element_text(size = 6),
          legend.key.size = unit(0.1, "inches"),
          panel.border=element_rect(size=0.5),
          strip.text=element_text(size=8),
          strip.background=element_blank())+
    scale_color_gradientn(colours = c("#2166ac",'#abd1e5',"#f0f3f4",'#f9bfa3',"#b31b2c"))+
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
}

##FigS4C.
{
  plot_data <- left_join(object_macro_meta,AUC_res)
  plot_data$M1vsM2 <- res_t[match(plot_data$cellid,cellid)]$res1
  plot_data$CellType_n <- factor(plot_data$CellType_n,levels = c('Macro_CXCL10','Macro_NFKB1','Macro_LGALS3','Macro_C1QC','Macro_MT1H'))
  
  p1 <- ggplot(plot_data,aes(CellType_n,ISG_score,fill=response,color=response))+
    scale_fill_manual(values = color_pdsd_fill)+
    scale_color_manual(values = color_pdsd)+
    geom_boxplot(width=0.7,linewidth=0.235,outlier.shape = NA)+
    ylab('ISG score')+xlab('')+theme_classic()+ggtitle('Macrophage (n = 19,463)')+
    stat_compare_means(aes(group = response,label = paste0(..p.format..)),method = "wilcox.test",hide.ns = F,size=2)+
    theme(
      axis.text.x=element_blank(),
      axis.title.y = element_text(size = 6),
      axis.title.x = element_text(size = 6),
      axis.text.y = element_text(size = 6,colour = 'black'),
      axis.line = element_line(linewidth  = 0.235), 
      axis.ticks = element_line(linewidth  = 0.235),
      plot.title = element_text(hjust=0.5,size = 6),
      legend.position = 'none',
      strip.background=element_blank()
    )
  p2 <- ggplot(plot_data,aes(CellType_n,M1vsM2,fill=response,color=response))+
    scale_fill_manual(values = color_pdsd_fill)+
    scale_color_manual(values = color_pdsd)+
    geom_boxplot(width=0.7,linewidth=0.235,outlier.shape = NA)+
    ylab('M1 likeness score')+xlab('')+theme_classic()+ggtitle('')+
    stat_compare_means(aes(group = response,label = paste0(..p.format..)),method = "wilcox.test",hide.ns = F,size=2)+
    theme(
      axis.text.x=element_text(angle=45,vjust=1,hjust=1,size = 6,colour = 'black'),
      axis.title.y = element_text(size = 6),
      axis.title.x = element_text(size = 6),
      axis.text.y = element_text(size = 6,colour = 'black'),
      axis.line = element_line(linewidth  = 0.235), 
      axis.ticks = element_line(linewidth  = 0.235),
      plot.title = element_text(hjust=0.5,size = 6),
      legend.position = 'none',
      strip.background=element_blank()
    )
  p1+p2+plot_layout(ncol = 1)
}

##FigS4D.
{
  plot_data <- object_Myeloid_meta[,.(count=.N),by=.(sample,response,CellType_n)] %>% arrange(sample)
  plot_data <- plot_data[,.(ratio=count/sum(count),allcount=sum(count),count,CellType_n,response),by=.(sample)]
  ggplot(plot_data[CellType_n=='Macro_CXCL10'],aes(response, ratio))+
    geom_boxplot(aes(fill=response,color=response),width=0.4,outlier.shape = NA)+
    scale_fill_manual(values=color_pdsd_fill)+
    scale_color_manual(values=color_pdsd)+
    geom_jitter(fill='white',color='black',shape=21,width =0.2,size=1,stroke = 0.2)+
    ylab('% of Myeloid cells')+xlab('')+ggtitle('Macro_CXCL10 Cells')+
    stat_compare_means(comparisons=list(c('R','NR')),method ="wilcox.test",size=2)+
    theme_classic() +scale_y_continuous(labels = scales::percent)+
    theme(
      axis.text.x=element_text(size = 6,colour = 'black'),
      axis.title.y = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      axis.text.y = element_text(size = 6,colour = 'black'),
      axis.line = element_line(linewidth  = 0.235),
      axis.ticks = element_line(linewidth  = 0.235),
      plot.title = element_text(hjust=0.5,size = 8),
      legend.position = 'none'
    )
}

##FigS4E.
{
  object_macro <- subset(object_Myeloid, subset=CellType_umap1=='Macrophage')
  object_macro$use_label <- 'other'
  object_macro@meta.data[object_macro$CellType_n=='Macro_CXCL10',]$use_label <- 'Macro_CXCL10'
  
  use_diff <- FindMarkers(object_macro,ident.1 = 'Macro_CXCL10',ident.2 = 'other',group.by = 'use_label',test.use = "wilcox",
                          min.pct = 0.1,return.thresh = 0.01,logfc.threshold = 0.25)
  use_diff$gene <- rownames(use_diff)
  use_diff <- data.table(use_diff)
  setorder(use_diff,-avg_log2FC)
  
  bp1 <- enrichGO(use_diff[avg_log2FC>0,gene],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
  bp2 <- enrichGO(use_diff[avg_log2FC<0,gene],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
  use_go_1 <- bp1@result[bp1@result$p.adjust < 0.05,] %>% data.table()
  use_go_2 <- bp2@result[bp2@result$p.adjust < 0.05,] %>% data.table()
  
  up_select <- c('cytokine-mediated signaling pathway',
                 'response to interferon-gamma',
                 'positive regulation of cytokine production',
                 'response to interferon-beta',
                 'cellular response to chemokine',
                 'myeloid leukocyte migration')
  down_select <- c('myeloid cell differentiation',
                   'receptor-mediated endocytosis',
                   'cellular transition metal ion homeostasis',
                   'negative regulation of immune system process')
  plot_data <- rbind(data.table(use_go_1[Description%in%up_select,],label='Macro_CXCL10'),
                     data.table(use_go_2[Description%in%down_select,],label='Macro_other'))
  plot_data$q <- -log10(plot_data$qvalue)
  plot_data$label <- factor(plot_data$label,levels = c('Macro_CXCL10','Macro_other'))
  plot_data$Description <- factor(plot_data$Description,levels = rev(plot_data$Description))
  plot_data[label=='Macro_other','q'] <- -(plot_data[label=='Macro_other',q])
  up <- plot_data[label=='Macro_CXCL10',]
  down <- plot_data[label=='Macro_other',]
  
  ggplot(plot_data,aes(x =q, y = Description, fill = label))+ 
    geom_col(alpha=0.6)+theme_bw()+
    scale_fill_manual(values = c('#EFC675','#63B5CB'))+
    geom_text(data = up,aes(x = -0.2, y = Description, label = Description),size = 2,hjust = 1)+ 
    geom_text(data = down,aes(x = 0.2, y = Description, label = Description),size = 2,hjust = 0)+
    xlim(-30,30)+labs(x = '-log10(qvalue)', y = '', title = '')+
    theme(axis.text.x=element_text(size = 6,colour = 'black'),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(color = 'black',linewidth  = 0.235),
          axis.title.y = element_text(size = 6),
          axis.title.x = element_text(size = 6),
          plot.title = element_text(hjust=0.5,size = 8)
    )
}



