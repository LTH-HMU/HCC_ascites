
##Fig5A.
{
  load(file = file.path('T R vs NR diff gene.rda'))
  T_top50 <- use_diff[1:50,gene]
  
  load(file = file.path('NK R vs NR diff gene.rda'))
  NK_top50 <- use_diff[1:50,gene]
  
  load(file = file.path('All B R vs NR diff gene.rda'))
  B_top50 <- use_diff[1:50,gene]
  
  load(file = file.path('DC R vs NR diff gene.rda'))
  DC_top50 <- use_diff[1:50,gene]
  
  load(file = file.path('macro R vs NR diff.rda'))
  macro_top50 <- use_diff[1:50,gene]
  
  library(VennDiagram)
  plotdata <- list(T_top50,NK_top50,B_top50,DC_top50,macro_top50)
  col1 <- c('#5C79B3','#B4C2DA','#D6CEB0','#A5BA94','#D7AA8C')
  venn.diagram(x=plotdata,scaled = F,alpha= 0.5,lwd=1,lty=1,col=col1, 
               label.col = 'black',cex = 1,fontface = "bold",fill=col1,
               category.names = c("Up in Response T", "Up in Response NK","Up in Response B","Up in Response DC","Up in Response macro_CXCL10") , #标签名
               cat.dist = c(0.2, 0.2, 0.2, 0.2, 0.2),
               cat.pos = c(0, -10, 240, 120, 120),
               cat.cex = 1,cat.fontface = "bold",cat.col=col1,
               cat.default.pos = "outer",output=TRUE,filename=NULL
  )
}

##Fig5B.
{
  mat <- GetAssayData(object)[sig_ISG_select,]
  object.meta$tmp_label <- paste0(object.meta$CellType_n,'-',object.meta$response)
  plot_data <- groupMeans(mat,object.meta$tmp_label,na.rm = TRUE, sparse = T)
  plot_data <- t(scale(t(plot_data)))
  col_table <- data.table(colnames(plot_data))
  col_table$label <- str_split_fixed(col_table$V1,'-',n=2)[,2]
  col_table$label <- factor(col_table$label,levels = c('R','NR'))
  bk <- c(seq(-2,-0.1,by=0.02),seq(0,2,by=0.02))
  colour_bk <- c(colorRampPalette(c("#2166ac","#d1e5f0"))(83),
                 colorRampPalette(c("#d1e5f0","#f7f7f7"))(15),
                 colorRampPalette(c("#f7f7f7","#fddbc7"))(15),
                 colorRampPalette(c("#fddbc7","#b2182b"))(84))
  
  col_bar = HeatmapAnnotation(sample_label = col_table$label,height = unit(2, "cm"),col = list('sample_label'=color_pdsd))
  Heatmap(plot_data,
          show_row_names = T,
          cluster_rows = F,
          clustering_method_rows='ward.D',
          clustering_method_columns='ward.D',
          row_names_side = 'right',
          col = circlize::colorRamp2(bk, colour_bk),
          column_names_gp = gpar(fontsize = 6),
          row_names_gp = gpar(fontsize = 6),
          top_annotation = col_bar,
          border = T,
          border_gp = gpar(col = 'grey80', lwd = 0.5),
          rect_gp  = gpar(col = 'grey80', lwd = 0.5),
          width = unit(12, "cm"), height = unit(4, "cm")
  )
}

##Fig5C.
{
  DimPlot(object_B,cols=color_B,group.by = 'CellType_n',label =F,label.size = 3,pt.size = 1,raster=F,shuffle=T)+NoAxes()
  
  plot_data <- left_join(object_B_meta,AUC_res)
  plot_data$ISG_score_scale <- scale(plot_data$ISG_score)
  plot_data$ISG_label <- 'ISG_high'
  plot_data[ISG_score_scale<0,]$ISG_label <- 'ISG_low'
  plot_data$ISG_label <- factor(plot_data$ISG_label,levels = c('ISG_high','ISG_low'))
  
  #ISG score
  {
    ggplot(plot_data,aes(CellType_n,ISG_score,fill=response,color=response))+
      scale_fill_manual(values = color_pdsd_fill)+
      scale_color_manual(values = color_pdsd)+
      geom_boxplot(width=0.8,linewidth=0.235,outlier.shape = NA)+
      ylab('ISG score')+xlab('')+theme_classic()+ggtitle('B cell lineage (n = 1,843)')+
      stat_compare_means(aes(group = response,label = paste0(..p.format..)),method = "wilcox.test",hide.ns = F,size=1.8)+
      scale_x_discrete(limits= names(color_B))+
      theme(
        axis.text.x=element_text(angle=45,vjust=1,hjust=1,size = 6,colour = 'black'),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 6,colour = 'black'),
        axis.line = element_line(linewidth  = 0.235),
        axis.ticks = element_line(linewidth  = 0.235),
        plot.title = element_text(hjust=0.5,size = 8),
        legend.position = 'none',
        strip.background=element_blank()
      )
  }
  
  #function enrichment
  {
    use_obj <- subset(object_B,cells=plot_data$cellid)
    use_obj$ISG_label <- plot_data$ISG_label
    use_diff <- FindMarkers(use_obj,ident.1 = 'ISG_high',ident.2 = 'ISG_low',group.by = 'ISG_label',test.use = "wilcox",
                            min.pct = 0.2,logfc.threshold = 0.25)
    use_diff$gene <- rownames(use_diff)
    use_diff <- data.table(use_diff)
    setorder(use_diff,-avg_log2FC)
    
    bp1 <- enrichGO(use_diff[avg_log2FC>0.25&p_val_adj<0.05,gene],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
    bp2 <- enrichGO(use_diff[avg_log2FC<(-0.25)&p_val_adj<0.05,gene],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
    up_select <- c('I-kappaB kinase/NF-kappaB signaling',
                   'regulation of innate immune response',
                   'pattern recognition receptor signaling pathway',
                   'response to type I interferon',
                   'B cell activation',
                   'regulation of B cell proliferation')
    down_select <- c(
      'ribosome biogenesis',
      'rRNA processing',
      'ribosome assembly','ribonucleoprotein complex assembly')
    plot_data <- rbind(data.table(data.table(bp1@result)[Description%in%up_select,],label='ISG_high'),
                       data.table(data.table(bp2@result)[Description%in%down_select,],label='ISG_low'))
    plot_data$q <- -log10(plot_data$qvalue)
    plot_data$Description <- factor(plot_data$Description,levels = rev(plot_data$Description))
    plot_data[label=='ISG_low','q'] <- -(plot_data[label=='ISG_low',q])
    up <- plot_data[label=='ISG_high',]
    down <- plot_data[label=='ISG_low',]
    
    ggplot(plot_data,aes(x =q, y = Description, fill = label))+
      geom_col(alpha=0.6)+theme_bw()+
      scale_fill_manual(values = color_ISG)+
      geom_text(data = up,aes(x = -0.2, y = Description, label = Description),size = 2,hjust = 1)+
      geom_text(data = down,aes(x = 0.2, y = Description, label = Description),size = 2,hjust = 0)+
      labs(x = '-log10(qvalue)', y = ' ', title = 'B lineage ISG_high vs ISG_low Enriched Pathway') +
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
  
  #B cell activation
  {
    ggplot(plot_data,aes(ISG_label, b_activate))+
      geom_boxplot(aes(fill=ISG_label),width=0.4,outlier.shape = NA,linewidth  = 0.235)+
      scale_fill_manual(values=color_ISG)+
      geom_jitter(fill='white',color='black',shape=21,width =0.1,size=0.4,stroke = 0.1)+
      ylab('B cell activation (GO:0042113)')+xlab('')+ggtitle('B lineage cells')+
      stat_compare_means(comparisons=list(c('ISG_high','ISG_low')),method ="wilcox.test",size=2)+
      theme_classic()+
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
}

##Fig5D.
{
  DimPlot(object_NK,cols=color_NK,group.by = 'CellType_n',label =F,label.size = 3,pt.size = 2,raster=T,shuffle=T)+NoAxes()
  
  plot_data <- left_join(object_NK_meta,AUC_res)
  plot_data$cytotoxicity_score_scale <- scale(plot_data$cytotoxicity_score)
  plot_data$ISG_score_scale <- scale(plot_data$ISG_score)
  plot_data$ISG_label <- 'ISG_high'
  plot_data[ISG_score_scale<0,]$ISG_label <- 'ISG_low'
  
  #ISG score
  {
    ggplot(plot_data,aes(CellType_n,ISG_score,fill=response,color=response))+
      scale_fill_manual(values = color_pdsd_fill)+
      scale_color_manual(values = color_pdsd)+
      geom_boxplot(width=0.8,linewidth=0.235,outlier.shape = NA)+
      ylab('ISG score')+xlab('')+theme_classic()+ggtitle('NK cells (n = 15,615)')+
      stat_compare_means(aes(group = response,label = paste0(..p.format..)),method = "wilcox.test",hide.ns = F,size=1.8)+
      theme(
        axis.text.x=element_text(angle=45,vjust=1,hjust=1,size = 6,colour = 'black'),
        axis.title.y = element_text(size = 6),
        axis.title.x = element_text(size = 6),
        axis.text.y = element_text(size = 6,colour = 'black'),
        axis.line = element_line(linewidth  = 0.235),
        axis.ticks = element_line(linewidth  = 0.235),
        plot.title = element_text(hjust=0.5,size = 8),
        legend.position = 'none',
        strip.background=element_blank()
      )
  }
  
  #function enrichment
  {
    use_obj <- subset(object_NK,cells=plot_data$cellid)
    use_obj$ISG_label <- plot_data$ISG_label
    use_diff <- FindMarkers(use_obj,ident.1 = 'ISG_high',ident.2 = 'ISG_low',group.by = 'ISG_label',test.use = "wilcox",
                            min.pct = 0.2,logfc.threshold = 0.25)
    use_diff$gene <- rownames(use_diff)
    use_diff <- data.table(use_diff)
    
    bp1 <- enrichGO(use_diff[avg_log2FC>0.25&p_val_adj<0.05,gene],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
    bp2 <- enrichGO(use_diff[avg_log2FC<(-0.25)&p_val_adj<0.05,gene],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
    up_select <- c('response to interferon-gamma',
                   'regulation of innate immune response',
                   'response to type I interferon',
                   'cytokine-mediated signaling pathway',
                   'response to interferon-beta',
                   'lymphocyte differentiation')
    down_select <- c(
      'homotypic cell-cell adhesion',
      'postsynapse organization')
    plot_data <- rbind(data.table(data.table(bp1@result)[Description%in%up_select,],label='ISG_high'),
                       data.table(data.table(bp2@result)[Description%in%down_select,],label='ISG_low'))
    plot_data$q <- -log10(plot_data$qvalue)
    plot_data$Description <- factor(plot_data$Description,levels = rev(plot_data$Description))
    plot_data[label=='ISG_low','q'] <- -(plot_data[label=='ISG_low',q])
    up <- plot_data[label=='ISG_high',]
    down <- plot_data[label=='ISG_low',]
    
    ggplot(plot_data,aes(x =q, y = Description, fill = label)) +
      geom_col(alpha=0.6) +
      theme_bw()+
      scale_fill_manual(values = color_ISG)+
      geom_text(data = up,aes(x = -0.2, y = Description, label = Description),size = 2,hjust = 1)+ 
      geom_text(data = down,aes(x = 0.2, y = Description, label = Description),size = 2,hjust = 0)+ 
      labs(x = '-log10(qvalue)', y = ' ', title = 'NK cells ISG_high vs ISG_low Enriched Pathway') + 
      theme(axis.text.x=element_text(size = 6,colour = 'black'),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line.x = element_line(color = 'black',linewidth  = 0.235),
            # axis.text = element_text(size = 12)
            axis.title.y = element_text(size = 6),
            axis.title.x = element_text(size = 6),
            plot.title = element_text(hjust=0.5,size = 8)
      )
  }
  
  #Cytotoxicity
  {
    ggplot(plot_data,aes(ISG_label, cytotoxicity_score))+
      geom_boxplot(aes(fill=ISG_label),width=0.4,outlier.shape = NA,linewidth  = 0.235)+
      scale_fill_manual(values=color_ISG)+
      geom_jitter(fill='white',color='black',shape=21,width =0.1,size=0.4,stroke = 0.1)+
      ylab('Cytotoxicity score')+xlab('')+ggtitle('NK cells')+
      stat_compare_means(comparisons=list(c('ISG_high','ISG_low')),method ="wilcox.test",size=2)+
      theme_classic()+
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
}

##Fig5E.
{
  DimPlot(object_DC,cols=color_all_celltype1,group.by = 'CellType_n',label =T,label.size = 3,pt.size = 1,raster=F,shuffle=T)+NoAxes()
  plot_data <- left_join(object_DC_meta,AUC_res)
  plot_data$ISG_score_scale <- scale(plot_data$ISG_score)
  plot_data$ISG_label <- 'ISG_high'
  plot_data[ISG_score_scale<0,]$ISG_label <- 'ISG_low'
  
  #ISG score
  {
    ggplot(plot_data,aes(CellType_n,ISG_score,fill=response,color=response))+
      scale_fill_manual(values = color_pdsd_fill)+
      scale_color_manual(values = color_pdsd)+
      geom_boxplot(width=0.8,linewidth=0.235,outlier.shape = NA)+
      ylab('ISG score')+xlab('')+theme_classic()+ggtitle('DC (n = 974 )')+
      stat_compare_means(aes(group = response,label = paste0(..p.format..)),method = "wilcox.test",hide.ns = F,size=1.8)+
      scale_x_discrete(limits= c('cDC1','cDC2','cDC3_LAMP3','pDC'))+
      theme(
        axis.text.x=element_text(angle=45,vjust=1,hjust=1,size = 6,colour = 'black'),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 6,colour = 'black'),
        axis.line = element_line(linewidth  = 0.235),
        axis.ticks = element_line(linewidth  = 0.235),
        plot.title = element_text(hjust=0.5,size = 8),
        legend.position = 'none',
        strip.background=element_blank()
      )
  }
  
  #function enrichment
  {
    use_obj <- subset(object_DC,cells=plot_data$cellid)
    use_obj$ISG_label <- plot_data$ISG_label
    use_diff <- FindMarkers(use_obj,ident.1 = 'ISG_high',ident.2 = 'ISG_low',group.by = 'ISG_label',test.use = "wilcox",
                            min.pct = 0.2,logfc.threshold = 0.25)
    use_diff$gene <- rownames(use_diff)
    use_diff <- data.table(use_diff)
    
    bp1 <- enrichGO(use_diff[avg_log2FC>0.25&p_val_adj<0.05,gene],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
    bp2 <- enrichGO(use_diff[avg_log2FC<(-0.25)&p_val_adj<0.05,gene],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
    up_select <- c(
      'positive regulation of cytokine production',
      'mononuclear cell differentiation',
      'I-kappaB kinase/NF-kappaB signaling',
      'pattern recognition receptor signaling pathway',
      'leukocyte cell-cell adhesion',
      'antigen receptor-mediated signaling pathway')
    down_select <- c(
      'actin filament organization',
      'postsynaptic actin cytoskeleton organization',
      'regulation of cell morphogenesis')
    plot_data <- rbind(data.table(data.table(bp1@result)[Description%in%up_select,],label='ISG_high'),
                       data.table(data.table(bp2@result)[Description%in%down_select,],label='ISG_low'))
    plot_data$q <- -log10(plot_data$qvalue)
    plot_data$Description <- factor(plot_data$Description,levels = rev(plot_data$Description))
    plot_data[label=='ISG_low','q'] <- -(plot_data[label=='ISG_low',q])
    up <- plot_data[label=='ISG_high',]
    down <- plot_data[label=='ISG_low',]
    
    ggplot(plot_data,aes(x =q, y = Description, fill = label)) + 
      geom_col(alpha=0.6)+theme_bw()+
      scale_fill_manual(values = color_ISG)+
      geom_text(data = up,aes(x = -0.2, y = Description, label = Description),size = 2,hjust = 1)+ 
      geom_text(data = down,aes(x = 0.2, y = Description, label = Description),size = 2,hjust = 0)+ 
      labs(x = '-log10(qvalue)', y = ' ', title = 'DC ISG_high vs ISG_low Enriched Pathway') + 
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
  
  #Antigen processing and presentation
  {
    ggplot(plot_data,aes(ISG_label, msigdb_go_bp_antigen))+
      geom_boxplot(aes(fill=ISG_label),width=0.4,outlier.shape = NA,linewidth  = 0.235)+
      scale_fill_manual(values=color_ISG)+
      geom_jitter(fill='white',color='black',shape=21,width =0.1,size=0.4,stroke = 0.1)+
      ylab('Antigen processing and presentation (GO:0019882)')+xlab('')+ggtitle('DCs')+
      stat_compare_means(comparisons=list(c('ISG_high','ISG_low')),method ="wilcox.test",size=2)+
      theme_classic()+
      theme(
        axis.text.x=element_text(size = 6,colour = 'black'),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 6,colour = 'black'),
        axis.line = element_line(linewidth  = 0.235),
        axis.ticks = element_line(linewidth  = 0.235),
        plot.title = element_text(hjust=0.5,size = 8),
        legend.position = 'none')
  }
}

##FigS5A.
{
  markers <- list( t1=c(c("CD79A","CD79B","MS4A1","CD19","CD37"),
                        c('IL4R','TCL1A','FCER2'),
                        c('CD27','CD24','CD44'),
                        c("JCHAIN","MZB1",'CD38'),
                        c("MKI67","TOP2A","TUBB","STMN1","TYMS"),
                        c('NKG7'),
                        c('FCGR3A','KLRB1'),
                        c('NCAM1','SELL','IL7R','KLRC2'),
                        c('CADM1','CLEC9A','XCR1','DBN1','FKBP1B'),
                        c('CD1C','FCER1A','CLEC10A','FCGR2B'),
                        c('LAMP3','FSCN1','CD40'),
                        c('IL3RA',"TCF4","IRF7",'LILRA4')))
  
  DotPlot(use_obj, features = markers,group.by = "CellType_n",dot.scale=2)+
    theme_bw()+
    scale_y_discrete(limits=rev(c(names(color_B),names(color_NK),names(color_DC))))+
    theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5,size=6,color = 'black'),
          axis.text.y=element_text(angle=0,hjust = 1,vjust=0.5,size=6,color = 'black'),
          legend.title=element_text(size = 6),
          legend.text=element_text(size = 6),
          legend.key.size = unit(0.1, "inches"),
          panel.border=element_rect(size=0.5),
          strip.text=element_text(size=6),
          strip.background=element_blank())+
    scale_color_gradientn(colours = c("#2166ac",'#abd1e5',"#f0f3f4",'#f9bfa3',"#b31b2c"))+
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
}

##FigS5B.
{
  #Part 0:
  {
    use_data <- all.metadata
    use_data$use_celltype <- use_data$CellType_umap1
    use_data[celltype_main=='B Lineage',]$use_celltype <- 'B Lineage'
    use_data <- left_join(use_data,AUC_res)
    ISG_hi_meta <- data.table()
    ISG_low_meta <- data.table()
    for (i in unique(use_data$use_celltype)) {
      tmp <- use_data[use_celltype==i,]
      tmp$ISG_score_scale <- scale(tmp$ISG_score)
      tmp$ISG_label <- 'ISG_high'
      tmp[ISG_score_scale<0,]$ISG_label <- 'ISG_low'
      ISG_hi_meta <- rbind(ISG_hi_meta,tmp[ISG_label=='ISG_high',])
      ISG_low_meta <- rbind(ISG_low_meta,tmp[ISG_label=='ISG_low',])
    }
    obj_ISG_hi <- subset(object,cells=ISG_hi_meta$cellid)
    obj_ISG_hi@meta.data <- ISG_hi_meta
    obj_ISG_low <- subset(object,cells=ISG_low_meta$cellid)
    obj_ISG_low@meta.data <- ISG_low_meta
  }
  
  #Part I: 
  {
    cellchat_isg_hi <- createCellChat(object = obj_ISG_hi@assays$RNA@data,meta = obj_ISG_hi@meta.data,group.by = 'use_celltype')
    cellchat_isg_low <- createCellChat(object = obj_ISG_low@assays$RNA@data,meta = obj_ISG_low@meta.data,group.by = 'use_celltype')
    cellchat_step1 <- function(cellchat){
      cellchat@DB <- CellChatDB.human
      cellchat <- subsetData(cellchat)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)
      cellchat <- computeCommunProb(cellchat)
      cellchat <- filterCommunication(cellchat, min.cells = 10)
      cellchat <- computeCommunProbPathway(cellchat)
      cellchat <- aggregateNet(cellchat)
      return(cellchat)
    }
    cellchat_isg_hi <- cellchat_step1(cellchat_isg_hi)
    cellchat_isg_low <- cellchat_step1(cellchat_isg_low)
  }
  
  #Part II:
  {
    object.list <- list(isg_hi = cellchat_isg_hi, isg_low = cellchat_isg_low)
    cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
  }
  
  cellchat_color <- c('CD4_Naive'='#c8e6c9',
                      'CD4_Mem'='#9EBEE2',
                      'CD4_GZMK'= '#BCBD22',
                      'CD4_EFF'='#ffa58d',
                      'CD4_Tfh'='#9E76C1',
                      'CD4_Treg'='#FF9332',
                      'CD8_Naive'="#0a7b37",
                      'CD8_Tcm'="#3BBEB2",
                      'CD8_GZMK'="#5D9FCF",
                      'CD8_Active'='#f9be25',
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
                      "B Lineage"='#d5c8a0',
                      "Epithelial"='#9668BC')
  netVisual_aggregate(cellchat_isg_hi, signaling = 'CXCL', layout = "chord",color.use=cellchat_color,cell.order = names(cellchat_color),
                      edge.width.max = 1,vertex.size.max = 0.5,edge.label.cex = 0.7,vertex.label.cex=0.8,
                      arrow.size = 0.1)
  
  netVisual_aggregate(cellchat_isg_low, signaling = 'CXCL', layout = "chord",color.use=cellchat_color,cell.order = names(cellchat_color),
                      edge.width.max = 1,vertex.size.max = 0.5,edge.label.cex = 0.7,vertex.label.cex=0.8,
                      arrow.size = 0.1)
}