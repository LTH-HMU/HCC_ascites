
# fig1B
DimPlot(object, reduction = "umap",cols=color_umap1,pt.size = 1.2,group.by = "CellType_umap1",label =F,raster=T)+ggtitle("")+NoLegend()

# fig1C
{
  bk <- c(seq(-2,-0.1,by=0.02),seq(0,2,by=0.02))
  colour_bk <- c(colorRampPalette(c("#2166ac","#d1e5f0"))(83),
                 colorRampPalette(c("#d1e5f0","#f7f7f7"))(15),
                 colorRampPalette(c("#f7f7f7","#fddbc7"))(15),
                 colorRampPalette(c("#fddbc7","#b2182b"))(84))
  object@active.ident <- factor(object$CellType_umap1,levels = names(color_umap1))
  diff_umap1 <- FindAllMarkers(object, only.pos = T, min.pct = 0.2,test.use = "wilcox",return.thresh = 0.01,logfc.threshold = 0.25)
  diff_top50 <- use_diff %>% group_by(cluster) %>% slice_head(n=50) %>% data.table()
  plot_data <- GetAssayData(object,slot = "data")[diff_top50$gene,]
  plot_data <- groupMeans(plot_data,object$CellType_umap1,na.rm = T,sparse = T)
  plot_data <- plot_data[,names(color_umap1)]
  plot_data <- t(scale(t(plot_data)))
  Heatmap(plot_data,
          cluster_rows=F,cluster_columns = F,
          show_column_names = T,show_row_names = F,
          row_names_gp = gpar(fontsize = 6),
          row_names_side = 'right',column_names_gp = gpar(fontsize = 6),
          col = circlize::colorRamp2(bk, colour_bk),
          row_split = anno_table$type,column_split = col_ano$type,
          row_gap = unit(0.5, 'mm'),column_gap = unit(0.5, 'mm'),
          width = unit(7, "cm"), height = unit(8, "cm"))
}

# fig1D
{
  DimPlot(object, reduction = "umap",cols=color_umap1,pt.size = 2,group.by = "CellType_umap1",label =F,raster=T,split.by = 'combine_label')+ggtitle("")+NoLegend()+NoAxes()
  
  library(viridis)
  metaData <- all.metadata
  umaps <- Embeddings(object,reduction = "umap")
  metaData <- cbind(metaData,umaps)
  for (i in unique(metaData$combine_label)) {
    x <- metaData[metaData$combine_label == i,]
    dim(x)[1]
    if (dim(x)[1]>8000) {
      g.overlay <- ggplot(data = x,aes(x = UMAP_1, y = UMAP_2)) +
        stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F)+
        geom_point(data = x[sample(seq_along(x$cellid),8000,replace = F),],
                   aes(x = UMAP_1, y = UMAP_2),color = 'white',size = .005)+
        scale_fill_viridis(option="A")+theme_classic()+NoAxes()+NoLegend()
    }else{
      g.overlay <- ggplot(data = x,aes(x = UMAP_1, y = UMAP_2)) +
        stat_density_2d(aes(fill = ..density..), geom = "raster",contour = F)+
        geom_point(data = x,
                   aes(x = UMAP_1, y = UMAP_2),color = 'white',size = .005)+
        scale_fill_viridis(option="A")+theme_classic()+NoAxes()+NoLegend()
    }
  }
}

# fig1E
{
  R_oe <- calTissueDist(object.meta,colname.cluster = "CellType_umap1",colname.tissue = "combine_label")
  or_data_plot <- R_oe[names(color_umap1),c('pre_R','post_R','pre_NR','post_NR')]
  
  bk <- c(seq(0,1,by=0.02),seq(1,2.2,by=0.02))
  colour_bk <- c(colorRampPalette(c("#2166ac","#d1e5f0"))(41),
                 colorRampPalette(c("#d1e5f0","#f7f7f7"))(10),
                 colorRampPalette(c("#f7f7f7","#fddbc7"))(10),
                 colorRampPalette(c("#fddbc7","#b2182b"))(51))
  
  anno_data <- object.meta[,.(count=.N),by=.(combine_label,CellType_umap1)]
  anno_data <- anno_data[,.(precent=count/sum(count),combine_label),by=.(CellType_umap1)] %>% spread(combine_label,precent) %>% data.frame()
  anno_data <- anno_data[rownames(or_data_plot),]
  row_anno = rowAnnotation('Sample\ncomposition' = anno_barplot(anno_data,border=F,bar_width = 0.9,width = unit(1, "cm"),
                                                                gp = gpar(fill = color_combine_fill, col = color_combine)),
                           annotation_name_gp = gpar(fontsize = 6),annotation_name_side = 'top')
  left_ann = rowAnnotation(Celltype = rownames(or_data_plot),col = list(Celltype = color_umap1))
  
  Heatmap(or_data_plot,
          cluster_rows=F,cluster_columns = F,
          col = circlize::colorRamp2(bk, colour_bk),
          row_names_side = 'left',
          row_names_gp = gpar(fontsize = 6),column_names_gp = gpar(fontsize = 6),
          right_annotation = row_anno,left_annotation = left_ann,
          border = T,border_gp = gpar(col = 'grey80', lwd = 0.5), 
          rect_gp  = gpar(col = 'white', lwd = 1), 
          width = unit(2, "cm"), height = unit(5, "cm")
  )
  
}

# fig1F
{
  sc_diff <- function(use_obj,savefile=file.path(fig_all,'celltyp_n_diff','T R vs NR diff.rda')){
    use_diff <- FindMarkers(use_obj,ident.1 = 'R',ident.2 = 'NR',group.by = 'response',test.use = "wilcox",min.pct = 0.05,logfc.threshold = 0)
    use_diff$gene <- rownames(use_diff)
    use_diff <- data.table(use_diff)
    setorder(use_diff,p_val_adj)
    save(use_diff,file = savefile)
    return(use_diff)
  }
  for (i in unique(object$CellType_n)) {
    use_obj <- subset(object,cells=all.metadata[CellType_n==i,cellid])
    use_diff <- sc_diff(use_obj,file.path('celltyp_n_diff',paste0(i,' R vs NR diff.rda')))
  }

  DEG_table_cell <- data.table()
  for (i in unique(object.meta$CellType_n)) {
    print(paste0(i,' Star...'))
    
    tmp_cell <- object.meta[CellType_n==i,cellid]
    
    load(file.path('celltyp_n_diff',paste0(i,' R vs NR diff.rda')))
    diff_gene <- use_diff[p_val_adj<0.05,gene]
    
    tmp_mat <- GetAssayData(object)[diff_gene,tmp_cell]
    tmp_index <- colSums(tmp_mat>0) 
    DEG_table_cell <- rbind(DEG_table_cell,data.table(cellid=names(tmp_index),DEG=tmp_index))
  }
  
  object$DEG <- DEG_table_cell[match(colnames(object),cellid)]$DEG
  FeaturePlot(object,features = 'DEG',pt.size = 1.5,raster = T)+
    scale_color_gradientn(colours = c(colorRampPalette(c("#fdf2b5","#F67B51"))(10),
                                      colorRampPalette(c("#F67B51","#A30023"))(10)))+NoAxes()
  
}

# fig1G, 1H, 1I
{
  ##T
  {
    load('TCGA_LIHC.rda')
    
    use_diff <- FindMarkers(object_T,only.pos=T,ident.1 = 'R',ident.2 = 'NR',group.by = 'response',
                            test.use = "wilcox",min.pct = 0.05,logfc.threshold = 1)
    use_diff$gene <- rownames(use_diff)
    use_diff <- data.table(use_diff)
    setorder(use_diff,p_val_adj)
    tmp <- use_diff[p_val_adj==0,] %>% arrange(-avg_log2FC)
    selectgene <- tmp[1:50,gene]
    use_gsva <- gsva(as.matrix(TPM),list(score=selectgene),method='ssgsea',kcdf='Gaussian',parallel.sz=8)
    use_gsva <- data.table(sampleid=colnames(use_gsva),use_score=use_gsva[1,])
    
    use_cli <- TCGA_LIHC_cli[sample_type=='Primary Tumor']
    use_cli$use_score <- use_gsva[match(use_cli$sampleID,sampleid),use_score]
    
    coxmodel <- coxph(Surv(PFI,PFI_event=='recurrence')~use_score,data=use_cli)
    smoothCoxph(use_cli$PFI,use_cli$PFI_event,use_cli$use_score,xlab="Response T cell like score",xlim=c(-1.5,1.5))
    cutoff <- surv_cutpoint(use_cli,time = 'PFI',event = 'PFI_event' ,variables =c('use_score'))$cutpoint$cutpoint
    
    myfit <- survfit(Surv(PFI/30,PFI_event=='recurrence')~use_score<cutoff,data=use_cli)
    ggsurvplot(fit=myfit,data=use_cli,pval=T,risk.table=F,
               ylab="Probability of PFS", xlab="Survival time (months)",
               palette = c('#F5A439','#0D9AD8'),
               linetype=c(1,1),legend=c(0.6,0.9),
               legend.labs = c("Response T cell like","Response T cell unlike"))
    
    plot_data <- use_diff
    plot_data$color_label <- 'no'
    plot_data[avg_log2FC<(-0.25)&p_val_adj<0.01,'color_label'] <- 'left'
    plot_data[avg_log2FC>0.25&p_val_adj<0.01,'color_label'] <- 'right'
    ggplot(plot_data,aes(avg_log2FC,-log10(p_val_adj)))+
      geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
      geom_vline(xintercept = log2(1.2), linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
      geom_vline(xintercept = -log2(1.2), linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
      geom_point(aes(color=color_label),size=0.5)+
      scale_color_manual(values = c('no'='#7F7D7D','left'='#4B38AD','right'='#E16047'))+
      ggrepel::geom_text_repel(aes(label = show_gene),color="red",
                               size=2,segment.size=0.02,segment.color='red',max.overlaps=10000)+
      theme_classic()+xlab("Log2(Fold Change))")+
      ylab("-Log10(adj.P.Value)")+ggtitle('T cells Response vs. Noresponse')+
      theme(plot.title = element_text(hjust=0.5),legend.position = 'none')
  }
}

# figS1A
FeaturePlot(object,features = c('KRT18','PTPRC','CD3D','KLRD1','CD79A','CD68','CSF3R','CD1C','IL3RA','MKI67'),pt.size = 2,ncol=5)

# figS1B
DimPlot(object, reduction = "umap",label =F,group.by = "TCR",pt.size = 1.5,raster=T,shuffle=T,cols=c("TCR+"="#d62728","TCR-" = "#C7C6C6"))+NoAxes()

# figS1C
{
  plot_data <- all.metadata[,.(count=.N),by=.(sample,CellType_umap1)] %>% arrange(sample)
  plot_data <- plot_data[,.(ratio=count/sum(count),allcount=sum(count),count,CellType_umap1),by=.(sample)]
  plot_data$CellType_umap1 <- factor(plot_data$CellType_umap1,levels = names(color_umap1))
  
  ggplot(plot_data, aes(x = sample, y = count, fill = CellType_umap1))+
    geom_bar(position = 'fill',stat = "identity",na.rm=T,width=0.8)+
    scale_y_continuous(labels = scales::percent)+
    scale_x_discrete(limits=sample_order)+
    scale_fill_manual(values = color_umap1)+
    theme_classic()+ggtitle('')+ylab('Cell proportion(%)')+xlab('')+
    theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size = 6,colour = 'black'),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.y = element_text(size = 6,colour = 'black'),
          axis.line = element_line(linewidth  = 0.235),
          axis.ticks = element_line(linewidth  = 0.235),
          plot.title = element_text(hjust=0.5,size = 8),
          legend.text = element_text(size = 6),
          legend.title = element_blank(),
          legend.key.size = unit(0.1, "inches")
    )
}

# figS1D
{
  table(all.metadata$CellType_2)
  use_data <- all.metadata
  use_data$use_celltype <- use_data$CellType_2
  use_data[use_celltype%in%c('pDC'),use_celltype:='Myeloid']
  
  plot_data <- use_data[,.(count=.N),by=.(combine_label,use_celltype)]
  plot_data <- plot_data[,.(sum_count=sum(count),ratio=count/sum(count),count,use_celltype),by=.(combine_label)]
  
  ggplot(plot_data,aes(combine_label,ratio,group = use_celltype,color=use_celltype))+
    geom_point(size=0.5)+geom_line(linewidth=0.1)+
    scale_color_manual(values = color_tmp)+
    scale_y_continuous(labels = scales::percent)+
    theme_classic()+ylab('Fraction')+xlab('')+
    theme(axis.text.x=element_text(angle=30,vjust=1,hjust=1,size = 6,colour = 'black'),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.y = element_text(size = 6,colour = 'black'),
          axis.line = element_line(linewidth  = 0.235), #保证坐标轴描边为0.5
          axis.ticks = element_line(linewidth  = 0.235),
          plot.title = element_text(hjust=0.5,size = 8),
          legend.text = element_text(size = 6),
          legend.title = element_blank(),
          panel.border = element_rect(linewidth  = 0.235,fill = NA)
    )
}




