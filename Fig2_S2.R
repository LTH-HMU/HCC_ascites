


##Fig2A.
{
  all.metadata$CD8A <- GetAssayData(object)['CD8A',]
  all.metadata[CD8A>0,'gene_label'] <- 'CD8A+'
  all.metadata$CD4 <- GetAssayData(object)['CD4',]
  all.metadata[CD4>0,'gene_label'] <- 'CD4+'
  
  plot_data <- all.metadata[!CellType_2=='Epithelial',.(count=.N),by=.(sample,response,gene_label)] %>% arrange(sample)
  plot_data <- plot_data[,.(ratio=count/sum(count),allcount=sum(count),count,gene_label,response),by=.(sample)]
  
  ggplot(plot_data[gene_label=='CD8A+'],aes(response, ratio))+
    geom_boxplot(aes(fill=response,color=response),width=0.4,outlier.shape = NA,linewidth  = 0.235)+
    scale_fill_manual(values=color_pdsd_fill)+
    scale_color_manual(values=color_pdsd)+
    geom_jitter(fill='white',color='black',shape=21,width =0.15,size=1,stroke = 0.2)+
    ylab('% of immune cells')+xlab('')+ggtitle('CD8A+ cells')+
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
  
  ggplot(plot_data[gene_label=='CD4+'],aes(response, ratio))+
    geom_boxplot(aes(fill=response,color=response),width=0.4,outlier.shape = NA,linewidth  = 0.235)+
    scale_fill_manual(values=color_pdsd_fill)+
    scale_color_manual(values=color_pdsd)+
    geom_jitter(fill='white',color='black',shape=21,width =0.15,size=1,stroke = 0.2)+
    ylab('% of immune cells')+xlab('')+ggtitle('CD4+ cells')+
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

##Fig2B.
DimPlot(object_T,cols=color_T,group.by = 'CellType_n',label =T,label.size = 3,pt.size = 2,raster=T,shuffle=T)+NoAxes()

##Fig2C.
{
  DEG_table_cell <- data.table()
  for (i in unique(object_T.meta$CellType_n)) {
    tmp_cell <- object_T.meta[CellType_n==i,cellid]
    load(file.path('celltyp_n_diff',paste0(i,' R vs NR diff.rda')))
    diff_gene <- use_diff[p_val_adj<0.05,gene]
    tmp_mat <- GetAssayData(object_T)[diff_gene,tmp_cell]
    tmp_index <- colSums(tmp_mat>0) 
    DEG_table_cell <- rbind(DEG_table_cell,data.table(cellid=names(tmp_index),DEG=tmp_index))
  }
  object_T$DEG <- DEG_table_cell[match(colnames(object_T),cellid)]$DEG
  FeaturePlot(object_T,features = 'DEG',pt.size = 1.5,raster = T)+
    scale_color_gradientn(colours = c(colorRampPalette(c("#fdf2b5","#F67B51"))(10),
                                      colorRampPalette(c("#F67B51","#A30023"))(10)))+NoAxes()
}

##Fig2D.
{
  use_data <- left_join(object_T.meta[,'cellid'],TCR_meta)
  use_data[is.na(Patient_cloneID),]$Freq.Patient <- 0
  object_T$Frequency_plot <- log10(use_data$Freq.Patient+1)
  
  FeaturePlot(object_T,cells=object_T.meta[response=='R',cellid],features = 'Frequency_plot',pt.size = 1.2,raster = F)+
    scale_color_gradientn(colours = c('#ebeee8',colorRampPalette(c("#313695","#abd9e9"))(10),
                                      colorRampPalette(c("#abd9e9","#fee090"))(10),
                                      colorRampPalette(c("#fee090","#a50026"))(10)),
                          limits=c(0,3.5))+NoAxes()+ggtitle('clone size in Responders')
  
  FeaturePlot(object_T,cells=object_T.meta[response=='NR',cellid],features = 'Frequency_plot',pt.size = 1.2,raster = F)+
    scale_color_gradientn(colours = c('#ebeee8',colorRampPalette(c("#313695","#abd9e9"))(10),
                                      colorRampPalette(c("#abd9e9","#fee090"))(10),
                                      colorRampPalette(c("#fee090","#a50026"))(10)),
                          limits=c(0,3.4))+NoAxes()+ggtitle('clone size in NonResponders')
  
}

##Fig2E.
{
  use_diff <- FindMarkers(object_cd8,ident.1 = 'R',ident.2 = 'NR',group.by = 'response',test.use = "wilcox",min.pct = 0.2,logfc.threshold = 0.25)
  use_diff$gene <- rownames(use_diff)
  use_diff <- data.table(use_diff)
  
  bp1 <- enrichGO(use_diff[avg_log2FC>0.25&p_val_adj<0.05,gene],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
  bp2 <- enrichGO(use_diff[avg_log2FC<(-0.25)&p_val_adj<0.05,gene],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
  use_go_1 <- bp1@result[bp1@result$p.adjust < 0.05,] %>% data.table()
  use_go_2 <- bp2@result[bp2@result$p.adjust < 0.05,] %>% data.table()
  
  up_select <- c('response to type I interferon',
                 'T cell differentiation',
                 'positive regulation of cytokine production',
                 'I-kappaB kinase/NF-kappaB signaling',
                 'type I interferon production',
                 'type I interferon signaling pathway',
                 "response to interferon-gamma")
  down_select <- c('MHC class II protein complex assembly',
                   'regulation of leukocyte cell-cell adhesion',
                   'lymphocyte proliferation',
                   'negative regulation of immune system process')
  
  plot_data <- rbind(data.table(use_go_1[Description%in%up_select,],label='R'),
                     data.table(use_go_2[Description%in%down_select,],label='NR'))
  plot_data$q <- -log10(plot_data$qvalue)
  
  plot_data$Description <- factor(plot_data$Description,levels = rev(plot_data$Description))
  plot_data[label=='NR','q'] <- -(plot_data[label=='NR',q])
  
  up <- plot_data[label=='R',]
  down <- plot_data[label=='NR',]
  
  ggplot(plot_data,aes(x =q, y = Description, fill = label))
    geom_col(alpha=0.6)
    theme_bw()+
    scale_fill_manual(values = color_pdsd)+
    geom_text(data = up,aes(x = -0.2, y = Description, label = Description),size = 2,hjust = 1)
    geom_text(data = down,aes(x = 0.2, y = Description, label = Description),size = 2,hjust = 0)
    labs(x = '-log10(qvalue)', y = ' ', title = 'CD8T R vs NR Enriched Pathway')
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

##Fig2F.
{
  plot_data <- object_T.meta[CellType_2=='CD8+ cells',]
  plot_data$CD8A <- GetAssayData(object)['CD8A',plot_data$cellid]
  plot_data$IFIT3 <- GetAssayData(object)['IFIT3',plot_data$cellid]
  
  ggplot(plot_data[response=='R'],aes(x=CD8A,y=IFIT3))+
    geom_point(color='#88ABCB',size=0.5)+
    geom_hline(yintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    geom_vline(xintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    theme_classic()+xlim(0,4)+ylim(0,5)+
    ggtitle('')+ylab(paste0('IFIT3+'))+xlab('CD8A+')+
    theme(axis.text.x=element_text(size = 6,colour = 'black'),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.y = element_text(size = 6,colour = 'black'),
          axis.line = element_line(linewidth  = 0.235),
          axis.ticks = element_line(linewidth  = 0.235),
          plot.title = element_text(hjust=0.5,size = 8)
    )
  
  ggplot(plot_data[response=='NR'],aes(x=CD8A,y=IFIT3))+
    geom_point(color='#F29EA1',size=0.5)+
    geom_hline(yintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    geom_vline(xintercept = 0.15, linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
    theme_classic()+xlim(0,4)+ylim(0,5)+
    ggtitle('')+ylab(paste0('IFIT3+'))+xlab('CD8A+')+
    theme(axis.text.x=element_text(size = 6,colour = 'black'),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.y = element_text(size = 6,colour = 'black'),
          axis.line = element_line(linewidth  = 0.235),
          axis.ticks = element_line(linewidth  = 0.235),
          plot.title = element_text(hjust=0.5,size = 8)
    )
}

##Fig2G.
{
  use_data <- TCR_meta
  use_data$IFIT3 <- GetAssayData(object_T)['IFIT3',use_data$cellid]
  use_data$IFIT3_label <- ifelse(use_data$IFIT3>0,'IFIT3+','IFIT3-')
  
  tmp_plot_pie <- function(plotdata=use_data[response=='R'&Freq.Patient==1,],title='1'){
    ggpie(data = plotdata, group_key = 'IFIT3_label', count_type = "full",
          border_size=0.2,border_color = NA,
          label_info = "all", label_type = "none",
          label_size = 3, label_pos = "in",label_split = NULL)+
      ggtitle(title)+
      theme( plot.title = element_text(hjust=0.5,size = 8),
             legend.text = element_text(size = 6),
             legend.title = element_text(size = 6),
             legend.key.height = unit(0.25,"cm"),
             legend.key.width = unit(0.25,"cm"),
             legend.position = "right")+
      scale_fill_manual(values=c('#EBECED','#EC7B73'),name='Celltype')
  }
  
  p1 <- tmp_plot_pie(plotdata=use_data[response=='R'&Freq.Patient==1,],title='1')
  p2 <- tmp_plot_pie(plotdata=use_data[response=='R'&Freq.Patient>1&Freq.Patient<=9,],title='2~9')
  p3 <- tmp_plot_pie(plotdata=use_data[response=='R'&Freq.Patient>=10&Freq.Patient<=30,],title='10~30')
  p4 <- tmp_plot_pie(plotdata=use_data[response=='R'&Freq.Patient>30,],title='n>30')
  p5 <- tmp_plot_pie(plotdata=use_data[response=='NR'&Freq.Patient==1,],title='1')
  p6 <- tmp_plot_pie(plotdata=use_data[response=='NR'&Freq.Patient>1&Freq.Patient<=9,],title='2~9')
  p7 <- tmp_plot_pie(plotdata=use_data[response=='NR'&Freq.Patient>=10&Freq.Patient<=30,],title='10~30')
  p8 <- tmp_plot_pie(plotdata=use_data[response=='NR'&Freq.Patient>30,],title='n>30')
  (p1|p2|p3|p4)/(p5|p6|p7|p8)
}

#Fig2H-I.
{
  CD8_CTL <- subset(object_cd8_NKT,cells=object_cd8_NKT_meta[CellType_n%in%c('CD8_CTL'),cellid])
  CD8_CTL <- NormalizeData(CD8_CTL)
  CD8_CTL <- FindVariableFeatures(CD8_CTL,nfeatures=2000)
  
  TCR.genes <- grep("^TR[ABGD][VJ]",rownames(CD8_CTL),value = T)# keep GD TCR genes for GDT cells
  BCR.genes <- c(grep("^IG[KHL][VJC]",rownames(CD8_CTL),value = T),
                 grep("^AC[0-9]",rownames(CD8_CTL),value = T))# some RNA genes are also excluded.
  MT.genes <- c(grep("^MT-|MTRN",rownames(CD8_CTL),value = T))
  RP.genes <- c(grep("^RP[SL]",rownames(CD8_CTL),value = T))
  HSP.genes <- c(grep("^HSP|DNAJ",rownames(CD8_CTL),value = T))
  dissociation <- c("ATF3","BTG2","CEBPB","CEBPB-AS1","CEBPD","CXCL1","EGR1","FOS","FOSB",
                    "FOSL1","FOSL1P1","FOSL2","ID3","IER2","JUN","JUNB","JUND","MT1A",
                    "MT1B","MT1E","MT1F","MT1G","MT1H","MT1L","MT1M","MT1X","MT2A",
                    "NFKBIA","NR4A1","PPP1R15A","SOCS3","UBC","ZFP36")
  var.genes <- VariableFeatures(CD8_CTL)
  var.genes <- setdiff(var.genes,c(TCR.genes,BCR.genes,MT.genes,RP.genes,HSP.genes))
  VariableFeatures(CD8_CTL) <- var.genes
  
  expression_data <- CD8_CTL@assays$RNA@data
  cell_metadata <- CD8_CTL@meta.data
  gene_annotation <- data.frame(gene_short_name=rownames(expression_data),
                                stringsAsFactors = F,row.names = rownames(expression_data))
  pd <- new("AnnotatedDataFrame", data = cell_metadata)
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  object.cds <- newCellDataSet(cellData=expression_data,phenoData = pd,
                               featureData = fd,expressionFamily=negbinomial.size())
  object.cds <- estimateSizeFactors(object.cds)
  object.cds <- estimateDispersions(object.cds)
  object.cds <- detectGenes(object.cds, min_expr = 0.1)
  
  ordering_genes <- intersect(VariableFeatures(CD8_CTL),allgene)[1:500]
  object.cds <- setOrderingFilter(object.cds, ordering_genes)
  object.cds <- reduceDimension(object.cds, max_components = 2,reduction_method = 'DDRTree',pseudo_expr=1)
  object.cds <- orderCells(object.cds)

  # metadata
  {
    object.cds_meta <- data.table(cellid=colnames(object.cds),pData(object.cds))
    object.cds_meta$combine_label <- object_T.meta[match(object.cds_meta$cellid,cellid),combine_label]
    object.cds_meta$response <- object_T.meta[match(object.cds_meta$cellid,cellid),response]
    df <- data.frame(t(object.cds@reducedDimS))
    object.cds_meta$x <- df[,1]
    object.cds_meta$y <- df[,2]
    
    object.cds_meta <- left_join(object.cds_meta,AUC_res)
    object.cds$ISG_score <- tmp$ISG_score
    object.cds$cytotoxicity_score <- tmp$cytotoxicity_score
    object.cds$neoantigen_reactivity_cd8_score <- tmp$NeoTCR_CD8
  }
  # base plot
  {
    #response
    {
      plot_cell_trajectory(object.cds, color_by = "response",cell_size=2,show_branch_points = F,show_tree = F)+
        scale_color_manual(values = color_pdsd)+
        theme(
          axis.text = element_text(size = 12,colour = 'black'),
          axis.title = element_text(size = 12),
          legend.position = c(0.5,0.8),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.key.width = unit(0.5, "cm"),
          legend.key.height = unit(0.5, "cm"),
          axis.line = element_line(linewidth  = 0.235), 
          axis.ticks = element_line(linewidth  = 0.235)
        )+guides(color = guide_legend(override.aes = list(size = 2)))
    }
    
    #Pseudotime
    {
      plot_cell_trajectory(object.cds, color_by = "Pseudotime",cell_size=2,show_branch_points = F)+
        scale_color_gradientn(colours = c(colorRampPalette(c("#363d99","#a0cae1"))(10),
                                          colorRampPalette(c("#a0cae1","#d1dcbf"))(5),
                                          colorRampPalette(c("#d1dcbf","#fedf90"))(5),
                                          colorRampPalette(c("#fedf90","#a8092a"))(10)))+
        theme(
          axis.text = element_text(size = 12,colour = 'black'),
          axis.title = element_text(size = 12),
          legend.position = c(0.5,0.8),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.key.width = unit(0.5, "cm"),
          legend.key.height = unit(0.5, "cm"),
          axis.line = element_line(linewidth  = 0.235), 
          axis.ticks = element_line(linewidth  = 0.235)
        )+ guides(color = guide_colorbar(direction = "horizontal"))
    }
    
    
    ##nFeature_RNA
    {
      object.cds$nFeature_RNA_plot <- object.cds_meta$nFeature_RNA_plot
      plot_cell_trajectory(object.cds, color_by = "nFeature_RNA_plot",cell_size=2,show_branch_points = F,show_tree = F)+
        scale_color_gradientn(colours = c(colorRampPalette(c("seagreen","#ccdd8b"))(7),
                                          colorRampPalette(c("#ccdd8b","#fdf2b5"))(5),
                                          colorRampPalette(c("#fdf2b5","#882601"))(10)))+
        theme(
          axis.text = element_text(size = 12,colour = 'black'),
          axis.title = element_text(size = 12),
          legend.position = c(0.5,0.8),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.key.width = unit(0.5, "cm"), 
          legend.key.height = unit(0.5, "cm"),
          axis.line = element_line(linewidth  = 0.235), 
          axis.ticks = element_line(linewidth  = 0.235)
        )+ guides(color = guide_colorbar(direction = "horizontal")) 
    }
    
    ##ISG_score
    {
      ggplot(object.cds_meta,aes(x=x,y=y,color=ISG_score))+
        geom_point(size=2)+
        scale_color_gradientn(colours = c("#f0f0f0",colorRampPalette(c("#313695","#abd9e9"))(5),
                                          colorRampPalette(c("#abd9e9","#fee090"))(5),
                                          colorRampPalette(c("#fee090","#a50026"))(10)))+
        theme_classic()+NoAxes()+ggtitle('ISG score')+
        theme(
          plot.title = element_text(hjust=0.5,size = 8),
          axis.text = element_text(size = 12,colour = 'black'),
          axis.title = element_text(size = 12),
          legend.position = c(0.5,0.8),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.key.width = unit(0.5, "cm"), 
          legend.key.height = unit(0.5, "cm"),
          axis.line = element_line(linewidth  = 0.235), 
          axis.ticks = element_line(linewidth  = 0.235)
        )+ guides(color = guide_colorbar(direction = "horizontal"))
    }
    
    ##cytotoxicity
    {
      ggplot(object.cds_meta,aes(x=x,y=y,color=cytotoxicity_score))+
        geom_point(size=2)+
        scale_color_gradientn(
          breaks = c(0,0.3, 0.6, 0.9),                    
          labels = c("0", "0.3", "0.6", "0.9"),          
          colours = c("#f0f0f0",colorRampPalette(c("#313695","#abd9e9"))(14),
                      colorRampPalette(c("#abd9e9","#fee090"))(6),
                      colorRampPalette(c("#fee090","#a50026"))(6)))+
        theme_classic()+NoAxes()+ggtitle('Cytotoxicity score')+
        theme(
          plot.title = element_text(hjust=0.5,size = 8),
          axis.text = element_text(size = 12,colour = 'black'),
          axis.title = element_text(size = 12),
          legend.position = c(0.5,0.8),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 6),
          legend.key.width = unit(0.5, "cm"), 
          legend.key.height = unit(0.5, "cm"),
          axis.line = element_line(linewidth  = 0.235), 
          axis.ticks = element_line(linewidth  = 0.235)
        )+ guides(color = guide_colorbar(direction = "horizontal")) 
    }
    
    ##neoantigen_reactivity
    {
      ggplot(object.cds_meta,aes(x=x,y=y,color=neoantigen_reactivity_cd8_score))+
        geom_point(size=2)+
        scale_color_gradientn(colours = c("#f0f0f0",colorRampPalette(c("#313695","#abd9e9"))(10),
                                          colorRampPalette(c("#abd9e9","#fee090"))(10),
                                          colorRampPalette(c("#fee090","#a50026"))(20)))+
        theme_classic()+NoAxes()+ggtitle('Neoantigen reactivity')+
        theme(
          plot.title = element_text(hjust=0.5,size = 8),
          axis.text = element_text(size = 12,colour = 'black'),
          axis.title = element_text(size = 12),
          legend.position = c(0.5,0.8),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.key.width = unit(0.5, "cm"),
          legend.key.height = unit(0.5, "cm"),
          axis.line = element_line(linewidth  = 0.235), 
          axis.ticks = element_line(linewidth  = 0.235)
        )+ guides(color = guide_colorbar(direction = "horizontal"))
    }
    
    #gene
    {
      gene_plot <- function(plot_data,use_gene){
        plot_data[[use_gene]] <- GetAssayData(object_T)[use_gene,plot_data$cellid]
        p1 <- ggplot(data=plot_data,mapping = aes(x=Pseudotime, y=get(use_gene))) + 
          geom_smooth(method = 'loess',color='#DCB129', se=T,linewidth  = 0.5,span = 1)+
          theme_light() +
          xlab('Pseudotime') + ylab(use_gene)+ggtitle('')+
          theme(axis.text.x=element_text(size = 6,colour = 'black'),
                axis.title.y = element_text(size = 6),
                axis.title.x = element_text(size = 6),
                axis.text.y = element_text(size = 6,colour = 'black'),
                axis.line = element_line(linewidth  = 0.235),
                axis.ticks = element_line(linewidth  = 0.235),
                plot.title = element_text(hjust=0.5,size = 8)
          )
        return(p1)
      }
      
      p1 <- gene_plot(object.cds_meta,'IFIT3')
      p2 <- gene_plot(object.cds_meta,'ISG15')
      p3 <- gene_plot(object.cds_meta,'OASL')
      p4 <- gene_plot(object.cds_meta,'LAG3')
      p1/p2/p3/p4
      
      p1 <- gene_plot(object.cds_meta,'GZMB')
      p2 <- gene_plot(object.cds_meta,'IFNG')
      p3 <- gene_plot(object.cds_meta,'TNF')
      p4 <- gene_plot(object.cds_meta,'HLA-B')
      p1/p2/p3/p4
    }
  }
}

##FigS2A.
DimPlot(object_cd8, cols=color_T,group.by = 'CellType_n',label =T,label.size = 3,pt.size = 2,raster=T,shuffle=T)+NoAxes()

##FigS2B.
DimPlot(object_cd4, cols=color_T,group.by = 'CellType_n',label =T,label.size = 3,pt.size = 2,raster=T,shuffle=T)+NoAxes()

##FigS2C.
{
  markers <- list(  `T cell`=c(c('CD4',"CD8A"),
                               c("CCR7",'SELL','LEF1','TCF7'),
                               c('CD44','IL7R','LTB','GZMK'),
                               c('TOX',"ICOS","PDCD1",'TOX2','LAG3'),
                               c("IL2RA",'FOXP3','IKZF2'),
                               c('GZMB','GZMH',"GNLY","NKG7",'CST7','PRF1'),
                               c("KIR2DL3","KIR3DL2","FCGR3A"),
                               c('TRDV2','TRGV9'),
                               c("RORC",'SLC4A10','TRAV1-2','CCR6'),
                               c('MKI67',"TOP2A",'TYMS')))
  DotPlot(object_T, features = markers,group.by = "CellType_n",dot.scale=3)+
    scale_y_discrete(limits=rev(names(color_T)))+theme_bw()+
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

#FigS2D.
{
  use_gene <- c("IFIT2","IFIT3","IFIT5","ISG15","IRF1","IRF7","IRF9","OASL","ADAR","EIF2AK2","EPSTI1","MX1","OGFR","SP110","STAT1")
  VlnPlot(object_cd8,features = use_gene,cols=color_pdsd,
          group.by="response",fill.by="ident",
          stack = T,flip=T,adjust=2)+
    ylab("Expression")+xlab("")+ggtitle("CD8T")+
    theme(title=element_text(size = 6),
          axis.text.x=element_text(size = 6,colour = 'black',angle = 0,hjust = 0.5),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.y = element_text(size = 6,colour = 'black'),
          axis.line = element_line(linewidth  = 0.235),
          axis.ticks = element_line(linewidth  = 0.235),
          strip.text=element_text(size=6)
    )+NoLegend()+
    stat_signif(comparisons = list(c("R","NR")),textsize=2,y_position=0.5,map_signif_level=T,tip_length=0)
}

##FigS2E.
{
  use_data <- TCR_meta[response=='R',]
  use_data$IFIT3 <- GetAssayData(object_T)['IFIT3',use_data$cellid]
  use_data$use_label <- NULL
  use_data[Freq.Patient==1,use_label:='1']
  use_data[Freq.Patient>1&Freq.Patient<=9,use_label:='2~9']
  use_data[Freq.Patient>9&Freq.Patient<=30,use_label:='10~30']
  use_data[Freq.Patient>30,use_label:='>30']
  use_data$use_label <- factor(use_data$use_label,levels = c('1','2~9','10~30','>30'))

  ggplot(use_data,aes(use_label, IFIT3,fill=use_label))+
    geom_boxplot(width=0.6,outlier.shape = NA,linewidth  = 0.235)+
    scale_fill_manual(values=c('#D9D9D9','#3288BD','#ABDDA4','#D53E4F'))+
    ylab('Clonal IFIT3 exp.')+
    xlab('Clone size')+ggtitle('T cell clones in Responders')+
    stat_compare_means(comparisons=list(c('1','>30'),c('2~9','>30'),c('10~30','>30')),method ="wilcox.test",size=2)+
    theme_classic()+
    theme(
      axis.text.x=element_text(size = 6,colour = 'black'),
      axis.title.y = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      axis.text.y = element_text(size = 6,colour = 'black'),
      axis.line = element_line(linewidth  = 0.235),
      axis.ticks = element_line(linewidth  = 0.235),
      plot.title = element_text(hjust=0.5,size = 8)
    )
}

##FigS2F.
{
  use_meta <- all.metadata[CellType_n!='Epithelial',]
  plot_data <- use_meta[,.(count=.N),by=.(sample,response,CellType_n)] %>% arrange(sample)
  plot_data <- plot_data[,.(ratio=count/sum(count),allcount=sum(count),count,CellType_n,response),by=.(sample)]
  plot_data <- plot_data[CellType_n%in%c('CD8_GZMK','CD8_CTL','NKT'),]
  
  ggplot(plot_data,aes(response, ratio))+
    geom_boxplot(aes(fill=response,color=response),width=0.4,outlier.shape = NA,linewidth  = 0.235)+
    scale_fill_manual(values=color_pdsd_fill)+
    scale_color_manual(values=color_pdsd)+
    geom_jitter(fill='white',color='black',shape=21,width =0.1,size=1,stroke = 0.2)+
    ylab('% of immune cells')+xlab('')+ggtitle('')+
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
      legend.position = 'none',
      strip.text=element_text(size=6),
      strip.background=element_blank()
    )+facet_wrap(~CellType_n,scales = 'free_y',nrow = 1)
}
