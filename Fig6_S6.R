
##BayesPrism
{
  ##scData
  {
    DC_cell <- all.metadata[CellType_umap1%in%c('DC','pDC'),cellid]
    Mo_Mac_cell <- sample(all.metadata[CellType_umap1%in%c('Monocyte','Macrophage'),cellid],5000)
    CD4T_cell <- sample(all.metadata[CellType_2%in%c('CD4+ cells'),cellid],5000)
    CD8T_cell <- sample(all.metadata[CellType_2%in%c('CD8+ cells'),cellid],5000)
    NK_cell <- sample(all.metadata[celltype_main%in%c('NK'),cellid],5000)
    B_cell <- all.metadata[celltype_main%in%c('B Lineage'),cellid]
    
    use_obj <- subset(object, cells = c(DC_cell,Mo_Mac_cell,CD4T_cell,CD8T_cell,B_cell,NK_cell))
    use_obj$use_celltype <- use_obj$celltype_main
    use_obj@meta.data[use_obj$CellType_2=='CD4+ cells',]$use_celltype <- 'CD4T'
    use_obj@meta.data[use_obj$CellType_2=='CD8+ cells',]$use_celltype <- 'CD8T'
    use_obj@meta.data[use_obj$CellType_umap1%in%c('DC','pDC'),]$use_celltype <- 'DC'
    use_obj@meta.data[use_obj$CellType_umap1%in%c('Monocyte','Macrophage'),]$use_celltype <- 'Mono_Macro'
    
    sc.dat <- GetAssayData(use_obj,slot ="count") %>% as.matrix() %>% t()
    cell.type.labels <- use_obj$use_celltype
    cell.state.labels <- use_obj$CellType_n
    plot.cor.phi (input=sc.dat,input.labels=cell.type.labels,title="cell type correlation",cexRow=0.5, cexCol=0.5)
    #Filter outlier genes
    sc.stat <- plot.scRNA.outlier(input=sc.dat,cell.type.labels=cell.type.labels,species="hs",return.raw=TRUE )
    sc.dat.filtered <- cleanup.genes (input=sc.dat,input.type="count.matrix",
                                      species="hs",gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY"),
                                      exp.cells=5)
  }
  
  ##bulk data GO30140&IMbrave150 cohorts
  {
    load(file = 'GO30140_IMbrave150.rda')
    use_cli <- EGA_cli[drug!='Sorafenib'&treat=='Pre-treatment',]
    bk.dat <- t(EGA_count[,use_cli$sampleid])
  }
  
  ##bulk data TCGA-LIHC
  {
    load(file.path('TCGA_LIHC.rda'))
    use_cli <- TCGA_LIHC_cli[sample_type=='Primary Tumor',]
    bk.dat <- t(count[,use_cli$sampleID])
  }
  
  ##bulk data ICGC
  {
    load(file.path('ICGC_LIRI.rda'))
    use_cli <- ICGC_LIRI_cli
    bk.dat <- t(ICGC_LIRI_mat_count[,use_cli$icgc_specimen_id])
  }
  
  ## process
  {
    sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,gene.type = "protein_coding")
    myPrism <- new.prism(
      reference=sc.dat.filtered.pc,
      mixture=bk.dat,
      input.type="count.matrix",
      cell.type.labels = cell.type.labels,
      cell.state.labels = cell.state.labels,
      key=NULL,
      outlier.cut=0.01,
      outlier.fraction=0.1,
    )
    bp.res <- run.prism(prism = myPrism, n.cores=20)
  }
}

##Fig6A-C; FigS6D TCGA_LIHC
{
  load(file.path('TCGA_LIHC.rda'))
  use_gene <- c('IFIT1','IFIT2','IFIT3','IFIT5','IFI35','ISG15','IRF1','IRF3','IRF7',
                'IRF9','OASL','ADAR','EIF2AK2','EPSTI1','GBP4','MX1','OGFR','SP110','STAT1','STAT2')
  use_gsva <- gsva(as.matrix(TPM),list(score=use_gene),method='ssgsea',kcdf='Gaussian',parallel.sz=8)
  use_gsva <- data.table(sampleid=colnames(use_gsva),use_score=use_gsva[1,])
  use_cli <- TCGA_LIHC_cli[sample_type=='Primary Tumor']
  use_cli$use_score <- use_gsva[match(use_cli$sampleID,sampleid),use_score]
  
  ##FigS6D
  {
    cutoff <- surv_cutpoint(use_cli,time = 'PFI',event = 'PFI_event' ,variables =c('use_score'))$cutpoint$cutpoint
    myfit <- survfit(Surv(PFI/30,PFI_event=='recurrence')~use_score<cutoff,data=use_cli)
    ggsurvplot(fit=myfit,data=use_cli,pval=T,risk.table=F,
                     ylab="Probability of PFS", xlab="Survival time (months)",
                     palette = c('#F5A439','#0D9AD8'),
                     linetype=c(1,1),
                     legend=c(0.7,0.9),
                     legend.labs = c("ISG high","ISG low"))+ggtitle('TCGA LIHC dataset (n=370)')
    #heatmap
    {
      setorder(use_cli,-use_score)
      use_cli$gourp <- ifelse(use_cli$use_score>cutoff,'ISG_high','ISG_low')
      use_mat <- TPM[sig_ISG_select,use_cli$sampleID]
      plot_data <- t(scale(t(use_mat)))
      bk <- c(seq(-2,-0.1,by=0.02),seq(0,2,by=0.02))
      colour_bk <- c(colorRampPalette(c("#2166ac","#d1e5f0"))(83),
                     colorRampPalette(c("#d1e5f0","#f7f7f7"))(15),
                     colorRampPalette(c("#f7f7f7","#fddbc7"))(15),
                     colorRampPalette(c("#fddbc7","#b2182b"))(84))
      col_bar = HeatmapAnnotation(
        ssgsea_ISG= anno_barplot(use_cli$use_score, bar_width = 1, baseline =min(use_cli$use_score),border =F,
                                 gp = gpar(fill = "#cad7c5", col = NA),height=unit(1, "cm")),
        sample_label = use_cli$gourp,col = list('sample_label'=c('high'='#EFBF58','low'="#316AA7"))
      )
      Heatmap(plot_data,show_row_names = T,show_column_names = F,
              cluster_columns = F,cluster_rows = F,
              clustering_method_rows='ward.D',
              clustering_method_columns='ward.D',
              row_names_side = 'right',
              col = circlize::colorRamp2(bk, colour_bk),
              column_names_gp = gpar(fontsize = 6),
              row_names_gp = gpar(fontsize = 6),
              top_annotation = col_bar,
              border = T,border_gp = gpar(col = 'grey80', lwd = 0.5),
              width = unit(8, "cm"), height = unit(4, "cm")
      )
    }
  }
  
  #Fig6A
  {
    ciber_res <- fread(file.path(fig_6,'CIBERSORTx_TCGA-LIHC_abs.txt'))
    plot_data <- ciber_res
    plot_data$M1_M2 <- plot_data$`Macrophages M1`/plot_data$`Macrophages M2`
    plot_data <- melt(plot_data,id.vars=1, measure.vars=2:ncol(plot_data),variable.name="immuCell",value.name="Proportion")
    plot_data <- left_join(plot_data,use_cli[,c('sampleid','ISG_label')])
    plot_data$ISG_label <- factor(plot_data$ISG_label,levels = c('ISG_high','ISG_low'))

    ggplot(plot_data[immuCell=='T cells CD8'],aes(ISG_label,Proportion))+
      geom_boxplot(aes(fill=ISG_label),width=0.5,color='black',outlier.shape = NA,linewidth  = 0.235)+
      geom_jitter(fill='white',color='black',shape=21,width =0.15,size=0.5,stroke = 0.1)+
      scale_fill_manual(values=color_ISG)+
      stat_compare_means(comparisons=list(c('ISG_high','ISG_low')),method ="wilcox.test",size=2)+
      theme_classic()+ggtitle('CD8+ T cells')+ylab('cell proportion')+xlab('')+
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
    
    ggplot(plot_data[immuCell=='M1_M2'],aes(ISG_label,Proportion))+
      geom_boxplot(aes(fill=ISG_label),width=0.5,color='black',outlier.shape = NA,linewidth  = 0.235)+
      geom_jitter(fill='white',color='black',shape=21,width =0.15,size=0.5,stroke = 0.1)+
      scale_fill_manual(values=color_ISG)+
      stat_compare_means(comparisons=list(c('ISG_high','ISG_low')),method ="wilcox.test",size=2,label.y=1.9)+
      theme_classic()+ggtitle('Macrophages')+ylab('Fraction of M1 vs. M2 ')+xlab('')+
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
  
  #Fig6B-C
  {
    #The results of the BayesPrism analysis for the TCGA_LIHC dataset.
    load(file = file.path(fig_6,"TCGA_LIHC_BayesPrism_imm.rda"))
    theta <- get.fraction (bp=bp.res,which.theta="final",state.or.type="type")
    
    ##Fig6B ISG score
    {
      tmp_res <- data.table()
      for (i in colnames(theta)) {
        Z_mat <- get.exp (bp=bp.res,state.or.type="type",cell.name=i)
        bp_data <- vst(round(t(Z_mat)) )
        selectgene <- sig_ISG_select
        use_gsva <- gsva(as.matrix(bp_data),list(score=selectgene),method='ssgsea',kcdf='Gaussian',parallel.sz=8)
        use_gsva_new <- data.table(sampleID=colnames(use_gsva),use_score=use_gsva[1,],celltype=i)
        tmp_res <- rbind(tmp_res,use_gsva_new)
      }
      TCGA_LIHC_BayesPrism_imm_ISG <- left_join(tmp_res,use_cli[,c('sampleID','ISG_label')])
      
      ggplot(TCGA_LIHC_BayesPrism_imm_ISG,aes(celltype,use_score))+
        geom_boxplot(aes(fill=ISG_label),width=0.8,color='black',outlier.shape = NA,linewidth  = 0.235,)+
        scale_fill_manual(values=color_ISG)+
        stat_compare_means(aes(group = ISG_label,label = paste0(..p.format..)),method = "wilcox.test",hide.ns = F,size=1.8)+
        scale_x_discrete(limits=c('CD4T',"CD8T","NK","B Lineage","DC","Mono_Macro"))+
        theme_classic()+ggtitle('TCGA_LIHC BayesPrism')+ylab('ISG score')+xlab('')+
        theme(
          axis.text.x=element_text(angle=45,vjust=1,hjust=1,size = 6,colour = 'black'),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.y = element_text(size = 6,colour = 'black'),
          axis.line = element_line(linewidth  = 0.235),
          axis.ticks = element_line(linewidth  = 0.235),
          plot.title = element_text(hjust=0.5,size = 8),
          legend.position = 'none'
        )
    }
    
    ##Fig6C
    {
      #T cells Cytotoxicity
      {
        Z_mat <- get.exp (bp=bp.res,state.or.type="type",cell.name='CD8T')
        bp_data <- vst(round(t(Z_mat)) )
        selectgene <- c("GZMA","GZMB","GZMH","GZMM","GZMK","GNLY","PRF1","CTSW")
        use_gsva <- gsva(as.matrix(bp_data),list(score=selectgene),method='ssgsea',kcdf='Gaussian',parallel.sz=8)
        use_gsva <- data.table(sampleID=colnames(use_gsva),use_score=use_gsva[1,])
        use_cli$fun_socre <- use_gsva$use_score
        
        ggplot(use_cli,aes(ISG_label,fun_socre))+
          geom_boxplot(aes(fill=ISG_label),width=0.5,color='black',outlier.shape = NA,linewidth  = 0.235)+
          geom_jitter(fill='white',color='black',shape=21,width =0.15,size=0.5,stroke = 0.1)+
          scale_fill_manual(values=color_ISG)+
          stat_compare_means(comparisons=list(c('ISG_high','ISG_low')),method ="wilcox.test",size=2)+
          theme_classic()+ggtitle('CD8+ T cell')+ylab('Cytotoxicity score')+xlab('')+
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
      #NK cells Cytotoxicity
      {
        Z_mat <- get.exp (bp=bp.res,state.or.type="type",cell.name='NK')
        bp_data <- vst(round(t(Z_mat)) )
        selectgene <- c("GZMA","GZMB","GZMH","GZMM","GZMK","GNLY","PRF1","CTSW")
        use_gsva <- gsva(as.matrix(bp_data),list(score=selectgene),method='ssgsea',kcdf='Gaussian',parallel.sz=8)
        use_gsva <- data.table(sampleID=colnames(use_gsva),use_score=use_gsva[1,])
        use_cli$fun_socre <- use_gsva$use_score
        
        ggplot(use_cli,aes(ISG_label,fun_socre))+
          geom_boxplot(aes(fill=ISG_label),width=0.5,color='black',outlier.shape = NA,linewidth  = 0.235)+
          geom_jitter(fill='white',color='black',shape=21,width =0.15,size=0.5,stroke = 0.1)+
          scale_fill_manual(values=color_ISG)+
          stat_compare_means(comparisons=list(c('ISG_high','ISG_low')),method ="wilcox.test",size=2)+
          theme_classic()+ggtitle('NK cells')+ylab('Cytotoxicity')+xlab('')+
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
      
      #B cell activation
      {
        Z_mat <- get.exp (bp=bp.res,state.or.type="type",cell.name='B Lineage')
        bp_data <- vst(round(t(Z_mat)) )
        selectgene <- fread(file = file.path('msigDB_GO_B_activation.txt'))
        use_gsva <- gsva(as.matrix(bp_data),list(score=selectgene),method='ssgsea',kcdf='Gaussian',parallel.sz=8)
        use_gsva <- data.table(sampleID=colnames(use_gsva),use_score=use_gsva[1,])
        use_cli$fun_socre <- use_gsva$use_score
        
        ggplot(use_cli,aes(ISG_label,fun_socre))+
          geom_boxplot(aes(fill=ISG_label),width=0.5,color='black',outlier.shape = NA,linewidth  = 0.235)+
          geom_jitter(fill='white',color='black',shape=21,width =0.15,size=0.5,stroke = 0.1)+
          scale_fill_manual(values=color_ISG)+
          stat_compare_means(comparisons=list(c('ISG_high','ISG_low')),method ="wilcox.test",size=2)+
          theme_classic()+ggtitle('B cells')+ylab('B cell activation (GO:0042113)')+xlab('')+
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
      
      #DCs Antigen processing
      {
        Z_mat <- get.exp (bp=bp.res,state.or.type="type",cell.name='DC')
        bp_data <- vst(round(t(Z_mat)) )
        selectgene <- msigdb_go_bp[gs_name=='GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION',gene_symbol]
        use_gsva <- gsva(as.matrix(bp_data),list(score=selectgene),method='ssgsea',kcdf='Gaussian',parallel.sz=8)
        use_gsva <- data.table(sampleID=colnames(use_gsva),use_score=use_gsva[1,])
        use_cli$fun_socre <- use_gsva$use_score
        
        ggplot(use_cli,aes(ISG_label,fun_socre))+
          geom_boxplot(aes(fill=ISG_label),width=0.5,color='black',outlier.shape = NA,linewidth  = 0.235)+
          geom_jitter(fill='white',color='black',shape=21,width =0.15,size=0.5,stroke = 0.1)+
          scale_fill_manual(values=color_ISG)+
          stat_compare_means(comparisons=list(c('ISG_high','ISG_low')),method ="wilcox.test",size=2)+
          theme_classic()+ggtitle('DCs')+ylab('Antigen processing')+xlab('')+
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

    
  }
}

##Fig6D
{
  files <- list.files(TCGA_pancancer)
  cox_TCGA_res <- data.table()
  for (i in files) {
    tcga_name <- gsub('.rda','',i)
    load(file = file.path(tcga_path,i))
    use_tpm <- merge$tpm
    use_cli <- data.table(merge$clinical)
    use_cli <- use_cli[sample_type=='Primary Tumor']
    use_tpm <- use_tpm[,use_cli$barcode]
    use_gsva <- gsva(as.matrix(use_tpm),list(score=sig_ISG_select),method='ssgsea',kcdf='Gaussian',parallel.sz=8)
    use_gsva <- data.table(sampleid=colnames(use_gsva),use_score=use_gsva[1,])
    use_cli$use_score <- use_gsva$use_score
    use_cli$use_score <- scale(use_cli$use_score)[,1]
    coxmodel <- coxph(Surv(PFI,PFI_event=='1')~use_score,data=use_cli)
    tmp_dat <- data.table(
      `Dataset`= tcga_name,
      HR=signif(summary(coxmodel)$conf.int[1],3),
      `95%CI`=paste0(signif(summary(coxmodel)$conf.int[3],3),"-",signif(summary(coxmodel)$conf.int[4],3)),
      `p_value`=signif(summary(coxmodel)$coefficients[5],3),
    )
    cox_TCGA_res <- rbind(cox_TCGA_res,tmp_dat)
  }
  setorder(cox_TCGA_res,HR)
  
  # plot
  {
    plot_cox_res <- cox_TCGA_res
    plot_cox_res <- plot_cox_res[,c('Dataset','p_value','HR','low','upper')]
    plot_cox_res$` ` <- paste(rep(" ", dim(plot_cox_res)[1]), collapse = " ")
    
    tm <- forest_theme(base_size = 10,
                       ci_pch = 15,ci_col = "blue4",ci_fill = "blue4",
                       ci_alpha = 0.8,ci_lty = 1,ci_lwd = 1.5,ci_Theight = 0.2,
                       refline_lwd = 1,refline_lty = "dashed",refline_col = "grey20",
                       vertline_lwd = 1,vertline_lty = "dashed",vertline_col = "grey20",
                       footnote_cex = 0.6,footnote_fontface = "italic",footnote_col = "red4")
    p <- forest(data = as.matrix(plot_cox_res[,c(1:2,6)]), 
                est = plot_cox_res[, HR] %>% as.numeric(),
                lower = plot_cox_res[, low] %>% as.numeric(),
                upper = plot_cox_res[, upper] %>% as.numeric(),
                sizes = 0.4,ci_column = 3,ref_line = 0,
                arrow_lab = c("Low risk", "High Risk"),theme = tm) 
    pp <- edit_plot(p, row = c(1,2), gp = gpar(col = "red"))
    pp <- add_border(pp, part = "header", where = 'bottom')
    pp <- add_border(pp, part = "header", where = 'top')
  }
}

##Fig6E
{
  tmp_res <- data.table()
  for (i in colnames(theta)) {
    Z_mat <- get.exp (bp=bp.res,state.or.type="type",cell.name=i)
    bp_data <- vst(round(t(Z_mat)))
    use_gene <- c('IFIT1','IFIT2','IFIT3','IFIT5','IFI35','ISG15','IRF1','IRF3','IRF7',
                  'IRF9','OASL','ADAR','EIF2AK2','EPSTI1','GBP4','MX1','OGFR','SP110','STAT1','STAT2')
    use_gsva <- gsva(as.matrix(bp_data),list(score=use_gene),method='ssgsea',kcdf='Gaussian',parallel.sz=8)
    use_gsva_new <- data.table(sampleid=colnames(use_gsva),use_score=use_gsva[1,],celltype=i)
    tmp_res <- rbind(tmp_res,use_gsva_new)
  }
  BayesPrism_imm_ISG <- left_join(tmp_res,use_cli[,c('sampleid','response')])
  BayesPrism_imm_ISG$response <- factor(BayesPrism_imm_ISG$response,levels = c('CR','PR','SD','PD'))
  BayesPrism_imm_ISG$celltype <- factor(BayesPrism_imm_ISG$celltype,levels = c('CD4T',"CD8T","NK","B Lineage","DC","Mono_Macro"))
  
  ggplot(BayesPrism_imm_ISG,aes(response,use_score))+
    geom_boxplot(aes(fill=response),width=0.8,color='black',outlier.shape = NA,linewidth  = 0.235,)+
    geom_jitter(fill='white',color='black',shape=21,width =0.2,size=0.5,stroke = 0.1)+
    scale_fill_manual(values=c("#3D86B1","#A7CBDF",'#6FB56B','#9267A9'))+
    stat_compare_means(label = 'p.format',method = "anova",hide.ns = F,size=2)+
    theme_classic()+ggtitle('GO30140 BayesPrism')+ylab('ISG score')+xlab('')+
    theme(
      axis.text.x=element_text(size = 6,colour = 'black'),
      axis.title.y = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      axis.text.y = element_text(size = 6,colour = 'black'),
      axis.line = element_line(linewidth  = 0.235),
      axis.ticks = element_line(linewidth  = 0.235),
      plot.title = element_text(hjust=0.5,size = 8),
      strip.text=element_text(size=8),
      strip.background=element_blank(),
      legend.position = 'none'
    )+facet_wrap(~celltype,ncol = 6)
}