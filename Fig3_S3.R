
##Fig3A-B.
{
  curve <- estimateAbundance(TCR_meta,group = 'response',ci=0.95, nboot=100, clone="Patient_cloneID")
  plot(curve, colors = color_pdsd, legend_title="Sample")+
    scale_x_log10(breaks = c(1, 10, 100), labels = c("1", "10", "100"),limits=c(1,150))+
    ylab('Fraction of T cells')+
    xlab('Rank of top 150 clones')+
    ggtitle('All T cells')+
    theme(axis.text.x=element_text(size = 6,colour = 'black'),
          axis.title.y = element_text(size = 6),
          axis.title.x = element_text(size = 6),
          axis.text.y = element_text(size = 6,colour = 'black'),
          plot.title = element_text(hjust=0.5,size = 6)
    )
  
  isotype_test <- alphaDiversity(curve,group="response", min_q=0, max_q=2, step_q=1, nboot=200)
  plot(isotype_test, 2, colors=color_pdsd, main_title='',legend_title="")+
    ylab('TCR clonotype diversity \n (1/Simpson\'s index,mean±SD)')+
    scale_x_discrete(limits=names(color_pdsd))+
    theme(axis.text.x=element_text(size = 6,colour = 'black'),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.y = element_text(size = 6,colour = 'black'),
          axis.line = element_line(linewidth  = 0.235),
          axis.ticks = element_line(linewidth  = 0.235)
    )
}

#Fig3C.
{
  use_tcr <- TCR_meta
  clone_table <- use_tcr[,.(Freq.new=.N),by=.(Patient_cloneID)] %>% arrange(-Freq.new)
  use_tcr[Patient_cloneID%in%c(clone_table[1,Patient_cloneID]),use_label:='1']
  use_tcr[Patient_cloneID%in%c(clone_table[2:30,Patient_cloneID]),use_label:='2:30']
  use_tcr[Patient_cloneID%in%c(clone_table[31:100,Patient_cloneID]),use_label:='31:100']
  use_tcr[Patient_cloneID%in%c(clone_table[100:clone_table[,.N],Patient_cloneID]),use_label:='>100']
  use_tcr$use_label <- factor(use_tcr$use_label,levels = c('1','2:30','31:100','>100'))
  use_tcr$CellType_n <- factor(use_tcr$CellType_n,levels = names(color_T))
  
  use_tcr$use_celltype <- use_tcr$CellType_n
  use_tcr[CellType_2=='CD4+ cells',]$use_celltype <- 'CD4 T'
  use_tcr[CellType_2=='Other T',]$use_celltype <- 'Other T'
  plot_data <- use_tcr[,.(count=.N),by=.(use_celltype,use_label)]
  plot_data <- plot_data[,.(ratio=count/sum(count),allcount=sum(count),count,use_label),by=.(use_celltype)]
  
  ggplot(plot_data,aes(x=use_celltype,y=ratio,fill=use_label))+
    geom_bar(stat = 'identity',position = 'fill',width = 0.8)+
    scale_y_continuous(labels = scales::percent)+
    scale_fill_manual(values = c('#D04848','#F3B95F','#FDE767','#6895D2'))+
    theme_classic()+ggtitle('')+ylab('Cell proportion(%)')+xlab('')+labs(fill = "Clonal Indices")+
    theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1,size = 6,colour = 'black'),
          axis.title.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.y = element_text(size = 6,colour = 'black'),
          axis.line = element_line(linewidth  = 0.235),
          axis.ticks = element_line(linewidth  = 0.235),
          plot.title = element_text(hjust=0.5,size = 8),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          legend.key.size = unit(0.1, "inches")
    )
}

##Fig3D.
{
  use_data <- TCR_meta[CellType_2%in%c('CD8+ cells')]
  use_data <- left_join(use_data,AUC_res)
  plot_data <- use_data[,.(ISG_score = mean(ISG_select_score)),by=.(Patient_cloneID,response)]
  
  ggplot(plot_data,aes(response, ISG_score))+
    geom_boxplot(aes(fill=response,color=response),color='black',width=0.4,outlier.shape = NA,linewidth  = 0.235)+
    scale_fill_manual(values=color_pdsd_fill)+
    geom_jitter(fill='white',color='black',shape=21,width =0.1,size=0.5,stroke = 0.1)+
    ylab('Clonal ISG socre')+xlab('')+ggtitle('CD8+ clones')+
    stat_compare_means(comparisons=list(c('R','NR')),method ="wilcox.test",size=2)+
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

##Fig3E.
{
  use_data <- TCR_meta[CellType_2%in%c('CD8+ cells'),]
  pre_post_clone <- intersect(use_data[treat=='Pre',Patient_cloneID],use_data[treat=='Post',Patient_cloneID]) %>% unique()
  use_data <- use_data[Patient_cloneID%in%pre_post_clone,]
  use_data$use_label <- paste(use_data$Patient_cloneID,use_data$treat,sep = "_")
  
  use_mat <- GetAssayData(object_cd8)[,use_data$cellid]
  clone_mat <- groupMeans(use_mat, groups = use_data$use_label, na.rm = TRUE, sparse = T)
  clone_mat <- clone_mat[names(which(rowSums(clone_mat>0)>ncol(clone_mat)*0.25)),]
  plot_data_R <- clone_mat[,unique(use_data[response=='R',]$use_label)]
  plot_data_NR <- clone_mat[,unique(use_data[response=='NR',]$use_label)]

  #diff R
  {
    tmp_sample_pre <- paste(unique(use_data_R$Patient_cloneID),'Pre',sep = "_")
    tmp_sample_post <- paste(unique(use_data_R$Patient_cloneID),'Post',sep = "_")
    diff_result <- apply(plot_data_R,1,function(x){
      p <- wilcox.test(x[tmp_sample_post], x[tmp_sample_pre],paired=T)$p.value
      logfc <- mean(x[tmp_sample_post])-mean(x[tmp_sample_pre])
      mean1 <- mean(x[tmp_sample_post])
      mean2 <- mean(x[tmp_sample_pre])
      return(c('p'=p,'logFC'=logfc,'mean_post'=mean1,'mean_pre'=mean2))
    })
    diff_result <- t(diff_result)
    diff_result <- data.table(gene=rownames(diff_result),diff_result)
    diff_result$fdr <- p.adjust(diff_result$p,method = 'BH')
    setorder(diff_result,fdr)
    diff_result_R <- diff_result
  }
  
  #diff NR
  {
    tmp_sample_pre <- paste(unique(use_data_NR$Patient_cloneID),'Pre',sep = "_")
    tmp_sample_post <- paste(unique(use_data_NR$Patient_cloneID),'Post',sep = "_")
    diff_result <- apply(plot_data_NR,1,function(x){
      p <- wilcox.test(x[tmp_sample_post], x[tmp_sample_pre],paired=T)$p.value
      logfc <- mean(x[tmp_sample_post])-mean(x[tmp_sample_pre])
      mean1 <- mean(x[tmp_sample_post])
      mean2 <- mean(x[tmp_sample_pre])
      return(c('p'=p,'logFC'=logfc,'mean_post'=mean1,'mean_pre'=mean2))
    })
    diff_result <- t(diff_result)
    diff_result <- data.table(gene=rownames(diff_result),diff_result)
    diff_result$fdr <- p.adjust(diff_result$p,method = 'BH')
    setorder(diff_result,fdr)
    diff_result_NR <- diff_result
  }
  
  ##plot R
  {
    plot_data <- diff_result_R
    plot_data$color_label <- 'no'
    plot_data[logFC<0&fdr<0.01,'color_label'] <- 'left'
    plot_data[logFC>0&fdr<0.01,'color_label'] <- 'right'
    
    ggplot(plot_data,aes(logFC,-log10(fdr)))+
      geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
      geom_point(aes(color=color_label),size=1)+
      scale_color_manual(values = c('no'='#7F7D7D','left'='#7AB774','right'='#7f91c9'))+
      ggrepel::geom_text_repel(aes(label = show_gene),color="red",
                               size=2,segment.size=0.1,segment.color='red',max.overlaps=10000)+
      theme_classic()+xlab("Log2(Fold change)")+
      xlim(-1.3,1)+ylab("-Log10(adj.P.Value)")+ggtitle('Responders CD8+ clone-match Post vs. Pre')+
      theme(axis.text.x=element_text(colour = 'black'),
            axis.text.y = element_text(colour = 'black'),
            axis.line = element_line(linewidth  = 0.235),
            axis.ticks = element_line(linewidth  = 0.235),
            plot.title = element_text(hjust=0.5),
            legend.position = 'none'
      )
  }
  
  ##plot NR
  {
    plot_data <- diff_result_NR
    plot_data$color_label <- 'no'
    plot_data[logFC<0&fdr<0.01,'color_label'] <- 'left'
    plot_data[logFC>0&fdr<0.01,'color_label'] <- 'right'
    
    ggplot(plot_data,aes(logFC,-log10(fdr)))+
      geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "#999999",linewidth  = 0.235)+  
      geom_point(aes(color=color_label),size=1)+
      scale_color_manual(values = c('no'='#7F7D7D','left'='#7AB774','right'='#7f91c9'))+
      ggrepel::geom_text_repel(aes(label = show_gene),color="red",
                               size=2,segment.size=0.1,segment.color='red',max.overlaps=10000)+
      theme_classic()+xlab("Log2 (Fold Change)")+
      xlim(-0.5,0.7)+ylab("-Log10 (adj.P.Value)")+ggtitle('Non-responders CD8+ clone-match Post vs. Pre')+
      theme(axis.text.x=element_text(colour = 'black'),
            axis.text.y = element_text(colour = 'black'),
            axis.line = element_line(linewidth  = 0.235),
            axis.ticks = element_line(linewidth  = 0.235),
            plot.title = element_text(hjust=0.5),
            legend.position = 'none'
      )
  }
}

##Fig3F.
{
  use_data <- TCR_meta[CellType_2%in%c('CD8+ cells'),]
  pre_post_clone <- intersect(use_data[treat=='Pre',Patient_cloneID],use_data[treat=='Post',Patient_cloneID]) %>% unique()
  use_data <- use_data[Patient_cloneID%in%pre_post_clone,]
  use_data <- left_join(use_data,AUC_res)
  
  plot_data_R <- use_data[response=='R',.(NeoTCR_CD8=mean(NeoTCR_CD8)),by=.(Patient_cloneID,treat)]
  plot_data_R$NeoTCR_CD8 <- scale(plot_data_R$NeoTCR_CD8)
  
  plot_data_NR <- use_data[response=='NR',.(NeoTCR_CD8=mean(NeoTCR_CD8)),by=.(Patient_cloneID,treat)]
  plot_data_NR$NeoTCR_CD8 <- scale(plot_data_NR$NeoTCR_CD8)

  p1 <- ggplot(plot_data_R,aes(treat,NeoTCR_CD8))+
    geom_boxplot(alpha =0.5,linewidth  = 0.235,outlier.shape = NA, mapping = aes(fill = treat),color='black',width=0.5)+
    geom_point(size=1, shape=16,aes(group = Patient_cloneID,col = treat),alpha = 0.9,stroke = 0.2)+
    scale_color_manual(values = color_treat)+scale_fill_manual(values = color_treat)+
    geom_line(aes(group = Patient_cloneID), color = 'gray', lwd = 0.1)+
    stat_compare_means(comparisons=list(c('Pre','Post')),method ="wilcox.test",size=2.5,paired=T)+
    theme_classic()+ggtitle('Responders')+ylab('Clonal neoantigen reactivity')+xlab('')+
    theme(
      axis.text.x=element_text(size = 6,colour = 'black'),
      axis.title.y = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      axis.text.y = element_text(size = 6,colour = 'black'),
      axis.line = element_line(linewidth  = 0.235),
      axis.ticks = element_line(linewidth  = 0.235),
      plot.title = element_text(hjust=0.5,size = 8),
      legend.position = 'none')
  
  p2 <- ggplot(plot_data_NR,aes(treat,NeoTCR_CD8))+
    geom_boxplot(alpha =0.5,linewidth  = 0.235,outlier.shape = NA, mapping = aes(fill = treat),color='black',width=0.5)+
    geom_point(size=1, shape=16,aes(group = Patient_cloneID,col = treat),alpha = 0.9,stroke = 0.2)+
    scale_color_manual(values = color_treat)+
    scale_fill_manual(values = color_treat)+
    geom_line(aes(group = Patient_cloneID), color = 'gray', lwd = 0.1)+
    stat_compare_means(comparisons=list(c('Pre','Post')),method ="wilcox.test",size=2.5,paired=T)+
    theme_classic()+ggtitle('Non-responders')+ylab('Clonal neoantigen reactivity')+xlab('')+
    theme(
      axis.text.x=element_text(size = 6,colour = 'black'),
      axis.title.y = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      axis.text.y = element_text(size = 6,colour = 'black'),
      axis.line = element_line(linewidth  = 0.235),
      axis.ticks = element_line(linewidth  = 0.235),
      plot.title = element_text(hjust=0.5,size = 8),
      legend.position = 'none')
  p1+p2

}

##Fig3G.
{
  use_data <- TCR_meta[CellType_2%in%c('CD8+ cells')]
  use_data <- left_join(use_data,AUC_res)
  use_data2 <- use_data[,.(ISG_score = mean(ISG_select_score),NeoTCR_CD8=mean(NeoTCR_CD8)),by=.(Patient_cloneID,response)]
  use_data2$ISG_select_score_scale <- scale(use_data2$ISG_score)
  use_data2$ISG_label <- 'ISG_high'
  use_data2[ISG_select_score_scale<0,]$ISG_label <- 'ISG_low'
  
  ggplot(use_data2,aes(ISG_label, NeoTCR_CD8))+
    geom_boxplot(aes(fill=ISG_label),width=0.4,outlier.shape = NA,linewidth  = 0.235)+
    scale_fill_manual(values=color_ISG)+
    geom_jitter(fill='white',color='black',shape=21,width =0.1,size=0.5,stroke = 0.1)+
    ylab('clonal neoantigen reactivity')+xlab('')+ggtitle('CD8+ clones')+
    stat_compare_means(comparisons=list(c('ISG_high','ISG_low')),method ="wilcox.test",size=2)+
    theme_classic()+scale_x_discrete(limits=(names(color_ISG)))+
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

##Fig3H. FigS3B.
{
  use_data <- TCR_meta
  use_data <- left_join(use_data,AUC_res)
  use_data2 <- use_data[,.(ISG_score = mean(ISG_select_score)),by=.(Patient_cloneID,response)]
  use_data2$ISG_select_score_scale <- scale(use_data2$ISG_score)
  use_data2$ISG_label <- 'ISG_high'
  use_data2[ISG_select_score_scale<0,]$ISG_label <- 'ISG_low'
  
  use_data <- TCR_meta[Freq.Patient>1]
  plot_data_isghi <- use_data[Patient_cloneID%in%use_data2[ISG_label=='ISG_high',Patient_cloneID],]
  plot_data_isglow <- use_data[Patient_cloneID%in%use_data2[ISG_label=='ISG_low',Patient_cloneID],]
  
  #ISG hi
  {
    tmp_data1 <- plot_data_isghi
    tmp_data1 <- tmp_data1[,.(count=.N),by=.(CellType_n,treat)]
    tmp_data1 <- tmp_data1[,.(CellType_n,ratio=count/sum(count),count,allcount=sum(count)),by=.(treat)]
    tmp_data1$CellType_n <- factor(tmp_data1$CellType_n,levels = names(color_T))
    
    p1 <- ggplot(tmp_data1, aes(x = treat, y = ratio, fill = CellType_n,stratum = CellType_n,alluvium = CellType_n))+
      geom_flow(width = 0.3,curve_type = "sine",alpha = 1,size = 0.1)+
      geom_col(width = 0.2,color=NA,size=0.1)+
      scale_y_continuous(labels = scales::percent)+
      scale_fill_manual(values = color_T)+
      theme_classic()+ggtitle('ISG high T cell clones')+ylab('% of clone expanded T cells')+xlab('')+
      theme(axis.text.x=element_text(size = 6,colour = 'black'),
            axis.title.y = element_text(size = 8),
            axis.title.x = element_text(size = 8),
            axis.text.y = element_text(size = 6,colour = 'black'),
            axis.line = element_line(linewidth  = 0.235),
            axis.ticks = element_line(linewidth  = 0.235),
            plot.title = element_text(hjust=0.5,size = 8),
            legend.text = element_text(size = 6),
            legend.title = element_blank(),
            legend.key.size = unit(0.1, "inches"),
            strip.text=element_text(size=8),
            strip.background=element_blank()
      )
  }
  
  #ISG low
  {
    tmp_data1 <- plot_data_isglow
    tmp_data1 <- tmp_data1[,.(count=.N),by=.(CellType_n,treat)]
    tmp_data1 <- tmp_data1[,.(CellType_n,ratio=count/sum(count),count,allcount=sum(count)),by=.(treat)]
    tmp_data1$CellType_n <- factor(tmp_data1$CellType_n,levels = names(color_T))
    
    p2 <- ggplot(tmp_data1, aes(x = treat, y = ratio, fill = CellType_n,stratum = CellType_n,alluvium = CellType_n))+
      geom_flow(width = 0.3,curve_type = "sine",alpha = 1,size = 0.1)+
      geom_col(width = 0.2,color=NA,size=0.1)+
      scale_y_continuous(labels = scales::percent)+
      scale_fill_manual(values = color_T)+
      theme_classic()+ggtitle('ISG low T cell clones')+ylab('% of clone expanded T cells')+xlab('')+
      theme(axis.text.x=element_text(size = 6,colour = 'black'),
            axis.title.y = element_text(size = 8),
            axis.title.x = element_text(size = 8),
            axis.text.y = element_text(size = 6,colour = 'black'),
            axis.line = element_line(linewidth  = 0.235), 
            axis.ticks = element_line(linewidth  = 0.235),
            plot.title = element_text(hjust=0.5,size = 8),
            legend.text = element_text(size = 6),
            legend.title = element_blank(),
            legend.key.size = unit(0.1, "inches"),
            strip.text=element_text(size=8),
            strip.background=element_blank()
      )
  }
  
  ##FigS3B R vs NR
  {
    ##R
    {
      tmp_data1 <- use_data[response=='R',]
      tmp_data1 <- tmp_data1[,.(count=.N),by=.(CellType_n,treat)]
      tmp_data1 <- tmp_data1[,.(CellType_n,ratio=count/sum(count),count,allcount=sum(count)),by=.(treat)]
      tmp_data1$CellType_n <- factor(tmp_data1$CellType_n,levels = names(color_T))
      
      p1 <- ggplot(tmp_data1, aes(x = treat, y = ratio, fill = CellType_n,stratum = CellType_n,alluvium = CellType_n))+
        geom_flow(width = 0.3,curve_type = "sine",alpha = 1,size = 0.1)+
        geom_col(width = 0.2,color=NA,size=0.1)+
        scale_y_continuous(labels = scales::percent)+
        scale_fill_manual(values = color_T)+
        theme_classic()+ggtitle('R')+ylab('% of clone expanded T cells')+xlab('')+
        theme(axis.text.x=element_text(size = 6,colour = 'black'),
              axis.title.y = element_text(size = 8),
              axis.title.x = element_text(size = 8),
              axis.text.y = element_text(size = 6,colour = 'black'),
              axis.line = element_line(linewidth  = 0.235),
              axis.ticks = element_line(linewidth  = 0.235),
              plot.title = element_text(hjust=0.5,size = 8),
              legend.text = element_text(size = 6),
              legend.title = element_blank(),
              legend.key.size = unit(0.1, "inches"),
              strip.text=element_text(size=8),
              strip.background=element_blank()
        )
    }
    ##NR
    {
      tmp_data2 <- use_data[response=='NR',]
      tmp_data2 <- tmp_data2[,.(count=.N),by=.(CellType_n,treat)]
      tmp_data2 <- tmp_data2[,.(CellType_n,ratio=count/sum(count),count,allcount=sum(count)),by=.(treat)]
      tmp_data2$CellType_n <- factor(tmp_data2$CellType_n,levels = names(color_T))
      
      p2 <- ggplot(tmp_data2, aes(x = treat, y = ratio, fill = CellType_n,stratum = CellType_n,alluvium = CellType_n))+
        geom_flow(width = 0.3,curve_type = "sine",alpha = 1,size = 0.1)+
        geom_col(width = 0.2,color=NA,size=0.1)+
        scale_y_continuous(labels = scales::percent)+
        scale_fill_manual(values = color_T)+
        theme_classic()+ggtitle('NR')+ylab('% of clone expanded T cells')+xlab('')+
        theme(axis.text.x=element_text(size = 6,colour = 'black'),
              axis.title.y = element_text(size = 8),
              axis.title.x = element_text(size = 8),
              axis.text.y = element_text(size = 6,colour = 'black'),
              axis.line = element_line(linewidth  = 0.235),
              axis.ticks = element_line(linewidth  = 0.235),
              plot.title = element_text(hjust=0.5,size = 8),
              legend.text = element_text(size = 6),
              legend.title = element_blank(),
              legend.key.size = unit(0.1, "inches"),
              strip.text=element_text(size=8),
              strip.background=element_blank()
        )
    }
  }
}

##Fig3I.
{
  use_data <- TCR_meta
  use_data <- left_join(use_data,AUC_res)
  use_data$ISG_select_score_scale <- scale(use_data$ISG_select_score)
  use_data$ISG_label <- 'ISG_high'
  use_data[ISG_select_score_scale<0,]$ISG_label <- 'ISG_low'
  selsect_type <- c("CD4_Mem","CD4_GZMK","CD4_EFF","CD4_Tfh",
                    "CD8_Tcm","CD8_GZMK","CD8_CTL",'NKT','MAIT')
  use_data <- use_data[CellType_n%in%selsect_type]
  
  clone_overlap_treat <- function(use_data){
    df <- use_data[,.(n=.N),by=.(Patient_cloneID)]
    use_data <- use_data[Patient_cloneID%in%df[n>1,Patient_cloneID],]
    use_data$clone_id <- use_data$Patient_cloneID
    mat <- matrix(data = NA,nrow = length(selsect_type),ncol = length(selsect_type),dimnames=list(selsect_type,selsect_type))
    for (i in 1:length(selsect_type)) {
      for (j in 1:length(selsect_type)) {
        x=unique(use_data[CellType_n%in%selsect_type[i]&treat=='Pre',]$clone_id)
        y=unique(use_data[CellType_n%in%selsect_type[j]&treat=='Post',]$clone_id)
        iOrd <- intersect(x,y)
        mat[i,j] <- length(iOrd)/length(y)
      }
    }
    return(mat)
  }
  
  plot_data_sighi_mat <- clone_overlap_treat(use_data[ISG_label=='ISG_high'])
  plot_data_siglow_mat <- clone_overlap_treat(use_data[ISG_label=='ISG_low'])
  bk <- c(seq(0,0.7,by=0.01))
  colour_bk <- c(colorRampPalette(c("#3288bd","#abdda4"))(17),
                 colorRampPalette(c("#abdda4","#fee08b"))(17),
                 colorRampPalette(c("#fee08b","#fdae61"))(17),
                 colorRampPalette(c("#fdae61","#9e0142"))(20))
  
  h1 <- Heatmap(plot_data_sighi_mat,
                cluster_rows=F,
                cluster_columns = F,
                col = circlize::colorRamp2(bk, colour_bk),
                row_names_side = 'left',
                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 6),
                row_title_gp = gpar(fontsize = 8),
                column_title_gp = gpar(fontsize = 8),
                border = T,
                border_gp = gpar(col = 'grey80', lwd = 0.5),
                rect_gp  = gpar(col = 'white', lwd = 0.1),
                width = unit(2.5, "cm"), height = unit(2.5, "cm"),
                column_title = "ISG high T cells",row_title = "Post-treatment"
  )
  h2 <- Heatmap(plot_data_siglow_mat,
                cluster_rows=F,
                cluster_columns = F,
                col = circlize::colorRamp2(bk, colour_bk),
                row_names_side = 'left',
                row_names_gp = gpar(fontsize = 6),
                column_names_gp = gpar(fontsize = 6),
                row_title_gp = gpar(fontsize = 8),
                column_title_gp = gpar(fontsize = 8),
                border = T,
                border_gp = gpar(col = 'grey80', lwd = 0.5),
                rect_gp  = gpar(col = 'white', lwd = 0.1),
                width = unit(2.5, "cm"), height = unit(2.5, "cm"),
                column_title = "ISG low T cells",row_title = "Post-treatment"
  )
  h1+h2
}

##FigS3A.
{
  TCR_meta$cdr3s_aa <- TCR.clonotypes_patient[TCR_meta$Patient_cloneID,]$cdr3s_aa
  use_data <- TCR_meta
  patient1 <- unique(use_data$patient)
  patient2 <- gsub('MA','PT',patient1)
  
  clone_overlap_patient <- function(use_data){
    mat <- matrix(data = NA,nrow = length(patient2),ncol = length(patient2),dimnames=list(patient2,patient2))
    mat_num <- matrix(data = NA,nrow = length(patient2),ncol = length(patient2),dimnames=list(patient2,patient2))
    for (i in 1:length(patient2)) {
      for (j in 1:length(patient2)) {
        x=unique(use_data$cdr3s_aa[use_data$patient==patient1[i]])
        y=unique(use_data$cdr3s_aa[use_data$patient==patient1[j]])
        iOrd <- intersect(x,y)
        mat[i,j] <- length(iOrd)/length(y)
        mat_num[i,j] <- length(iOrd)
      }
    }
    return(list(mat,mat_num))
  }
  
  res_list <- clone_overlap_patient(use_data)
  
  h1 <- Heatmap(res_list[[1]],
                cluster_rows=F,cluster_columns = F,
                col = circlize::colorRamp2(c(0,0.5,1), c('#3288BD','#74c0fc','#e7f5ff')),
                row_names_side = 'left',row_names_gp = gpar(fontsize = 6),column_names_gp = gpar(fontsize = 6),
                row_title_gp = gpar(fontsize = 8),column_title_gp = gpar(fontsize = 8),
                cell_fun = function(j, i, x, y, width, height, fill){
                grid.text(sprintf("%d", res_list[[2]][i, j]),x,y,gp = gpar(fontsize = 6))},
                border = T,border_gp = gpar(col = 'grey80', lwd = 0.5), 
                rect_gp  = gpar(col = 'white', lwd = 0.1),  
                width = unit(4, "cm"), height = unit(4, "cm"),
                column_title = "Patients",row_title = "")
}

##FigS3C.
{
  use_diff <- FindMarkers(object_cd4T,ident.1 = 'R',ident.2 = 'NR',group.by = 'response',test.use = "wilcox",min.pct = 0.2,return.thresh = 0.01,logfc.threshold = 0.25)
  use_diff$gene <- rownames(use_diff)
  use_diff <- data.table(use_diff)
  setorder(use_diff,p_val_adj)
  
  bp1 <- enrichGO(use_diff[p_val_adj<0.05&avg_log2FC>0,gene],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
  bp2 <- enrichGO(use_diff[p_val_adj<0.05&avg_log2FC<0,gene],keyType="SYMBOL",ont="BP",OrgDb="org.Hs.eg.db")
  use_go_1 <- bp1@result[bp1@result$p.adjust < 0.05,] %>% data.table()
  use_go_2 <- bp2@result[bp2@result$p.adjust < 0.05,] %>% data.table()
  up_select <- c('I-kappaB kinase/NF-kappaB signaling',
                 'positive regulation of cytokine production',
                 'T cell differentiation',
                 'interferon-beta production',
                 'response to type I interferon',
                 'T-helper cell differentiation',
                 "CD4-positive, alpha-beta T cell activation")
  down_select <- c('cytolysis',
                   'lymphocyte mediated immunity',
                   'antigen processing and presentation of peptide antigen',
                   'regulation of leukocyte apoptotic process')
  plot_data <- rbind(data.table(use_go_1[Description%in%up_select,],label='R'),
                     data.table(use_go_2[Description%in%down_select,],label='NR'))
  plot_data$q <- -log10(plot_data$qvalue)
  plot_data$Description <- factor(plot_data$Description,levels = rev(plot_data$Description))
  plot_data[label=='NR','q'] <- -(plot_data[label=='NR',q])
  up <- plot_data[label=='R',]
  down <- plot_data[label=='NR',]
  
  ggplot(plot_data,aes(x =q, y = Description, fill = label)) + #数据映射
    geom_col(alpha=0.6)+theme_bw()+
    scale_fill_manual(values = color_pdsd)+
    geom_text(data = up,aes(x = -0.2, y = Description, label = Description),size = 2,hjust = 1)+ 
    geom_text(data = down,aes(x = 0.2, y = Description, label = Description),size = 2,hjust = 0)+ 
    xlim(-7.2,7.2)+
    labs(x = '-log10(qvalue)', y = ' ', title = 'CD4T R vs. NR ')
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

##FigS3D.
{
  plot_data <- left_join(object_T.meta,AUC_res)
  ggplot(plot_data,aes(CellType_2,ISG_select_score,fill=response,color=response))+
    scale_fill_manual(values = color_pdsd_fill)+scale_color_manual(values = color_pdsd)+
    geom_boxplot(width=0.7,linewidth=0.235,outlier.shape = NA)+
    ylab('ISG score')+xlab('')+theme_classic()+
    stat_compare_means(aes(group = response,label = paste0(..p.format..)),method = "wilcox.test",hide.ns = F,size=2)+
    theme(
      axis.text.x=element_text(angle=45,vjust=1,hjust=1,size = 6,colour = 'black'),
      axis.title.y = element_text(size = 6),
      axis.title.x = element_text(size = 6),
      axis.text.y = element_text(size = 6,colour = 'black'),
      axis.line = element_line(linewidth  = 0.235),
      axis.ticks = element_line(linewidth  = 0.235),
      plot.title = element_text(hjust=0.5,size = 6),
      legend.position = 'none'
    )
}












