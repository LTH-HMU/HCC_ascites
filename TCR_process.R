
dirs <- list.dirs(sc_tcr_path,full.names=F,recursive=F)

TCR.contigs <- data.table()
TCR.clonotypes <- data.table()
for (i in 1:length(dirs)) {
  print(i)
  tmp_tcr <- fread(file.path(sc_tcr_path,dirs[i],"filtered_contig_annotations.csv"))
  tmp_tcr$barcode <- paste0(dirs[i],'_',tmp_tcr$barcode)
  tmp_tcr$sample <- sample_names[i]
  TCR.contigs <- rbind(TCR.contigs,tmp_tcr)
  clonotype.i <- fread(file.path(sc_tcr_path,dirs[i],"clonotypes.csv"))
  clonotype.i$sample <- sample_names[i]
  TCR.clonotypes <- rbind(TCR.clonotypes,clonotype.i)
}

# rename barcode:
TCR.contigs$sample_cloneID <- paste0(TCR.contigs$sample,"_",TCR.contigs$raw_clonotype_id)
TCR.clonotypes$sample_cloneID <- paste0(TCR.clonotypes$sample,"_",TCR.clonotypes$clonotype_id)
TCR.clonotypes <- unique(TCR.clonotypes) %>% as.data.frame()
rownames(TCR.clonotypes) <- TCR.clonotypes$sample_cloneID
TCR.contigs$sample_fraction <- TCR.clonotypes[TCR.contigs$sample_cloneID,"proportion"]
TCR.contigs$sample_freq <- TCR.clonotypes[TCR.contigs$sample_cloneID,"frequency"]
TCR.contigs <- TCR.contigs[TCR.contigs$productive,]
save(TCR.contigs,TCR.clonotypes,file = file.path("productive filtered TCR.contigs and clonotypes for all samples.rda") )


#### Patient level clones ####
{
  #1. include samples to clonotypes:
  load(file = file.path(fig_tcr,"productive filtered TCR.contigs and clonotypes for all samples.rda"))
  sample2patient <- unique(all.metadata[,c("sample","patient")]) %>% data.frame()
  rownames(sample2patient) <- sample2patient$sample
  TCR.clonotypes$Patient <- sample2patient[TCR.clonotypes$sample,"patient"]
  TCR.contigs$Patient <- sample2patient[TCR.contigs$sample,"patient"]
  save(TCR.contigs,TCR.clonotypes,file = file.path("productive filtered TCR.contigs and clonotypes for all samples.rda") )
  
  # patient level clonotypes:
  {
    df <- data.table(TCR.clonotypes)
    df <- df[,.(Freq.Patient=sum(frequency)),by=.(cdr3s_aa,Patient)]
    TCR.clonotypes_patient <- NULL
    sampleID2PatientID <- NULL
    for (i in unique(df$Patient)) {
      temp <- df[df$Patient==i,]
      temp <- temp[order(temp$Freq.Patient,decreasing = T),]
      temp$Patient_cloneID <- paste0(temp$Patient,"-",seq_along(temp$cdr3s_aa))
      temp$Frac.Patient <- temp$Freq.Patient/sum(temp$Freq.Patient)
      temp <- data.frame(temp,row.names = temp$cdr3s_aa)
      x <- TCR.clonotypes[TCR.clonotypes$Patient==i,c("cdr3s_aa","sample_cloneID")]
      x$Patient_cloneID <- temp[x$cdr3s_aa,"Patient_cloneID"]
      sampleID2PatientID <- rbind(sampleID2PatientID,x)
      TCR.clonotypes_patient <- rbind(TCR.clonotypes_patient,temp)
    }
    TCR.clonotypes$Patient_cloneID <- sampleID2PatientID[TCR.clonotypes$sample_cloneID,"Patient_cloneID"]
    rownames(TCR.clonotypes_patient) <- TCR.clonotypes_patient$Patient_cloneID
  }
  
  cell2clone <- unique(TCR.contigs[,c(1,32:35)])
  cell2clone <- cbind(cell2clone,TCR.clonotypes[cell2clone$sample_cloneID,c('Patient','Patient_cloneID')])
  cell2clone <- cbind(cell2clone,TCR.clonotypes_patient[cell2clone$Patient_cloneID,c(3,5)])
  save(cell2clone,TCR.contigs,TCR.clonotypes,TCR.clonotypes_patient,file = file.path("productive filtered TCR.contigs and clonotypes for all samples.rda"))
}

