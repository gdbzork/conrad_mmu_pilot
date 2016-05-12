pairs_in_single_fastq <- function(fn) {
  pilot_log$info("Counting fastq reads in %s...",fn)
  cmd <- sprintf("gunzip -c %s | wc -l",fn)
  s <- system(cmd,intern=TRUE)
  pairs <- as.numeric(s) / 4 # because 4 lines per Fastq read
  return(pairs)
}

pairs_in_fastq <- function(fns) {
  pairs <- bplapply(fns,pairs_in_single_fastq)
  return(as.numeric(as.vector(pairs)))
}

reads_in_single_bam <- function(fn) {
  pilot_log$info("Counting bam reads in %s...",fn)
  cmd <- sprintf("samtools flagstat %s",fn)
  s <- system(cmd,intern=TRUE)
  res <- strsplit(s," ")
  reads <- as.numeric(c(res[[5]][1],res[[9]][1]))
  return(reads)
}

reads_in_bam <- function(fns) {
  pilot_log$info("wd: %s",getwd())
  res <- bplapply(fns,reads_in_single_bam)
  return(Reduce(rbind,res))
}

plot_heatmap <- function(data) {
  sample_dist <- dist(t(assay(data)))
  sample_mat <- as.matrix(sample_dist)
  rownames(sample_mat) <- paste(data$SampleID,
                                data$Treatment,
                                data$Age,
                                data$Sex,
                                sep="_")
  colnames(sample_mat) <- NULL
  colours <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
  pheatmap(sample_mat,
           clustering_distance_rows=sample_dist,
           clustering_distance_cols=sample_dist,
           col=colours)
}

de_analysis <- function(data,model,subset_vec=NULL,alpha=0.05,contrast=NULL) {
  de_obj <- DESeqDataSet(data,design = model)
  if (!is.null(subset_vec)) {
    de_obj <- subset(de_obj,select=subset_vec)
  }
  de_diff <- DESeq(de_obj)
  if (is.null(contrast)) {
    de_results <- results(de_diff,alpha=0.05)
  } else {
    de_results <- results(de_diff,alpha=0.05,contrast=contrast)
  }
  return(de_results)
}

venn2 <- function(data1,data2,name1,name2,title) {
  universe <- union(data1,data2)
  in1 <- universe %in% data1
  in2 <- universe %in% data2
  mat <- cbind(in1,in2)
  colnames(mat) <- c(name1,name2)
  cmat <- vennCounts(mat)
  vennDiagram(cmat,main=title)
}

vennUpDown <- function(res1,res2,name1,name2,title,alpha) {
  res1up <- rownames(res1)[!is.na(res1$padj)&res1$padj<alpha&res1$log2FoldChange>0]
  res1down <- rownames(res1)[!is.na(res1$padj)&res1$padj<alpha&res1$log2FoldChange<0]
  res2up <- rownames(res2)[!is.na(res2$padj)&res2$padj<alpha&res2$log2FoldChange>0]
  res2down <- rownames(res2)[!is.na(res2$padj)&res2$padj<alpha&res2$log2FoldChange<0]
  title_up <- paste0(title," UP",collapse=" ")
  venn2(res1up,res2up,name1,name2,title_up)
  title_down <- paste0(title," DOWN",collapse=" ")
  venn2(res1down,res2down,name1,name2,title_down)
}

mkDataFrame <- function(data,fn) {
  entrez <- rownames(data)
  res <- data.frame(entrezID=rownames(data),
                    symbol=mapIds(org.Mm.eg.db,keys=entrez,column='SYMBOL',keytype='ENTREZID',multiVals='first'),
                    name=mapIds(org.Mm.eg.db,keys=entrez,column='GENENAME',keytype='ENTREZID',multiVals='first'),
                    ensembl=mapIds(org.Mm.eg.db,keys=entrez,column='ENSEMBL',keytype='ENTREZID',multiVals='first'),
                    refseq=mapIds(org.Mm.eg.db,keys=entrez,column='REFSEQ',keytype='ENTREZID',multiVals='first'),
                    logFC=data$log2FoldChange,
                    pvalue_adj=data$padj)
  write.table(res,file=fn,row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
  return(res)
}

