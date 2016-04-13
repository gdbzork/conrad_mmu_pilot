reads_in_single_fastq <- function(fn) {
  pilot_log$info("Counting fastq reads in %s...",fn)
  cmd <- sprintf("zcat %s | wc -l",fn)
  s <- system(cmd,intern=TRUE)
  reads <- as.numeric(s) / 4 # because 4 lines per Fastq read
  return(reads)
}

reads_in_fastq <- function(fns) {
  return(as.numeric(sapply(fns,reads_in_single_fastq)))
}

reads_in_single_bam <- function(fn) {
  pilot_log$info("Counting bam reads in %s...",fn)
  cmd <- sprintf("samtools flagstat %s",fn)
  s <- system(cmd,intern=TRUE)
  n <- strsplit(s[5]," ")[[1]][1]
  reads <- as.numeric(n) / 2 # pairs
  return(reads)
}

reads_in_bam <- function(fns) {
  pilot_log$info("wd: %s",getwd())
  return(as.numeric(sapply(fns,reads_in_single_bam)))
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

de_analysis <- function(data,model,subset_vec,alpha=0.05) {
  de_obj <- DESeqDataSet(data,design = model)
  de_obj <- subset(de_obj,select=subset_vec)
  de_diff <- DESeq(de_obj)
  de_results <- results(de_diff,alpha=0.05)
  return(de_results)
}
