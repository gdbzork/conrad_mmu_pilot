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
  pilot_log$info("Counting fastq reads in %s...",fn)
  cmd <- sprintf("zcat %s | wc -l",fn)
  s <- system(cmd,intern=TRUE)
  reads <- as.numeric(s) / 4 # because 4 lines per Fastq read
  return(reads)
}

reads_in_bam <- function(fns) {
  return(as.numeric(sapply(fns,reads_in_single_bam)))
}