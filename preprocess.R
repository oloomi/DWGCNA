detection_call <- function(data.cel, data.dir, present.percent) {
  # Finding detection calls
  data.mas5calls = mas5calls(data.cel)
  # Get the actual P/M/A calls
  data.mas5calls.pma = exprs(data.mas5calls)
  write.table(data.mas5calls.pma, file = paste(data.dir, "/mas5calls.txt", sep=''), 
              quote = FALSE, sep = "\t")
  # Removing control probesets
  control.prbs <- grep("AFFX",rownames(data.mas5calls.pma))
  data.mas5calls.pma <- data.mas5calls.pma[-control.prbs,]  
  # Counting P/M/A calls number
  calls.count <- pma_calls_dist(data.mas5calls.pma, data.dir)
  
  # Genes which are present in more than present.num samples, will be considered present
  present.num <- round(ncol(data.mas5calls.pma) * present.percent / 100)
  # Getting names of present genes
  present.genes <- row.names(data.mas5calls.pma)[(calls.count[,1] >= present.num)]
  write.table(present.genes, file = paste(data.dir,"/present-genes-names.txt", sep=''), 
              quote = FALSE, row.names=FALSE, col.names=FALSE, sep = "\t")
}

preprocess_MAS5 <- function(data.cel, data.dir) {
  # MAS5.0 Normalization
  data.eset = mas5(data.cel)
  data.exprs = exprs(data.eset)
#   data.exprs = log(data.exprs, 2)
  # Removing control probesets
  control.prbs <- grep("AFFX",rownames(data.exprs))
  data.exprs <- data.exprs[-control.prbs,]
  # Saving the expression value of all genes
  write.table(data.exprs, file = paste(data.dir,"/exprs-vals-mas5.txt", sep=''), 
              quote = FALSE, sep = "\t")
  
}

preprocess_RMA <- function(data.cel, data.dir) {
  # RMA
  data.eset = rma(data.cel)
  data.exprs = exprs(data.eset)
  # Removing control probesets
  control.prbs <- grep("AFFX",rownames(data.exprs))
  data.exprs <- data.exprs[-control.prbs,]
  # Saving the expression value of all genes
  write.table(data.exprs, file = paste(data.dir,"/exprs-vals-rma.txt", sep=''), 
              quote = FALSE, sep = "\t")
}

extract_present_exprs <- function(data.dir, preprocess.method) {
  present.genes.names <- read.table(file = 
                                    paste(data.dir,"/present-genes-names.txt", sep=''))
  all.genes.exprs <- read.table(file = 
                                paste(data.dir,"/exprs-vals-", preprocess.method, ".txt", sep=''))
  # Extracting expression values of present genes
  present.genes.exprs <- all.genes.exprs[present.genes.names$V1, ]
  write.table(present.genes.exprs, file = 
              paste(data.dir,"/exprs-vals-", preprocess.method, "-prsnt.txt", sep=''), 
              quote = FALSE, sep = "\t")
}

extract_union_exprs <- function(pre.data.dir, post.data.dir, preprocess.method) {
  
  # Find union of list of present genes in pre- and post-treatment data
  pre.present.genes <- read.table(file = 
                                    paste(pre.data.dir,"/present-genes-names.txt", sep=''))
  post.present.genes <- read.table(file = 
                                     paste(post.data.dir,"/present-genes-names.txt", sep=''))
  union.present.genes <- union(pre.present.genes$V1, post.present.genes$V1)
  write.table(union.present.genes, file = "common/union-genes-names.txt", 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Extracting expression values of union genes from pre- and post-treatment data
  pre.exprs <- read.table(file = 
                                 paste(pre.data.dir,"/exprs-vals-", preprocess.method, ".txt", sep=''))
  post.exprs <- read.table(file = 
                                  paste(post.data.dir,"/exprs-vals-", preprocess.method, ".txt", sep=''))
  pre.exprs.union <- pre.exprs[union.present.genes, ]
  post.exprs.union <- post.exprs[union.present.genes, ]
  write.table(pre.exprs.union, file = 
              paste(pre.data.dir,"/exprs-vals-", preprocess.method, "-union.txt", sep=''), 
              quote = FALSE, sep = "\t")
  write.table(post.exprs.union, file = 
              paste(post.data.dir,"/exprs-vals-", preprocess.method, "-union.txt", sep=''), 
              quote = FALSE, sep = "\t")
}

pma_calls_dist <- function(data.mas5calls.pma, data.dir) {
  library("ggplot2")  
  #  3 columns holding number of Present, Absent, Marginal calls
  calls.count = matrix(data = 0, nrow = nrow(data.mas5calls.pma), ncol = 3)
  
  # P/A/M (not pma)
  for(i in 1:nrow(data.mas5calls.pma))
    for(j in 1:ncol(data.mas5calls.pma)) {
      if(data.mas5calls.pma[i,j] == "P") {
        calls.count[i,1] <- calls.count[i,1] + 1
      } else if(data.mas5calls.pma[i,j] == "A") {
        calls.count[i,2] <- calls.count[i,2] + 1
      } else {
        calls.count[i,3] <- calls.count[i,3] + 1
      }
    }
  
  # Drawing histogram of number of calls
  calls.count.df <- rbind(data.frame(Call = "Present", count = calls.count[,1]),
                          data.frame(Call = "Absent", count = calls.count[,2]),
                          data.frame(Call = "Marginal", count = calls.count[,3]))
  calls.hist <- ggplot(calls.count.df, aes(x = count, fill = Call)) + 
    geom_histogram(position = "dodge", binwidth = 5) +
    scale_fill_manual(values = c("seagreen3", "firebrick1", "royalblue2")) +
    xlab("Detection Calls(No. of Samples)") + ylab("Frequency(No. of Genes)") + 
    theme(axis.text.x=element_text(size=12)) + theme(axis.text.y=element_text(size=12))
  print(calls.hist)
  ggsave(paste(data.dir,"/pma-calls-hist.png", sep=''), width=4, height=4, dpi=100)
  
  # Writing histogram data in file
  # NOTE! for binwidth = 1
  end.bin <- ncol(data.mas5calls.pma)
  p.hist <- hist(calls.count[,1], breaks = c(0:(end.bin + 1)), right = FALSE, plot = FALSE)
  a.hist <- hist(calls.count[,2], breaks = c(0:(end.bin + 1)), right = FALSE, plot = FALSE)
  m.hist <- hist(calls.count[,3], breaks = c(0:(end.bin + 1)), right = FALSE, plot = FALSE)
  pam.hist <- data.frame(cbind(0:end.bin, p.hist$counts, p.hist$density, a.hist$counts, 
                               a.hist$density, m.hist$counts, m.hist$density))
  colnames(pam.hist) <- c("Sample.counts", "P.calls", "P.density", "A.calls", "A.density", 
                          "M.calls", "M.density")
  write.table(pam.hist, file = paste(data.dir,"/mas5calls-summary.txt", sep=''), sep = "\t", 
              row.names = FALSE, quote = FALSE) 
  return(calls.count)
}
