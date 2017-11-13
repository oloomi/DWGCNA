find_diff_exprs <- function(pre.data.dir, post.data.dir, preprocess.method, class.label) {
  library(multtest)
  pre.exprs <- read.table(file = 
                            paste(pre.data.dir,"/exprs-vals-", preprocess.method, ".txt", sep=''))
  post.exprs <- read.table(file = 
                             paste(post.data.dir,"/exprs-vals-", preprocess.method, ".txt", sep=''))
  # Attaching pre- and post-data
  pre.post.exprs <- cbind(pre.exprs, post.exprs)
  # Check for normality
#   qqnorm(data.exprs)
  # Run differential expression analysis
  mtp.data <- MTP(X=pre.post.exprs, Y=class.label, typeone="fdr", B=1000)
  # Sort genes by p-value(FDR)
  order.ind <- order(mtp.data@adjp)
  diff.genes <- data.frame(row.names=rownames(pre.post.exprs)[order.ind], adjp=mtp.data@adjp[order.ind])
  
  write.table(diff.genes, file = 
              paste("common/diff-exprs-", preprocess.method, ".txt", sep=''), 
              quote = FALSE, sep = "\t")
  save(diff.genes, file="diff.genes.RData")
}