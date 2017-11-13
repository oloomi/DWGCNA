threshold_stat <- function(data.dir, exprs.type, preprocess.method, cut.vector) {
  # exprs.type: present, union, all
  # preprocess.method: rma, mas5
  data.exprs <- read.table(file=
                             paste(data.dir,"/exprs-vals-", preprocess.method, 
                                   "-", exprs.type, ".txt", sep=''))
  data.exprs <- t(data.exprs)
  # Similarity matrix by biweight mid-correlation
  simil <- bicor(data.exprs, maxPOutliers=0.02)
#  simil <- cor(data.exprs)
  simil.abs <- abs(simil)
#   simil.hist <- simil[(simil > 0.8) & (simil != 1)]
#   hist(simil.hist)  
  # Find statistics for different thresholds
  nSamples = nrow(data.exprs)
  threshold.stat <- pickHardThreshold.fromSimilarity(
    simil.abs,
    RsquaredCut = 0.7,
    cutVector = cut.vector,
    moreNetworkConcepts = TRUE,
    removeFirst = FALSE, nBreaks = 10)
  
  write.table(threshold.stat[2]$fitIndices, 
              file=paste(data.dir,"/threshold-stat-", preprocess.method, "-", 
                         exprs.type, ".txt", sep=''), 
              quote=FALSE, row.names=FALSE)
#   return(simil)
}

construct_network <- function(data.dir, exprs.type, preprocess.method, cutoff) {
  data.exprs <- read.table(file=
                             paste(data.dir,"/exprs-vals-", preprocess.method, 
                                   "-", exprs.type, ".txt", sep=''))
  data.exprs <- t(data.exprs)
  # Similarity matrix by biweight mid-correlation
  simil <- bicor(data.exprs, maxPOutliers=0.02)
#   simil <- cor(data.exprs)
  simil.abs <- abs(simil)
  exportNetworkToCytoscape(
    simil.abs, 
    edgeFile=paste(data.dir,"/edges-", preprocess.method, "-", exprs.type, "-", 
                   cutoff, ".txt", sep=''), 
    nodeFile=paste(data.dir,"/nodes-", preprocess.method, "-", exprs.type, "-", 
                   cutoff, ".txt", sep=''), 
    weighted=TRUE, threshold=cutoff)
  
  #   adjacency <- simil.abs
  #   adjacency[(adjacency < 0.6) | (adjacency == 1)] <- 0
  #   write.table(adjacency, "post/adjacency-0.6-cut.txt", sep='\t')
}

calculate_pvalues <- function() {
  cor.vals <- seq(0.1,0.9,0.05)
  pre.cor.fisher <- corPvalueFisher(cor.vals, 18, twoSided = TRUE)
  post.cor.fisher <- corPvalueFisher(cor.vals, 19, twoSided = TRUE)
  
  pvals <- data.frame(cor.vals=cor.vals, pre.cor.fisher=pre.cor.fisher, post.cor.fisher=post.cor.fisher)
  write.csv(pvals, file="common/pcor-fisher-pvalues.csv")
  
  fine.cor.vals <- seq(0.6,0.9,0.02)
  pre.fine.cor.fisher <- corPvalueFisher(fine.cor.vals, 18, twoSided = TRUE)
  post.fine.cor.fisher <- corPvalueFisher(fine.cor.vals, 19, twoSided = TRUE)
  
  fine.pvals <- data.frame(fine.cor.vals=fine.cor.vals, pre.fine.cor.fisher=pre.fine.cor.fisher, 
                      post.fine.cor.fisher=post.fine.cor.fisher)
  
  write.csv(fine.pvals, file="common/fine-pcor-fisher-pvalues.csv")
}
