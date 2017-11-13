create_net_graph <- function(data.dir, exprs.type, preprocess.method, cutoff) {
  library('igraph')  
  edge.list <- read.table(file=
                            paste(data.dir,"/edges-", preprocess.method, "-", 
                                  exprs.type, "-", cutoff, ".txt", sep=''), 
                          header=TRUE)
  edge.list <- edge.list[,1:3]
  edge.list <- as.data.frame(sapply(edge.list,gsub,pattern="_at",replacement=""))
  
  g <- graph.data.frame(edge.list, directed=FALSE)
  return(g)
}

find_node_neighbor <- function(net.graph, node.name, output.path) {
  subg <- graph.neighborhood(net.graph, 1, node.name)
  neig <- subg[[1]]
  max.deg <- max(degree(neig))
  min.deg <- min(degree(neig))
  V(neig)$size <- 10 +  (degree(neig)-min.deg)/(max.deg-min.deg) * 10
  V(neig)$frame.color <- "SkyBlue2"
  V(neig)[node.name]$color <- "yellow"
  E(neig)$color <- "gray"
  E(neig)[incident(neig, node.name)]$color <- "yellow"
  #tkplot(neig, layout=layout.kamada.kawai, canvas.width=700 , canvas.height= 700)
#   png(file=paste(data.dir,"/neighbors-", node.name, "-", preprocess.method, "-", 
#                  exprs.type, "-", cutoff, ".png", sep=''))
  png(file=paste(output.path, ".png", sep=''))
  plot.igraph(neig, layout=layout.kamada.kawai)
#   tkplot(neig, layout=layout.kamada.kawai, canvas.width=700 , canvas.height= 700)
  dev.off()
  
#   library("hgu133a2hsentrezg.db")
#   getSYMBOL(neighbors(neig, past(node.name, "_at", sep='')))[1], "hgu133a2hsentrezg.db")
}

find_net_diff <- function(pre.net, post.net, output.path) {
#   union.gene.names <- read.table("common/union-genes-names.txt")
  diff.exprs <- read.table("common/diff-exprs-rma.txt", header=TRUE)
  row.names(diff.exprs) <- sapply(rownames(diff.exprs), gsub, pattern="_at", replacement="")
#   diff.exprs <- cbind(gene.id=rownames(diff.exprs), diff.exprs)
#   
#   # --> either one of the two below
#   # All genes sorted by fdr
#   diff.exprs.fdr <- diff.exprs
#   # Find genes with FDR < 0.05
# #   diff.exprs.fdr <- diff.exprs[diff.exprs$adjp < 0.05, ]  
#   
#   # Find which have FDR < 0.05 and are present in networks
#   net.fdr.genes <- intersect(union.gene.names$V1, diff.exprs.fdr$gene.id) 
#   net.fdr.genes <- sapply(net.fdr.genes, gsub, pattern="_at", replacement="")
  
  
  # Get names of network vertices(names of genes present in network)
  pre.net.genes <- get.vertex.attribute(pre.net, "name")
  post.net.genes <- get.vertex.attribute(post.net, "name")
  net.fdr.genes <- union(pre.net.genes, post.net.genes)
  
  neighbors.count = matrix(data = 0, nrow = length(net.fdr.genes), ncol = 5)
  
  for(i in 1:length(net.fdr.genes)) {
    # Because maybe some genes are single in network and have no edges,
    # check if they exist in network
    # Pre-treatment neighbors
    if(net.fdr.genes[i] %in% pre.net.genes) {
        pre.neighbors <- neighborhood(pre.net, 1, net.fdr.genes[i])      
        neighbors.count[i, 1] <- length(pre.neighbors[[1]])
    } else {
      neighbors.count[i, 1] <- 0
    }
    # Post-treatment neighbors
    if(net.fdr.genes[i] %in% post.net.genes) {
      post.neighbors <- neighborhood(post.net, 1, net.fdr.genes[i])
      neighbors.count[i, 2] <- length(post.neighbors[[1]])
    } else {
      neighbors.count[i, 2] <- 0
    }
    # Pre and post common neighbors
    if( (net.fdr.genes[i] %in% pre.net.genes) & (net.fdr.genes[i] %in% post.net.genes) ) {
      common.neighbors <- intersect(pre.neighbors[[1]], post.neighbors[[1]])
      neighbors.count[i, 3] <- length(common.neighbors)
    } else {
      neighbors.count[i, 3] <- 0
    }
    # Attaching FDR of the current gene
    neighbors.count[i, 4] <- diff.exprs[net.fdr.genes[i], ]
    # Calculating network difference
    neighbors.count[i, 5] <- abs(neighbors.count[i, 1] - neighbors.count[i, 2]) / 
                                (neighbors.count[i, 1] + neighbors.count[i, 2])
  }
  # Write the results to file
  neighbors.count.df <- data.frame(gene.id=net.fdr.genes, pre=neighbors.count[,1],
                                   post=neighbors.count[,2], common=neighbors.count[,3], 
                                   fdr=neighbors.count[,4], net.diff = neighbors.count[,5])
  write.csv(neighbors.count.df, file=paste(output.path, ".csv", sep=''), 
            row.names=FALSE, quote=FALSE)
}
