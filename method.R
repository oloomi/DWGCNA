ptm <- proc.time()

library(affy)

source("preprocess.R")
source("differential-exprs.R")
source("coexpression-network.R")
source("integration.R")

# Pre-treatment data
data.dir <- "pre"
# Loading raw .CEL files and associating them with our Custom CDF library
data.cel <- ReadAffy(celfile.path=paste(data.dir,"/CEL", sep=''), cdfname="HGU133A2_HS_ENTREZG")
# data.cel <- ReadAffy(celfile.path=paste(data.dir,"/CEL", sep=''))
# MAS 5.0 detection cal
detection_call(data.cel, data.dir, 75)
# MAS 5.0 normalization
preprocess_MAS5(data.cel, data.dir)
# RMA normalization
preprocess_RMA(data.cel, data.dir)
# Extract expression values of present genes
extract_present_exprs(data.dir, "mas5")
extract_present_exprs(data.dir, "rma")

# Repeating the above steps for post-treatment data
data.dir <- "post"
data.cel <- ReadAffy(celfile.path=paste(data.dir,"/CEL", sep=''), cdfname="HGU133A2_HS_ENTREZG")
# data.cel <- ReadAffy(celfile.path=paste(data.dir,"/CEL", sep=''))
detection_call(data.cel, data.dir, 75)
preprocess_MAS5(data.cel, data.dir)
preprocess_RMA(data.cel, data.dir)
extract_present_exprs(data.dir, "mas5")
extract_present_exprs(data.dir, "rma")

# Find the union of pre- and post-treatment data present genes
extract_union_exprs("pre", "post", "mas5")
extract_union_exprs("pre", "post", "rma")

# Find list of differentially expressed genes
# class.label <- rep(c(0,1), each=22)
class.label <- rep(c(0,1), c(18,19))
find_diff_exprs("pre", "post", "rma", class.label)

# Find network statistics for threshold selection
library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads()
# exprs.type: present, union, all
# preprocess.method: rma, mas5
# nSamples <- 22

# Find threshold statistics
cut.vector <- seq(0.1, 0.9, by=0.05)
# cut.vector <- seq(0.6, 0.9, by=0.02)
nSamples <- 18
data.dir <- "pre"
threshold_stat(data.dir, "union", "mas5", cut.vector)
# threshold_stat(data.dir, "union", "rma", cut.vector)
nSamples <- 19
data.dir <- "post"
threshold_stat(data.dir, "union", "mas5", cut.vector)
# threshold_stat(data.dir, "union", "rma", cut.vector)

# Decide on threshold, then construct the network
# Networks for MAS5.0
nSamples <- 18
construct_network("pre", "union", "mas5", 0.72)
nSamples <- 19
construct_network("post", "union", "mas5", 0.7)
# Networks for RMA
# construct_network("pre", "union", "rma", 0.65)
# construct_network("post", "union", "rma", 0.65)

# elapsed <- proc.time() - ptm

# Create network graph
pre.net <- create_net_graph("pre", "union", "mas5", 0.72)
post.net <- create_net_graph("post", "union", "mas5", 0.7)

# Find differences between pre and post networks
find_net_diff(pre.net, post.net, "common/net-diff-mas5-0.72-0.7-order1")

# Compare E2F1 gene neighbors in MAS pre and post network
find_node_neighbor(pre.net, "1869", "pre/neighbors-1869-mas5-union-0.72-0.7")
find_node_neighbor(post.net, "1869", "post/neighbors-1869-mas5-union-0.72-0.7")

elapsed <- proc.time() - ptm
