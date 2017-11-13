
# Reading gene expression profiles
pre.data.dir <- "pre"
post.data.dir <- "post"
preprocess.method <- "rma"
pre.exprs <- read.table(file = 
                          paste(pre.data.dir,"/exprs-vals-", preprocess.method, ".txt", sep=''))
post.exprs <- read.table(file = 
                           paste(post.data.dir,"/exprs-vals-", preprocess.method, ".txt", sep=''))
# Attaching pre- and post-data
pre.post.exprs <- cbind(pre.exprs, post.exprs)
# Remove _at suffix
rownames(pre.post.exprs) <- sapply(rownames(pre.post.exprs),gsub,pattern="_at",replacement="")

# Read genes information file
markers <- read.csv(file="common/net-diff-mas5-0.72-0.7-order1.csv")


#---- Select one of these two: ----#
# Sort by lowest to highest FDR and get order indices
markers.sorted <- markers[order(markers$fdr), ]
# Sort by highest to lowest network difference and get order indices
markers.sorted <- markers[order(markers$net.diff, decreasing=TRUE), ]

# Get gene IDs of the order
marker.genes <- markers.sorted[, 1]

# Classification
library('CMA')

# class labels
class.label <- rep(c(0,1), c(18,19))

# Getting expression profiles of top markers
exprs.data <- subset(pre.post.exprs, subset=(rownames(pre.post.exprs) %in% marker.genes[1:200]))
exprs.data <- t(data.matrix(exprs.data))

exprs.data <- t(data.matrix(pre.post.exprs))

# Leaving-one-out cross-validation
# loo.dat <- GenerateLearningsets(y = class.label, method = "LOOCV")
# k-fold cross-validation
set.seed(321)
kfold.dat <- GenerateLearningsets(y = class.label, method = "CV", fold = 3, strat = TRUE)
 
varsel <- GeneSelection(X = exprs.data, y = class.label, learningsets = kfold.dat, method = "welch.test")

# pKNN
knn.cl.outs <- classification(X = exprs.data, y = class.label, learningsets = kfold.dat, 
                            classifier = pknnCMA, k = 3, genesel = varsel, nbgene = 2)

knn.cl.outs <- classification(X = exprs.data, y = class.label, learningsets = kfold.dat, 
                              classifier = pknnCMA, k = 3)

result <- join(knn.cl.outs)
ftable(result)
roc(result, sub="FDR - Top 2")

eval <- evaluation(knn.cl.outs, measure="auc")
show(eval)
boxplot(eval)

