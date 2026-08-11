library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)

sce <- readVisium('9w/outs/binned_outputs/square_002um/')

# Pre-processing
set.seed(42)
sce <- spatialPreprocess(sce, platform="VisiumHD", 
                         n.PCs=15, n.HVGs=2000, log.normalize=TRUE)

# Select the number of clusters
sce <- qTune(sce, qs=seq(2, 10), platform="VisiumHD", d=15, cores=8)
pdf("test.figures/BayesSpace_qTune.pdf")
qPlot(sce)
dev.off()