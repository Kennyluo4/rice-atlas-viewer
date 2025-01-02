##suspended, plot umap using Scanpy

## libraries========================================
library(Seurat)
library(scCustomize)
library(ggplot2)

# load arguments
arg <- commandArgs(T)
print(arg)

sobj_file <- as.character(arg[1])
##Late_seed.obj

feature_ID_file <- as.character(arg[2])
##geneIDs

output_figure<- as.character(arg[3])

## testing 
sobj_file <- 'Rplot/data/sobjs/Crown_root.obj'
feature_ID_file <- 'Rplot/temp/umap_features.txt'
output_figure <- 'Rplot/temp/test_umap.tiff'

## load the data ========================================

sobj <- readRDS(sobj_file)
feature_df <- read.delim(feature_ID_file)
feature_IDs <- 
  
## plot the figure ======================================
p <- FeaturePlot(sobj, features = feature_IDs, slot = 'count')

ggsave(output_figure, p, height = 10, width = 10)