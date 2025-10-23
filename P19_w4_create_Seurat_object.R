library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)

##### creating Seurat Objects ####
# reads the outs from CellRanger and creates the Seurat Obj
# it filters for genes detected at least in 3 cells
# and for cells that express at least 300 genes
data10x <- Read10X('/Users/majagabric/Documents/PhD/Data analysis/Popescu et al., 2019/Popescu et al., 2019/aggrData/Pop4_aggr/count/filtered_feature_bc_matrix')
seurat_obj <- CreateSeuratObject(data10x, 
                                 min.cells = 3, min.features = 200, 
                                 project = "Popescu_w4_2019")

##### Add sample_id metadata #####
aggr_metadata <- read.csv('/Users/majagabric/Documents/PhD/Data analysis/Popescu et al., 2019/Popescu et al., 2019/aggrData/Pop4_aggr/aggregation.csv')

barcode_suffix <- sub(".*-(\\d+)$", "\\1", colnames(seurat_obj))

# map barcode suffic to the sample name 
sample_map <- c(
  "1" = "YS_w4_P19_1",
  "2" = "YS_w4_P19_2",
  "3" = "YS_w4_P19_3",
  "4" = "YS_w4_P19_4",
  "5" = "YS_w4_P19_5"
)

# use the barcode suffix to get the correct sample name and create a vector sample_id
sample_id <- sample_map[barcode_suffix]

# make sure everything is matched correctly - set names that match seurat object columns
names(sample_id) <- colnames(seurat_obj)

# add sample_id to the metadata
seurat_obj <- AddMetaData(seurat_obj, metadata = sample_id, col.name = "sample_id")

##### Save Seurat object #####
saveRDS(seurat_obj, file="P19_w4_suerat_obj.rds")

##### Load Seurat object #####
seurat_obj <- readRDS('/Users/majagabric/Documents/PhD/Data analysis/Popescu et al., 2019/Popescu et al., 2019/P19_w4_suerat_obj.rds')

##### Add metadata #####
# Get the cell tags
barcodes <- colnames(seurat_obj)

# Cell tags are made of 16 bases, a dash and the sample number, so you can simply split by "-" and take whatever's on the right
# TGGCCAGGTTATCACG-1
# We use sapply to go through each barcode and apply the strsplit function, which splits a string based on a delimiter. 
# strsplit will return an array of length 2; the first element will be the actual cell tag, and the second (which we will return) is the sample number
# We save the result of that in a metadata called Sample
seurat_obj$Sample <- sapply(strsplit(barcodes, "-"), function(x){return(x[2])})

# We can now create the other metadata depending on Sample
# Define sex for each sample then create the Sex metadata accordingly
seurat_obj$Sample <- as.numeric(seurat_obj$Sample)
sex_metadata <- c("F", "F", "F", "F", "F")
seurat_obj$Sex <- sex_metadata[seurat_obj$Sample]
seurat_obj$Week_gestation <- 4
seurat_obj$Tissue <- "YS"


