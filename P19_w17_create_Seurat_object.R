library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)

##### creating Seurat Objects ####
# reads the outs from CellRanger and creates the Seurat Obj
# it filters for genes detected at least in 3 cells
# and for cells that express at least 300 genes
data10x <- Read10X('/Users/majagabric/Documents/PhD/Data analysis/Popescu et al., 2019/Popescu et al., 2019/aggrData/Pop17_aggr/count/filtered_feature_bc_matrix')
seurat_obj_w17 <- CreateSeuratObject(data10x, 
                                     min.cells = 3, min.features = 200, 
                                     project = "Popescu2019_w17")

##### Add metadata #####
aggr_metadata <- read.csv('/Users/majagabric/Documents/PhD/Data analysis/Popescu et al., 2019/Popescu et al., 2019/aggrData/Pop17_aggr/aggregation.csv')

barcode_suffix <- sub(".*-(\\d+)$", "\\1", colnames(seurat_obj_w17))

# map barcode suffic to the sample name 
sample_map <- c(
  "1" = "FL_w17_P19_1",
  "2" = "FL_w17_P19_2",
  "3" = "FL_w17_P19_3"
)

# use the barcode suffix to get the correct sample name and create a vector sample_id
sample_id <- sample_map[barcode_suffix]

# make sure everything is matched correctly - set names that match seurat object columns
names(sample_id) <- colnames(seurat_obj_w17)

# add sample_id to the metadata
seurat_obj_w17 <- AddMetaData(seurat_obj_w17, metadata = sample_id, col.name = "sample_id")

# Get the cell tags
barcodes <- colnames(seurat_obj_w17)

# Cell tags are made of 16 bases, a dash and the sample number, so you can simply split by "-" and take whatever's on the right
# TGGCCAGGTTATCACG-1
# We use sapply to go through each barcode and apply the strsplit function, which splits a string based on a delimiter. 
# strsplit will return an array of length 2; the first element will be the actual cell tag, and the second (which we will return) is the sample number
# We save the result of that in a metadata called Sample
seurat_obj_w17$Sample <- sapply(strsplit(barcodes, "-"), function(x){return(x[2])})

# We can now create the other metadata depending on Sample
# Define sex for each sample then create the Sex metadata accordingly
seurat_obj_w17$Sample <- as.numeric(seurat_obj_w17$Sample)
sex_metadata <- "F"
seurat_obj_w17$Sex <- sex_metadata[seurat_obj_w17$Sample]

# Add Week_gestation and Tissue
seurat_obj_w17$Week_gestation <- 17
seurat_obj_w17$Tissue <- "FL"

##### Save Seurat object #####
saveRDS(seurat_obj_w17, file="P19_w17_seurat_obj_preQC.rds")

##### Load Seurat object #####
seurat_obj_w17 <- readRDS('/Users/majagabric/Documents/PhD/Data analysis/Popescu et al., 2019/Popescu et al., 2019/P19_w17_seurat_obj_preQC.rds')


