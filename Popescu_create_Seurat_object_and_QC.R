library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(patchwork)

##### creating Seurat Objects ####
# reads the outs from CellRanger and creates the Seurat Obj
# it filters for genes detected at least in 3 cells
# and for cells that express at least 300 genes
data10x <- Read10X('/Users/majagabric/Documents/PhD/Data analysis/Popescu et al., 2019/Popescu et al., 2019/aggrData/Pop_aggr/count/filtered_feature_bc_matrix')
seurat_obj <- CreateSeuratObject(data10x, 
                                     min.cells = 3, min.features = 200, 
                                     project = "Popescu2019")

##### Add metadata #####
aggr_metadata <- read.csv('/Users/majagabric/Documents/PhD/Data analysis/Popescu et al., 2019/Popescu et al., 2019/aggrData/Pop_aggr/aggregation.csv')

barcode_suffix <- sub(".*-(\\d+)$", "\\1", colnames(seurat_obj))

# map barcode suffic to the sample name 
sample_map <- c(
  "1" = "YS_w4_P19_1",
  "2" = "YS_w4_P19_2",
  "3" = "YS_w4_P19_3",
  "4" = "YS_w4_P19_4",
  "5" = "YS_w4_P19_5",
  "6" = "YS_w7_P19_1",
  "7" = "YS_w7_P19_2",
  "8" = "FL_w7_P19_1",
  "9" = "FL_w7_P19_2",
  "10" = "FL_w7_P19_3",
  "11" = "FL_w7_P19_4",
  "12" = "FL_w7_P19_5",
  "13" = "FL_w7_P19_6",
  "14" = "FL_w8_P19_1",
  "15" = "FL_w8_P19_2",
  "16" = "FL_w8_P19_3",
  "17" = "FL_w8_P19_4",
  "18" = "FL_w8_P19_5",
  "19" = "FL_w8_P19_6",
  "20" = "FL_w8_P19_7",
  "21" = "FL_w9_P19_1",
  "22" = "FL_w9_P19_2",
  "23" = "FL_w9_P19_3",
  "24" = "FL_w9_P19_4",
  "25" = "FL_w9_P19_5",
  "26" = "FL_w9_P19_6",
  "27" = "FL_w9_P19_7",
  "28" = "FL_w9_P19_8",
  "29" = "FL_w9_P19_9",
  "30" = "FL_w9_P19_10",
  "31" = "FL_w9_P19_11",
  "32" = "FL_w11_P19_1",
  "33" = "FL_w11_P19_2",
  "34" = "FL_w12_P19_1",
  "35" = "FL_w12_P19_2",
  "36" = "FL_w12_P19_3",
  "37" = "FL_w12_P19_4",
  "38" = "FL_w12_P19_5",
  "39" = "FL_w12_P19_6",
  "40" = "FL_w13_P19_1",
  "41" = "FL_w13_P19_2",
  "42" = "FL_w13_P19_3",
  "43" = "FL_w13_P19_4",
  "44" = "FL_w13_P19_5",
  "45" = "FL_w14_P19_1",
  "46" = "FL_w14_P19_2",
  "47" = "FL_w14_P19_3",
  "48" = "FL_w16_P19_1",
  "49" = "FL_w16_P19_2",
  "50" = "FL_w16_P19_3",
  "51" = "FL_w16_P19_4",
  "52" = "FL_w16_P19_5",
  "53" = "FL_w16_P19_6",
  "54" = "FL_w17_P19_1",
  "55" = "FL_w17_P19_2",
  "56" = "FL_w17_P19_3"
)

# use the barcode suffix to get the correct sample name and create a vector sample_id
sample_id <- sample_map[barcode_suffix]

# make sure everything is matched correctly - set names that match seurat object columns
names(sample_id) <- colnames(seurat_obj)

# add sample_id to the metadata
seurat_obj <- AddMetaData(seurat_obj, metadata = sample_id, col.name = "sample_id")

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
sex_metadata <- c("F","F","F","F","F","F","F","F","F","F","F","F","F","M","M","F","F","F","F","F","M","M","M","F","F","F","F","F","F","F","F","M","M","F","F","F","F","F","F","M","M","M","M","M","M","M","M","M","M","F","F","F","F","F","F","F"
)
seurat_obj$Sex <- sex_metadata[seurat_obj$Sample]

# Add Week_gestation and Tissue
seurat_obj$Week_gestation <- sub(".*_w([0-9]+)_.*", "\\1", seurat_obj$sample_id)
# Convert to numeric
seurat_obj$Week_gestation <- as.numeric(seurat_obj$Week_gestation)

# FL samples
seurat_FL <- subset(seurat_obj, subset = grepl("^FL", sample_id))
# YS samples
seurat_YS <- subset(seurat_obj, subset = grepl("^YS", sample_id))
# add tissue metadata
seurat_obj$Tissue <- ifelse(grepl("^FL", seurat_obj$sample_id), "FL", "YS")

##### Save Seurat object #####
saveRDS(seurat_obj, file="Popescu_seurat_obj_preQC.rds")

##### Load Seurat object #####
seurat_obj <- readRDS('/Users/majagabric/Documents/PhD/Data analysis/Popescu et al., 2019/Popescu et al., 2019/Popescu_seurat_obj_preQC.rds')

##### QC #####
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Violin plots
p <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
p[[3]] <- p[[3]] + ylim(0, 25)
p

# Check correlation of nFeature_RNA,nCount_RNA, and percent.mt
# plot1
plot1 <- FeatureScatter(seurat_obj,
                        feature1 = "nFeature_RNA",
                        feature2 = "percent.mt")
threshold_x <- 200 
threshold_x2 <- 6000
threshold_y <- 0
threshold_y2 <- 15
plot1 +
  geom_vline(xintercept = threshold_x, color = "forestgreen") +
  geom_vline(xintercept = threshold_x2, color = "purple") +
  geom_hline(yintercept = threshold_y, color = "forestgreen") +
  geom_hline(yintercept = threshold_y2, color = "purple") +
  scale_x_continuous(breaks = c(scales::pretty_breaks()(plot1$data$nFeature_RNA), threshold_x, threshold_x2)) +
  scale_y_continuous(breaks = c(scales::pretty_breaks()(plot1$data$percent.mt), threshold_y, threshold_y2)) +
coord_cartesian(ylim = c(0, 30)) # this only zooms in on the Y-axis, use and adjust c(x, x) as needed

# plot2
plot2 <- FeatureScatter(seurat_obj,
                        feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA")
threshold2_y2 <- 6000
threshold2_y <- 200
plot2 +
  geom_hline(yintercept = threshold2_y2, color = "purple") +
  geom_hline(yintercept = threshold2_y, color = "forestgreen") +
  scale_y_continuous(breaks = c(scales::pretty_breaks()(plot2$data$nFeature_RNA),threshold2_y,threshold2_y2)) #+
coord_cartesian(clip = "off")  

# set thresholds
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)


# SCTansform
seurat_obj <- SCTransform(seurat_obj, vst.flavor = "v2", verbose = TRUE)

# PCA
seurat_obj <- RunPCA(seurat_obj, npcs = 50)

ElbowPlot(seurat_obj, ndims = ncol(Embeddings(seurat_obj, "pca")))


##### Load Seurat object #####
seurat_obj <- readRDS('/Users/majagabric/Documents/PhD/Data analysis/Calvanese et al., 2022/Calvanese et al., 2022/Calvanese_suerat_obj_postPCA.rds')


# UMAP and tSNE
seurat_obj <- RunTSNE(seurat_obj, dims = 1:28)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:28)

plot1 <- TSNEPlot(seurat_obj)
plot2 <- UMAPPlot(seurat_obj)
plot1 + plot2

##### Save Seurat object #####
saveRDS(seurat_obj, file="Calvanese_Seurat_object.rds")

