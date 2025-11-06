suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(ggrepel)
  library(gridExtra)
  library(patchwork)
  library(stringr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(tidyverse)
  library(dplyr)
  library(cowplot)
  library(CellChat)
})

seurat_obj_w7_2025 <- readRDS('/Users/majagabric/Documents/PhD/Data analysis/Popescu et al., 2019/Popescu et al., 2019/P19_w7_seurat_obj.rds')
cellchat_obj_w7_liver <- readRDS('/Users/majagabric/Documents/PhD/Data analysis/Popescu data/5_to_17_analysis_Pop_Wes/objects/Liver_wk7_cellchatobj.rds')
seurat_obj_stromal_w7 <- readRDS('/Users/majagabric/Documents/PhD/Data analysis/Popescu data/5_to_17_analysis_Pop_Wes/objects/Stromal_Wk7.rds')

##### Liver genes ####
identifyOverExpressedGenes(
  cellchat_obj_w7_liver,
  data.use = NULL,
  group.by = NULL,
  idents.use = NULL,
  invert = FALSE,
  group.dataset = NULL,
  pos.dataset = NULL,
  features.name = "features",
  only.pos = TRUE,
  features = NULL,
  return.object = TRUE,
  thresh.pc = 0,
  thresh.fc = 0,
  thresh.p = 0.05
)
cellchat_obj_w7_liver <- identifyOverExpressedGenes(cellchat_obj_w7_liver, thresh.p = 0.05, thresh.fc = 0.1, thresh.pc = 0.1) 
deg_results <- cellchat_obj_w7_liver@var.features[["features.info"]]
top15_deg <- deg_results %>%
  group_by(clusters) %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 15)
library(readr)
write_csv(top15_deg, "Top15 DEGs from w7 liver_cellchat.csv")

HSC_MPP <- c("CD34", "CD69", "SELL", "MDK", "TRH", "PRSS2", "CD200", "ANGPT1", "CD99", "NPR", "RAMP1", "LTA4H", "CXCL10", "ITGA4")

##### Stroma ####
