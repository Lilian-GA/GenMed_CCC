install.packages("Seurat")
library(Seurat)

library(SeuratDisk)

#Pancreas : HT224P1
obj <- readRDS("data/HT224P1-S1.rds")
DefaultAssay(obj)
DefaultAssay(obj) <- "RNA"
obj[["SCT"]] <- NULL
SaveH5Seurat(obj, filename = "data/HT224P1-S1.h5seurat")
Convert("data/HT224P1-S1.h5seurat", dest = "h5ad")



obj <- readRDS("data/snRNA_L4__HT425B1-S1H1_combo.rds")
DefaultAssay(obj)
DefaultAssay(obj) <- "RNA"
obj[["SCT"]] <- NULL
SaveH5Seurat(obj, filename = "data/snRNA_L4__HT425B1-S1H1_combo.h5seurat")
Convert("data/snRNA_L4__HT425B1-S1H1_combo.h5seurat", dest = "h5ad")

obj <- readRDS("data/data_rds/HT259P1-S1H1.rds")
DefaultAssay(obj)
DefaultAssay(obj) <- "RNA"
obj[["SCT"]] <- NULL
SaveH5Seurat(obj, filename = "data/HT259P1-S1H1.h5seurat")
Convert("data/HT259P1-S1H1.h5seurat", dest = "h5ad")

obj <- readRDS("data/data_rds/snRNA_L4__HT243B1-H3A2.rds")
DefaultAssay(obj)
DefaultAssay(obj) <- "RNA"
obj[["SCT"]] <- NULL
SaveH5Seurat(obj, filename = "data/snRNA_L4__HT243B1-H3A2.h5seurat")
Convert("data/snRNA_L4__HT243B1-H3A2.h5seurat", dest = "h5ad")

