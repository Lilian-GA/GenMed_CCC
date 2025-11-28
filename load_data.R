library(Seurat)
install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratData)

obj <- readRDS("data/HT224P1-S1.rds")
DefaultAssay(obj)
DefaultAssay(obj) <- "RNA"
obj[["SCT"]] <- NULL

SaveH5Seurat(obj, filename = "data/HT224P1-S1.h5seurat")
Convert("data/HT224P1-S1.h5seurat", dest = "h5ad")


obj <- readRDS("data/HT224P1-S1Fc2U1Z1Bs1-SeuratObj.rds")
SaveH5Seurat(obj, filename = "data/HT224P1-S1Fc2U1Z1Bs1-SeuratObj.h5seurat")
Convert("data/HT224P1-S1Fc2U1Z1Bs1-SeuratObj.h5seurat", dest = "h5ad")


obj <- readRDS("data/HT425B1-S1H2Fs1U1Bp1-SeuratObj.rds")
SaveH5Seurat(obj, filename = "data/HT425B1-S1H2Fs1U1Bp1-SeuratObj.h5seurat")
Convert("data/HT425B1-S1H2Fs1U1Bp1-SeuratObj.h5seurat", dest = "h5ad")


obj <- readRDS("data/snRNA_L4__HT425B1-S1H1_combo.rds")
DefaultAssay(obj)
DefaultAssay(obj) <- "RNA"
obj[["SCT"]] <- NULL
SaveH5Seurat(obj, filename = "data/snRNA_L4__HT425B1-S1H1_combo.h5seurat")
Convert("data/snRNA_L4__HT425B1-S1H1_combo.h5seurat", dest = "h5ad")
