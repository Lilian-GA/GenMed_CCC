library(Seurat)
library(SeuratDisk)
library(reticulate)

# Charger l'objet
obj <- readRDS("data/data_rds/HT259P1-S1H1Fc2U1Z1Bs1-SeuratObj.rds")
DefaultAssay(obj) <- "Spatial"

# Créer un nouvel objet Seurat minimal pour Scanpy
new_obj <- CreateSeuratObject(
  counts = obj@assays$Spatial@counts,
  meta.data = obj@meta.data
)

# Ajouter un assay Spatial minimal
new_obj[["Spatial"]] <- CreateAssayObject(counts = obj@assays$Spatial@counts)
DefaultAssay(new_obj) <- "Spatial"

# Ajouter coordonnées spatial dans un DimReduc (nom libre)
coords_mat <- as.matrix(obj@images$slice1@coordinates[, c("imagerow", "imagecol")])
colnames(coords_mat) <- c("spatial_1", "spatial_2")

new_obj[["spatial"]] <- CreateDimReducObject(
  embeddings = coords_mat,
  key = "spatial_",
  assay = "Spatial"
)

# Ajouter PCA/UMAP si existant
if (!is.null(obj@reductions$pca)) new_obj[["pca"]] <- obj@reductions$pca
if (!is.null(obj@reductions$umap)) new_obj[["umap"]] <- obj@reductions$umap

# Exporter h5ad directement via SeuratDisk (sans image)
SaveH5Seurat(new_obj, filename = "data/HT259P1-S1H1Fc2U1Z1Bs1-SeuratObj.h5seurat", overwrite = TRUE)
Convert("data/HT259P1-S1H1Fc2U1Z1Bs1-SeuratObj.h5seurat", dest = "h5ad")
