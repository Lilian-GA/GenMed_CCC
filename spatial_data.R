library(Seurat)
library(SeuratDisk)

convert_seurat_spatial_to_h5ad <- function(
  rds_path,
  spatial_image = "slice1",
  h5seurat_out
) {

  message("Chargement du fichier RDS...")
  obj <- readRDS(rds_path)
  DefaultAssay(obj) <- "Spatial"

  message("Extraction des données...")

  # Create minimal Seurat object
  new_obj <- CreateSeuratObject(
    counts = obj@assays$Spatial@counts,
    meta.data = obj@meta.data
  )

  # Add Spatial assay
  new_obj[["Spatial"]] <- CreateAssayObject(
    counts = obj@assays$Spatial@counts
  )
  DefaultAssay(new_obj) <- "Spatial"

  # Add spatial coordinates into DimReduc
  if (!spatial_image %in% names(obj@images)) {
    stop(paste0("L'image '", spatial_image, "' n'existe pas dans obj@images."))
  }

  coords <- obj@images[[spatial_image]]@coordinates[, c("imagerow", "imagecol")]
  coords_mat <- as.matrix(coords)
  colnames(coords_mat) <- c("spatial_1", "spatial_2")

  new_obj[["spatial"]] <- CreateDimReducObject(
    embeddings = coords_mat,
    key = "spatial_",
    assay = "Spatial"
  )

  # Add PCA / UMAP if present
  if (!is.null(obj@reductions$pca)) {
    message("Ajout PCA")
    new_obj[["pca"]] <- obj@reductions$pca
  }
  if (!is.null(obj@reductions$umap)) {
    message("Ajout UMAP")
    new_obj[["umap"]] <- obj@reductions$umap
  }

  # Save h5seurat
  message("Enregistrement en h5seurat : ", h5seurat_out)
  SaveH5Seurat(new_obj, filename = h5seurat_out, overwrite = TRUE)

  # Convert to h5ad
  message("Conversion en h5ad : ", h5seurat_out)
  Convert(h5seurat_out, dest = "h5ad", overwrite = TRUE)

  message("Conversion terminée avec succès !")
}


obj <- readRDS("data/data_rds/HT224P1-S1Fc2U1Z1Bs1-SeuratObj.rds")

DefaultAssay(obj) <- "Spatial"
png::writePNG(obj@images$slice1@image, "histology_pancreas_224.png")
