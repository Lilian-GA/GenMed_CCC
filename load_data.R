library(Seurat)
library(SeuratDisk)

#Pancreas : HT224P1
obj <- readRDS("data/HT224P1-S1.rds")
DefaultAssay(obj)
DefaultAssay(obj) <- "RNA"
obj[["SCT"]] <- NULL
SaveH5Seurat(obj, filename = "data/HT224P1-S1.h5seurat")
Convert("data/HT224P1-S1.h5seurat", dest = "h5ad")


convert_seurat_to_h5ad <- function(
  rds_path,
  h5seurat_path
) {

  message("Lecture du fichier RDS…")
  obj <- readRDS(rds_path)

  message("Assay par défaut actuel : ", DefaultAssay(obj))
  DefaultAssay(obj) <- "RNA"
  message("Assay par défaut fixé à 'RNA'")

  if ("SCT" %in% names(obj@assays)) {
    obj[["SCT"]] <- NULL
    message("Assay SCT supprimé")
  } else {
    message("Aucun assay SCT à supprimer")
  }

  message("Sauvegarde en h5seurat : ", h5seurat_path)
  SaveH5Seurat(obj, filename = h5seurat_path, overwrite = TRUE)

  message("Conversion en h5ad : ", h5ad_path)
  Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)

  message("Conversion terminée ! Fichiers créés :")
  message("   - ", h5seurat_path)
}

convert_seurat_to_h5ad()