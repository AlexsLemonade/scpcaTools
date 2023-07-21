#' Perform clustering on a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object that requires clustering
#' @param pca_name The name used to access the PCA results
#' @param BLUSPARAM A BlusterParam object specifying the clustering algorithm to
#' use
#' @param cluster_column_name The name of the column to store the clustering
#' results in the SCE
#' @param seed Seed for reproducibility of clustering results
#' @param ... Allows for additional arguments to be provided to `clusterRows()`
#'
#' @return SingleCellExperiment object containing clustering results
cluster_sce <- function(sce,
                        pca_name,
                        BLUSPARAM,
                        cluster_column_name,
                        seed = NULL,
                        ...) {
  # Set the seed
  set.seed(seed)

  # Check that the PCs are present in the SingleCellExperiment object
  if (!pca_name %in% reducedDimNames(sce)) {
    stop("The provided `pca_name` cannot be found in the reduced dimensions of
         the SingleCellExperiment object.")
  }

  # Extract PCs
  pcs <- reducedDim(sce, pca_name)
  clustering_results <- bluster::clusterRows(pcs, BLUSPARAM, ...)
  sce[[cluster_column_name]] <- clustering_results

  return(sce)
}
