#' Perform clustering on a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object that requires clustering
#' @param pca_name The name used to access the PCA results
#' @param BLUSPARAM A BlusterParam object specifying the clustering algorithm to
#' use - can be `KmeansParam` or `NNGraphParam`
#' @param cluster_name The type of clustering method to be performed - can be
#' "kmeans", "walktrap", or "louvain"
#' @param seed Seed for reproducibility of clustering results
#'
#' @return SingleCellExperiment object containing clustering results
cluster_sce <- function(sce,
                        pca_name,
                        BLUSPARAM,
                        cluster_name,
                        seed = NULL,
                        ...) {
  # Set the seed
  set.seed(seed)

  # Extract PCs
  pcs <- reducedDim(sce, pca_name)
  clustering_results <- bluster::clusterRows(pcs, BLUSPARAM, ...)
  sce[[cluster_name]] <- clustering_results

  return(sce)
}
