#' Perform clustering on a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object that requires clustering
#' @param seed Seed for reproducibility of clustering results
#' @param pca_column The SingleCellExperiment column name to reference the PCA results
#' @param cluster_type The type of clustering method to be performed - can be "kmeans", "walktrap", or "louvain"
#' @param k The desired number of centers
#'
#' @return SingleCellExperiment object containing clustering results
cluster_sce <- function(sce,
                        pca_column,
                        cluster_type,
                        seed = NULL,
                        k = NULL,
                        ...) {
  # Set the seed
  set.seed(seed)

  # Extract PCs
  pcs <- reducedDim(sce, pca_column)
  if (cluster_type == "kmeans") {
    # Perform k-means clustering
    clustering_results <- bluster::clusterRows(pcs, bluster::KmeansParam(centers = k))
  } else if (cluster_type %in% c("walktrap", "louvain")) {
    # Grab weighting type based on cluster type
    weighting_type <- ifelse(cluster_type == "louvain", "jaccard", "rank")
    # Perform graph-based clustering
    clustering_results <- bluster::clusterRows(
      pcs,
      bluster::NNGraphParam(
        k = k,
        type = weighting_type,
        cluster.fun = cluster_type,
      )
    )
  } else {
    stop("Clustering type is not valid. Please specify 'kmeans', 'walktrap' or 'louvain' clustering.")
  }

  # Set column name to store results in
  cluster_name <- paste(cluster_type, k, sep = "_")
  sce[[cluster_name]] <- clustering_results

  return(sce)
}
