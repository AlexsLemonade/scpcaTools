#' Perform clustering on a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object that requires clustering
#' @param seed Seed for reproducibility of clustering results
#' @param pca_column The SingleCellExperiment column name to reference the PCA results
#' @param k_range Range of k values to use for k-means clustering. Default: `seq(5, 25, 5)`
#'
#' @return SingleCellExperiment object containing clustering results
cluster_sce <- function(sce,
                        seed = 2023,
                        pca_column,
                        k_range) {

  for (k in k_range) {

    # Set cluster name
    cluster_name <- paste("kmeans", k, sep = "_")

    # Perform k-means clustering
    clustering_results <- kmeans(x = reducedDim(sce, pca_column), centers = k)
    sce[[cluster_name]] <- unname(clustering_results$cluster)

  }

  return(sce)
}
