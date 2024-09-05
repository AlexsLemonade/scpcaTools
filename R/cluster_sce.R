#' Perform clustering on a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object that requires clustering
#' @param pc_name The name used to access the PCA results stored in the
#'   SingleCellExperiment object; default is "PCA"
#' @param BLUSPARAM A BlusterParam object specifying the clustering algorithm to
#'   use and any additional clustering options
#' @param cluster_column_name The name of the column to store the clustering
#'   results in the SCE; for naming you may want to include the type of clustering
#'   and k value of centers, e.g., "kmeans_10"
#' @param seed Seed for reproducibility of clustering results
#' @param ... Additional arguments to provide to `bluster::clusterRows()`
#'
#' @import SingleCellExperiment
#'
#' @return SingleCellExperiment object containing clustering results
#'
#' @examples
#' \dontrun{
#' # Perform K-means clustering with 10 centers
#' cluster_sce(sce, "PCA", bluster::KmeansParam(centers = 10), "kmeans_10")
#' }
#'
cluster_sce <- function(
    sce,
    pc_name = "PCA",
    BLUSPARAM,
    cluster_column_name,
    seed = NULL,
    ...) {
  # Set the seed
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Check that provided sce is indeed an SCE
  if (!is(sce, "SingleCellExperiment")) {
    stop("`sce` must be a `SingleCellExperiment` object")
  }

  # Check that the PCs are present in the SingleCellExperiment object
  if (!pc_name %in% reducedDimNames(sce)) {
    stop("The provided `pc_name` cannot be found in the reduced dimensions of the SingleCellExperiment object.")
  }

  # Extract PCs
  pcs <- reducedDim(sce, pc_name)
  clustering_results <- bluster::clusterRows(pcs, BLUSPARAM, ...)
  sce[[cluster_column_name]] <- clustering_results

  return(sce)
}
