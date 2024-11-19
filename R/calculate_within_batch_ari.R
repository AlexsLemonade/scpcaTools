#' Calculate within-batch ARI scores from a merged SCE object.
#'
#'
#' @param individual_sce_list A named list of individual SCE objects. It is
#'  assumed these have a reduced dimension slot with principal components named "PCA".
#'  Column names of each object should match the contents of the `cell_id_column` for the `merged_sce.`
#' @param pc_names A list of names of the PCs in the merged SCE `reducedDim` slot
#'  object. Example: c("PCA", "fastMNN_PCA"). A within-batch ARI will be returned for each `pc_name`.
#' @param merged_sce The merged SCE object containing data from multiple batches
#' @param batch_column The variable in `merged_sce` indicating the grouping of interest.
#'  Generally this is either batches or cell types. Default is "library_id".
#' @param cell_id_column The variable in `merged_sce` indicating the original cell barcode.
#'  Default is "cell_id".
#' @param BLUSPARAM A BlusterParam object specifying the clustering algorithm to
#'   use and any additional clustering options. Default is
#'   `bluster::NNGraphParam(cluster.fun = "louvain", type = "jaccard")`
#' @param seed Seed for initializing random sampling
#' @param ... Additional arguments to provide to `bluster::clusterRows()` within `cluster_sce()`
#'
#' @return Tibble with three columns: `ari`, representing the calculated ARI values;
#'   `batch_id`, the batch id associated with each `ari`; `pc_name`, the name
#'   associated with the pc results
#'
#' @import SingleCellExperiment
#' @importFrom rlang .data
#'
#' @export
calculate_within_batch_ari <- function(
    individual_sce_list,
    pc_names,
    merged_sce,
    batch_column = "library_id",
    cell_id_column = "cell_id",
    BLUSPARAM = bluster::NNGraphParam(cluster.fun = "louvain", type = "jaccard"),
    seed = NULL,
    ...) {
  # Set the seed for subsampling and clustering
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Check that `batch_column` is in colData of SCE
  if (!batch_column %in% colnames(colData(merged_sce))) {
    stop("The specified batch column is missing from the colData of the `merged_sce`.")
  }

  # Check that `cell_id_column` is in colData of SCE
  if (!cell_id_column %in% colnames(colData(merged_sce))) {
    stop("The specified cell id column is missing from the colData of the `merged_sce`.")
  }


  # Check that `pc_name` is in merged SCE object
  if (!all(pc_names %in% reducedDimNames(merged_sce))) {
    stop("One or more of the PC names provided in `pc_names` cannot be found in the `merged_sce`.")
  }

  # Calculate within-batch ARI values across list of PCs
  within_batch_ari_tibble <- pc_names |>
    purrr::map(\(pcs) {
      within_batch_ari_from_pcs(
        individual_sce_list = individual_sce_list,
        merged_sce = merged_sce,
        batch_column = batch_column,
        cell_id_column = cell_id_column,
        pc_name = pcs,
        BLUSPARAM = BLUSPARAM
      )
    }) |>
    dplyr::bind_rows()

  return(within_batch_ari_tibble)
}

#' Calculate within-batch ARI using provided PCs for use in integration metric calculations
#'
#' @param individual_sce_list A named list of individual SCE objects. It is
#'  assumed these have a reduced dimension slot with principal components named "PCA".
#'  Column names of each object should match the contents of the `cell_id_column` for the `merged_sce.`
#' @param merged_sce The merged SCE object containing data from multiple batches
#' @param pc_name The name used to access the PCA results stored in the
#'   SingleCellExperiment object
#' @param batch_column The variable in `merged_sce` indicating the grouping of interest.
#'  Generally this is either batches or cell types. Default is "library_id".
#' @param cell_id_column The variable in `merged_sce` indicating the original cell barcode.
#'  Default is "cell_id".
#' @param BLUSPARAM A BlusterParam object specifying the clustering algorithm to
#'   use and any additional clustering options. Default is
#'   `bluster::NNGraphParam(cluster.fun = "louvain", type = "jaccard")`
#' @param ... Additional arguments to provide to `bluster::clusterRows()` within `cluster_sce()`
#'
#' @return Tibble with three columns: `ari`, representing the calculated ARI values;
#'   `batch_id`, the batch id associated with each `ari`; `pc_name`, the name
#'   associated with the pc results
within_batch_ari_from_pcs <- function(
    individual_sce_list,
    merged_sce,
    pc_name,
    batch_column = "library_id",
    cell_id_column = "cell_id",
    BLUSPARAM = bluster::NNGraphParam(cluster.fun = "louvain", type = "jaccard"),
    ...) {
  # Check that `pc_name` is in merged SCE object
  if (!(pc_name %in% reducedDimNames(merged_sce))) {
    stop("The PC name provided in `pc_names` cannot be found in the `merged_sce`.")
  }

  # Pull out the PCs or analogous reduction from merged object
  merged_pcs <- reducedDim(merged_sce, pc_name)

  # Cluster merged pcs only one time
  merged_sce <- merged_sce |>
    cluster_sce(
      pc_name = pc_name,
      BLUSPARAM = BLUSPARAM,
      cluster_column_name = "merged_clusters"
    )

  # Check that list of SCE objects is named
  batch_ids <- names(individual_sce_list)
  if (is.null(batch_ids)) {
    stop("Must provide a named list of SCE objects.")
  } else {
    # Make sure the batch ids provided match between the list and the merged object
    if (!(identical(
      sort(batch_ids),
      sort(unique(colData(merged_sce)[, batch_column]))
    ))
    ) {
      stop("Names of provided SCE objects included in the individual SCE object
           do not match batch IDs present in the batch_column of the merged object")
    }
  }

  # For every batch id, cluster and then calculate ARI for that batch
  all_ari <-
    purrr::map_dbl(
      batch_ids,
      \(batch) {
        # Cluster pc matrix for specified batch
        individual_clustering_result <-
          individual_sce_list[[batch]] |>
          cluster_sce(
            pc_name = "PCA",
            BLUSPARAM = BLUSPARAM,
            cluster_column_name = "individual_clusters"
          )

        # Extract cells from merged clustering for batch
        merged_clusters <- merged_sce$merged_clusters |>
          purrr::set_names(colData(merged_sce)[[cell_id_column]])
        cells_to_keep <- which(merged_sce[[batch_column]] == batch)
        batch_merged_clusters <- merged_clusters[cells_to_keep]

        # Check order of cells
        if (!all.equal(names(batch_merged_clusters), colnames(individual_clustering_result))) {
          stop(
            "The cell barcodes and/or order of the merged and individual clustering results are not identical."
          )
        }

        # Calculate ARI between pre-merged clustering and post-merged clustering for the given batch
        ari <-
          bluster::pairwiseRand(
            individual_clustering_result$individual_clusters,
            batch_merged_clusters,
            mode = "index"
          )

        return(ari)
      }
    )

  # Create tibble with ARI and batch id
  within_batch_ari_tibble <- tibble::tibble(
    ari = all_ari,
    batch_id = batch_ids,
    pc_name = pc_name
  )

  return(within_batch_ari_tibble)
}
