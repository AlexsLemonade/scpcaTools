#' Calculate within-batch ARI scores from a merged SCE object.
#'
#'
#' @param individual_sce_list A named list of individual SCE objects. It is
#'  assumed these have a reduced dimension slot with principal components named "PCA".
#' @param pc_list A list of names that allow access to the PCs in the merged SCE
#'  object. Example: c("PCA", "fastMNN_PCA"). A within-batch ARI will be returned for each `pc_name`.
#' @param merged_sce The merged SCE object containing data from multiple batches
#' @param batch_column The variable in `merged_sce` indicating the grouping of interest.
#'  Generally this is either batches or cell types. Default is "library_id".
#' @param seed Seed for initializing random sampling
#'
#' @return Tibble with three columns: `ari`, representing the calculated ARI values;
#'   `batch_id`, the batch id associated with each `ari`; `pc_name`, the name
#'   associated with the pc results
#'
#' @import SingleCellExperiment
#' @importFrom rlang .data
#'
#' @export
calculate_within_batch_ari <- function(individual_sce_list,
                                       pc_list,
                                       merged_sce,
                                       batch_column = "library_id",
                                       seed = NULL) {
  # Set the seed for subsampling and clustering
  set.seed(seed)

  # Check that `batch_column` is in colData of SCE
  if (!batch_column %in% colnames(colData(merged_sce))) {
    stop("The specified batch column is missing from the colData of the `merged_sce`.")
  }

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

  # Check that `pc_name` is in merged SCE object
  if (!any(pc_list %in% reducedDimNames(merged_sce))) {
    stop("One or more of the PC names provided in `pc_list` cannot be found in the `merged_sce`.")
  }

  # Calculate within-batch ARI values across list of PCs
  within_batch_ari_tibble <-
    purrr::map(
      pc_list,
      \(pcs)
      within_batch_ari_from_pcs(
        merged_sce = merged_sce,
        batch_column = batch_column,
        pc_name = pcs
      )
    ) |>
    dplyr::bind_rows()

  return(within_batch_ari_tibble)
}

#' Calculate within-batch ARI using provided PCs for use in integration metric calculations
#'
#' @param merged_sce The merged SCE object containing data from multiple batches
#' @param pc_name The name used to access the PCA results stored in the
#'   SingleCellExperiment object
#' @param batch_column The variable in `merged_sce` indicating the grouping of interest.
#'  Generally this is either batches or cell types. Default is "library_id".
#'
#' @return Tibble with three columns: `ari`, representing the calculated ARI values;
#'   `batch_id`, the batch id associated with each `ari`; `pc_name`, the name
#'   associated with the pc results

within_batch_ari_from_pcs <-
  function(merged_sce, pc_name, batch_column = "library_id") {
    # Check that `pc_name` is in merged SCE object
    if (!pc_name %in% reducedDimNames(merged_sce)) {
      stop("The provided `pc_name` cannot be found in the `merged_sce` object.")
    }

    # Pull out the PCs or analogous reduction from merged object
    merged_pcs <- reducedDim(merged_sce, pc_name)

    # Cluster merged pcs only one time
    merged_sce <- merged_sce |>
      cluster_sce(
        pc_name = pc_name,
        BLUSPARAM = bluster::NNGraphParam(cluster.fun = "louvain", type = "jaccard"),
        cluster_column_name = "merged_clusters"
      )

    merged_clusters <- merged_sce$merged_clusters |>
      purrr::set_names(merged_sce[[batch_column]])

    # For every batch id, cluster and then calculate ARI for that batch
    all_ari <- batch_ids |>
      purrr::map_dbl(\(batch) {
        # Cluster pc matrix for specified batch
        individual_clustering_result <-
          individual_sce_list[[batch]] |>
          cluster_sce(
            pc_name = "PCA",
            BLUSPARAM = bluster::NNGraphParam(cluster.fun = "louvain", type = "jaccard"),
            cluster_column_name = "individual_clusters"
          )

        # Extract clusters from merged clustering for batch
        clusters_to_keep <- grep(batch, merged_sce[[batch_column]])
        batch_merged_clusters <- merged_clusters[clusters_to_keep]

        # Calculate ARI between pre-merged clustering and post-merged clustering for the given batch
        ari <-
          bluster::pairwiseRand(
            individual_clustering_result$individual_clusters,
            batch_merged_clusters,
            mode = "index"
          )

        return(ari)
      })

    # Create tibble with ARI and batch id
    within_batch_ari_tibble <- tibble::tibble(
      ari = all_ari,
      batch_id = batch_ids,
      pc_name = pc_name
    )

    return(within_batch_ari_tibble)
  }
