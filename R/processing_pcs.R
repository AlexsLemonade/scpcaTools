#' Remove cells with NA batch label from PCs and label rownames
#'
#' @param pcs The PCs to be filtered and labeled
#' @param batches Vector of cell-wise batch information whose length is the number of
#'   rows in `pcs`
#'
#' @return PCs with NA batch cells removed and labeled rownames
filter_pcs <- function(pcs, batches) {
  # Remove NA batch cells
  retain_indices <- which(!is.na(batches))

  # Check that there are retained indices present
  if (length(retain_indices) == 0) {
    stop("There are no batch cells that are not NA present.")
  }

  # Filter batches and PCs
  batches <- batches[retain_indices]
  pcs <- pcs[retain_indices, ]

  # Check dimensions after filtering
  if (nrow(pcs) != length(batches)) {
    stop("Incompatable PC and batch information dimensions after removing NAs.")
  }

  # Set PC rownames to be the batches
  rownames(pcs) <- batches

  return(pcs)
}



#' Downsample PCs for use in integration metric calculations
#'
#' @param pcs The PCs to downsample, these PCs should contain batch labels as rownames
#' @param frac_cells The fraction of cells to downsample to
#' @param min_cells The minimum number of cells after downsampling. Default: 50
#'
#' @return The downsampled PCs
downsample_pcs <- function(pcs, frac_cells, min_cells = 50) {
  # Check that there is a minimum number of cells
  num_cells <- nrow(pcs)
  if (frac_cells * num_cells < min_cells) {
    stop("Downsampling would result in fewer cells than the `min_cells` threshold.")
  }

  # Check that frac_cells is in range
  if (frac_cells <= 0 | frac_cells >= 1) {
    stop("The fraction of cells to downsample should be between 0 and 1.")
  }

  # Determines rows to sample
  downsampled_indices <- sample(1:num_cells,
    frac_cells * num_cells,
    replace = FALSE
  )

  # Extract PCs for downsample
  downsampled_pcs <- pcs[downsampled_indices, , drop = FALSE]

  return(downsampled_pcs)
}

#' Calculate within-batch ARI using provided PCs for use in integration metric calculations
#'
#' @param merged_sce The merged SCE object containing data from multiple batches
#' @param pc_name The name used to access the PCA results stored in the
#'   SingleCellExperiment object; default is "PCA"
#' @param batch_column The variable in `merged_sce` indicating the grouping of interest.
#'  Generally this is either batches or cell types. Default is "library_id".
#'
#' @return Tibble with three columns: `ari`, representing the calculated ARI values;
#'   `batch_id`, the batch id associated with each `ari`; `pc_name`, the name
#'   associated with the pc results

within_batch_ari_pcs <- function(merged_sce, pc_name, batch_column = "library_id") {
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
  within_batch_ari_tibble <- tibble::tibble(ari = all_ari,
                                            batch_id = batch_ids,
                                            pc_name = pc_name)

  return(within_batch_ari_tibble)
}
