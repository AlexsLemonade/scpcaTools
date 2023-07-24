#' Calculate within-batch ARI scores from an integrated SCE object.
#'
#'
#' @param individual_sce_list A list of individual SCE objects
#' @param pc_name The name that allows access to the PCs. Example: fastMNN_PCA
#' @param integrated_sce The integrated SCE object
#' @param batch_column The variable in `integrated_sce` indicating the grouping of interest.
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
                                       pc_name,
                                       integrated_sce,
                                       batch_column = "library_id",
                                       seed = NULL) {
  # Set the seed for subsampling and clustering
  set.seed(seed)

  # Check that list of SCE objects is named
  batch_ids <- names(individual_sce_list)
  if(is.null(batch_ids)){
    stop("Must provide a named list of SCE objects.")
  } else {
    # Make sure the batch ids provided match between the list and the integrated object
    if(!(identical(
      sort(batch_ids),
      sort(unique(colData(integrated_sce)[, batch_column]))))
    ){
      stop("Names of provided SCE objects included in the individual SCE object
           do not match batch IDs present in the batch_column of the integrated object")
    }
  }

  # Pull out the PCs or analogous reduction from integrated object
  integrated_pcs <- reducedDim(integrated_sce, pc_name)

  # Cluster integrated pcs only one time
  integrated_clustering_result <- integrated_sce |>
    cluster_sce(
      pca_name = pc_name,
      BLUSPARAM = bluster::NNGraphParam(cluster.fun = "louvain", type = "jaccard"),
      cluster_column_name = "integrated_clusters"
    )

  integrated_clustering_result <- colData(integrated_sce)[["integrated_clusters"]] |>
    rlang::set_names(rownames(integrated_pcs)) # make sure to set the names with the batch ids


  # For every batch id, cluster and then calculate ARI for that batch
  all_ari <- batch_ids |>
    purrr::map_dbl(\(batch){

      # Cluster pc matrix for specified batch
      individual_clustering_result <- individual_sce_list[[batch]] |>
        cluster_sce(
          pca_name = pc_name,
          BLUSPARAM = bluster::NNGraphParam(cluster.fun = "louvain", type = "jaccard"),
          cluster_column_name = "individual_clusters"
        )

      # Extract clusters from integrated clustering for batch
      clusters_to_keep <- grep(batch, names(integrated_clustering_result))
      batch_integrated_clusters <- integrated_clustering_result[clusters_to_keep]

      # Calculate ARI between pre-integrated clustering and post-integrated clustering for the given batch
      ari <- bluster::pairwiseRand(individual_clustering_result,
                                   batch_integrated_clusters,
                                   mode = "index")

      return(ari)

    })

  # Create tibble with ARI and batch id
  within_batch_ari_tibble <- tibble::tibble(
    ari = all_ari,
    batch_id = batch_ids,
    pca_name = pca_name
  )

  return(within_batch_ari_tibble)
}
