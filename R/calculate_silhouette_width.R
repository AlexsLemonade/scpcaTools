#' Calculate silhouette width scores from an integrated SCE object
#'
#' @param integrated_sce The integrated SCE object
#' @param pc_name The name that allows access to the PCs. Example: fastMNN_PCA
#' @param frac_cells The fraction of cells to downsample to. Default: 0.8
#' @param nreps The number of times to repeat sub-sampling procedure. Default: 20
#' @param seed Seed for initializing random sampling
#' @param batch_column The variable in `integrated_sce` indicating the grouping of interest.
#'  Generally this is either batches or cell types. Default is "library_id".
#'
#' @return Tibble with six columns: `rep`, representing the given downsampling replicate;
#'   `silhouette_width`, the calculated silhouette width for the given `rep`; `silhouette_cluster`,
#'   the assigned cluster for the cell during silhouette calculation, i.e. the true identity;
#'   `other_cluster`, the other assigned for the cell during silhouette calculation;
#'   `pc_name`, the name associated with the pc results
#'
#' @import SingleCellExperiment
#' @import bluster
#'
#' @export
calculate_silhouette_width <- function(integrated_sce,
                                       pc_name,
                                       frac_cells = 0.8,
                                       nreps = 20,
                                       seed = 2023,
                                       batch_column = "library_id") {

  # Set the seed for subsampling
  set.seed(seed)

  # Check that provided `pc_name` is present in SingleCellExperiment object
  if(!pc_name %in% reducedDimNames(integrated_sce)) {
    stop("The provided `pc_name` cannot be found in the SingleCellExperiment object.")
  }

  # Check that `nreps` is an integer
  if(!(is.integer(nreps))) {
    stop("The provided `nreps` should be an integer.")
  }

  # Check that frac_cells is in range
  if(frac_cells < 0 | frac_cells > 1) {
    stop("The fraction of cells to downsample should be between 0 and 1.")
  }

  # Pull out the PCs or analogous reduction
  pcs <- reducedDim(integrated_sce, pc_name)


  # Perform calculations
  all_silhouette <- purrr::map(1:nreps, \(rep) {

    # Downsample PCs
    downsampled <- downsample_pcs(pcs, frac_cells)
    # Calculate batch ASW and add into final tibble
    rep_silhouette <- bluster::approxSilhouette(downsampled$pc_name, downsampled$batch_labels) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(rep = rep) %>%
      dplyr::select(rep,
                    silhouette_width = width,
                    silhouette_cluster = cluster,
                    other_cluster = other) |>
      dplyr::bind_rows()

  })

  # Add integration method into tibble
  all_silhouette <- dplyr::mutate(all_silhouette, pc_name = pc_name)

  # Return tibble with silhouette width results which can further be summarized downstream
  return(all_silhouette)

}
