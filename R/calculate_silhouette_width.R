#' Calculate silhouette width scores from an integrated SCE object.
#' This function performs replicated calculations on downsampled data.
#'
#' @param integrated_sce The integrated SCE object
#' @param pc_name The name that allows access to the PCs. Example: fastMNN_PCA
#' @param batch_column The variable in `integrated_sce` indicating the grouping of interest.
#'  Generally this is either batches or cell types. Default is "library_id".
#' @param frac_cells The fraction of cells to downsample to. Default: 0.8
#' @param nreps The number of times to repeat sub-sampling procedure. Default: 20
#' @param seed Seed for initializing random sampling
#'
#' @return Tibble with five columns: `rep`, representing the given downsampling replicate;
#'   `silhouette_width`, the calculated silhouette width for the given `rep`; `silhouette_cluster`,
#'   the assigned cluster for the cell during silhouette calculation, i.e. the true identity;
#'   `other_cluster`, the other assigned for the cell during silhouette calculation;
#'   `pc_name`, the name associated with the pc results
#'
#' @import SingleCellExperiment
#'
#' @export
calculate_silhouette_width <- function(integrated_sce,
                                       pc_name,
                                       batch_column = "library_id",
                                       frac_cells = 0.8,
                                       nreps = 20,
                                       seed = NULL) {
  # Set the seed for subsampling
  set.seed(seed)

  # Check that provided `pc_name` is present in SingleCellExperiment object
  if (!pc_name %in% reducedDimNames(integrated_sce)) {
    stop("The provided `pc_name` cannot be found in the SingleCellExperiment object.")
  }

  # Check that `nreps` is an integer
  if (!(nreps == round(nreps))) {
    stop("The provided `nreps` should be an integer.")
  }

  # Check that frac_cells is in range
  if (frac_cells <= 0 | frac_cells >= 1) {
    stop("The fraction of cells to downsample should be between 0 and 1.")
  }

  # Pull out the PCs or analogous reduction
  pcs <- reducedDim(integrated_sce, pc_name)

  # Check that `batch_column` is in colData of SCE
  if (!batch_column %in% colnames(colData(integrated_sce))) {
    stop("The specified batch column is missing from the colData of the SingleCellExperiment object.")
  }

  # Remove batch NAs from PCs and label rownames
  labeled_pcs <- filter_pcs(pcs, colData(integrated_sce)[, batch_column])

  # Perform calculations
  all_silhouette <- purrr::map(1:nreps, \(rep) {
    # Downsample PCs
    downsampled <- downsample_pcs(labeled_pcs, frac_cells)
    # Calculate batch ASW and add into final tibble
    bluster::approxSilhouette(downsampled, rownames(downsampled)) |>
      tibble::as_tibble() |>
      dplyr::mutate(rep = rep) |>
      dplyr::select(
        rep,
        silhouette_width = width,
        silhouette_cluster = cluster,
        other_cluster = other
      )
  }) |>
    dplyr::bind_rows() |>
    # Add integration method into tibble
    dplyr::mutate(pc_name = pc_name)

  # Return tibble with silhouette width results which can further be summarized downstream
  return(all_silhouette)
}
