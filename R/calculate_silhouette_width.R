#' Calculate silhouette width scores from a merged SCE object.
#'
#' This function performs replicated calculations on downsampled data.
#'
#' @param merged_sce The merged SCE object containing data from multiple batches
#' @param pc_names A list of names that allow access to the PCs in the merged SCE
#'  object. Example: c("PCA", "fastMNN_PCA").
#' @param batch_column The variable in `merged_sce` indicating the grouping of interest.
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
#' @importFrom rlang .data
#'
#' @export
calculate_silhouette_width <- function(merged_sce,
                                       pc_names,
                                       batch_column = "library_id",
                                       frac_cells = 0.8,
                                       nreps = 20,
                                       seed = NULL) {
  # Check that provided `pc_name` is present in SingleCellExperiment object
  if (!all(pc_names %in% reducedDimNames(merged_sce))) {
    stop("One or more of the PC names provided in `pc_names` cannot be found in the `merged_sce`.")
  }

  # Check that `nreps` is an integer
  if (!(nreps == round(nreps))) {
    stop("The provided `nreps` should be an integer.")
  }

  # Check that frac_cells is in range
  if (frac_cells <= 0 || frac_cells >= 1) {
    stop("The fraction of cells to downsample should be between 0 and 1.")
  }

  # Check that `batch_column` is in colData of SCE
  if (!batch_column %in% colnames(colData(merged_sce))) {
    stop("The specified batch column is missing from the colData of the SingleCellExperiment object.")
  }

  # Calculate the silhouette width values across list of PCs
  all_silhouette_df <- purrr::map(
    pc_names,
    \(pcs) {
      silhouette_width_from_pcs(
        merged_sce = merged_sce,
        batch_column = batch_column,
        pc_name = pcs,
        frac_cells = frac_cells,
        nreps = nreps,
        seed = seed
      )
    }
  ) |>
    dplyr::bind_rows()

  # Return tibble with silhouette width results which can further be summarized downstream
  return(all_silhouette_df)
}


#' Calculate silhouette width scores using PCs from a merged SCE object for use
#' in integration metric calculations
#'
#' @param merged_sce The merged SCE object containing data from multiple batches
#' @param pc_name The name used to access the PCA results stored in the
#'   SingleCellExperiment object
#' @param batch_column The variable in `merged_sce` indicating the grouping of interest.
#'  Generally this is either batches or cell types. Default is "library_id".
#' @param frac_cells The fraction of cells to downsample to. Default: 0.8
#' @param nreps The number of times to repeat sub-sampling procedure. Default: 20
#' @param seed Seed for initializing random sampling. Default is NULL.
#'
#' @return Tibble with five columns: `rep`, representing the given downsampling replicate;
#'   `silhouette_width`, the calculated silhouette width for the given `rep`; `silhouette_cluster`,
#'   the assigned cluster for the cell during silhouette calculation, i.e. the true identity;
#'   `other_cluster`, the other assigned for the cell during silhouette calculation;
#'   `pc_name`, the name associated with the pc results
silhouette_width_from_pcs <-
  function(merged_sce,
           pc_name,
           batch_column = "library_id",
           frac_cells = 0.8,
           nreps = 20,
           seed = NULL) {
    # Set the seed for subsampling
    if (!is.null(seed)) {
      set.seed(seed)
    }

    # Pull out the PCs or analogous reduction
    pcs <- reducedDim(merged_sce, pc_name)

    # Remove batch NAs from PCs and label rownames
    labeled_pcs <-
      filter_pcs(pcs, colData(merged_sce)[[batch_column]])

    # Perform calculations
    all_silhouette <- purrr::map(1:nreps, \(rep) {
      # Downsample PCs
      downsampled <- downsample_pcs(labeled_pcs, frac_cells)
      # Calculate batch ASW and add into final tibble
      bluster::approxSilhouette(downsampled, rownames(downsampled)) |>
        tibble::as_tibble() |>
        dplyr::mutate(rep = rep) |>
        dplyr::select(
          "rep",
          silhouette_width = "width",
          silhouette_cluster = "cluster",
          other_cluster = "other"
        )
    }) |>
      dplyr::bind_rows() |>
      # Add integration method into tibble
      dplyr::mutate(pc_name = pc_name)

    return(all_silhouette)
  }
