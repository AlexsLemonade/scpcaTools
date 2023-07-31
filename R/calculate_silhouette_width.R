#' Calculate silhouette width scores from a merged SCE object.
#'
#' This function performs replicated calculations on downsampled data.
#'
#' @param merged_sce The merged SCE object containing data from multiple batches
#' @param pc_list A list of names that allow access to the PCs in the merged SCE
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
                                       pc_list,
                                       batch_column = "library_id",
                                       frac_cells = 0.8,
                                       nreps = 20,
                                       seed = NULL) {
  # Set the seed for subsampling
  set.seed(seed)

  # Check that provided `pc_name` is present in SingleCellExperiment object
  if (!any(pc_list %in% reducedDimNames(merged_sce))) {
    stop("One or more of the PC names provided in `pc_list` cannot be found in the `merged_sce`.")
  }

  # Check that `nreps` is an integer
  if (!(nreps == round(nreps))) {
    stop("The provided `nreps` should be an integer.")
  }

  # Check that frac_cells is in range
  if (frac_cells <= 0 | frac_cells >= 1) {
    stop("The fraction of cells to downsample should be between 0 and 1.")
  }

  # Check that `batch_column` is in colData of SCE
  if (!batch_column %in% colnames(colData(merged_sce))) {
    stop("The specified batch column is missing from the colData of the SingleCellExperiment object.")
  }

  # Calculate the silhouette width values across list of PCs
  all_silhouette_df <- purrr::map_df(
    pc_list,
    ~ calculate_silhouette_width_pcs(
      merged_sce = merged_sce,
      batch_column = batch_column,
      pc_name = .
    )
  )

  # Return tibble with silhouette width results which can further be summarized downstream
  return(all_silhouette_df)
}
