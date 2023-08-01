#' Calculate iLISI (integration Local Inverse Simpson's Index) scores
#'
#'
#' @param merged_sce The merged SCE object containing data from multiple batches
#' @param pc_name The name that allows access to the PCs. Example: fastMNN_PCA
#' @param batch_column The variable in `merged_sce` indicating the grouping of interest.
#'  Generally this is either batches or cell types. Default is "library_id".
#'
#'
#' @return Data frame with three columns, one row per cell. Columns are `ilisi_score`,
#'   `cell_barcode`, and `batch_id`
#'
#' @import SingleCellExperiment
#'
#' @export
calculate_ilisi <- function(merged_sce,
                            pc_name,
                            batch_column = "library_id") {
  # Check that provided `pc_name` is present in SingleCellExperiment object
  if (!pc_name %in% reducedDimNames(merged_sce)) {
    stop("The provided `pc_name` cannot be found in the SingleCellExperiment object.")
  }

  # Check that `batch_column` is in colData of SCE
  if (!batch_column %in% colnames(colData(merged_sce))) {
    stop("The specified batch column is missing from the colData of the SingleCellExperiment object.")
  }

  # Pull out the PCs or analogous reduction
  pcs <- reducedDim(merged_sce, pc_name)

  # Check that PCs have rownames already, since we won't be renaming them
  if (is.null(rownames(pcs))) {
    stop("PCs are missing rownames for iLISI calculation.")
  }

  # Remove batch NAs from PCs and label rownames
  labeled_pcs <- filter_pcs(pcs,
    colData(merged_sce)[, batch_column],
    # don't rename with batch labels
    rename_pcs = FALSE
  )

  # Create data frame for input to lisi calculation
  batch_df <- data.frame(batch = rownames(labeled_pcs))

  # Calculate iLISI score
  ilisi_df <- lisi::compute_lisi(
    labeled_pcs,
    # define the batches
    batch_df,
    # which variables in `batch_df` to compute lisi for
    "batch"
  ) |>
    dplyr::rename(ilisi_score = batch) |>
    dplyr::mutate(
      cell_barcode = rownames(labeled_pcs),
      batch_id = batch_df$batch
    )

  # Return data frame with cell-wise iLISI scores and associated cell & batch identifiers
  return(ilisi_df)
}
