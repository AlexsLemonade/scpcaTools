#' Downsample PCs for use in integration metric calculations 
#'  
#' @param integrated_pcs The full set of integrated pcs
#' @param pc_name The name that allows access to the PCs. Example: fastMNN_PCA
#' @param frac_cells The fraction of cells to downsample to
#'
#' @return List with two items: `pcs`, the downsampled PCs; `batch_labels`, the 
#'  corresponding batch labels as integers for the downsampled PCs
downsample_for_metrics <- function(integrated_pcs, pc_name, frac_cells) {
  
  num_cells <- nrow(integrated_pcs)
  
  # Determines rows to sample
  downsampled_indices <- sample(1:num_cells,
                                frac_cells*num_cells,
                                replace = FALSE) 
  
  # Extract PCs for downsample
  downsampled_integrated_pcs <- integrated_pcs[downsampled_indices, , drop = FALSE]
  
  return (
    list(
      pc_name = downsampled_integrated_pcs,
      batch_labels = rownames(downsampled_integrated_pcs)
    )
  )
  
}
