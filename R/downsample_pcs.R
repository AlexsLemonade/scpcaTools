#' Downsample PCs for use in integration metric calculations
#'
#' @param pcs The PCs to downsample
#' @param frac_cells The fraction of cells to downsample to
#'
#' @return List with two items: `pcs`, the downsampled PCs; `batch_labels`, the
#'  corresponding batch labels as integers for the downsampled PCs
downsample_pcs <- function(pcs, frac_cells) {
  num_cells <- nrow(pcs)

  # Check that frac_cells is in range
  if (frac_cells < 0 | frac_cells > 1) {
    stop("The fraction of cells to downsample should be between 0 and 1.")
  }

  # Determines rows to sample
  downsampled_indices <- sample(1:num_cells,
    frac_cells * num_cells,
    replace = FALSE
  )

  # Extract PCs for downsample
  downsampled_pcs <- pcs[downsampled_indices, , drop = FALSE]

  return(
    list(
      pc_name = downsampled_pcs,
      batch_labels = rownames(downsampled_pcs)
    )
  )
}
