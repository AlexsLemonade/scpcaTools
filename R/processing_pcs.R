#' Remove cells with NA batch label from PCs and label rownames
#'
#' @param pcs The PCs to be filtered and labeled
#' @param batches Vector of cell-wise batch information whose length is the number of
#'   rows in `pcs`
#' @param rename_pcs Whether to rename the PC rownames with the cell's batch label.
#'   Default is `TRUE`.
#'
#' @return PCs with NA batch cells removed and labeled rownames
filter_pcs <- function(pcs, batches, rename_pcs = TRUE) {
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

  # rename PCs to batches if specified
  if (rename_pcs) {
    rownames(pcs) <- batches
  }

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
