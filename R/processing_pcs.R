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
  if(length(retain_indices) == 0) {
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
#'
#' @return List with two items: `pcs`, the downsampled PCs; `batch_labels`, the
#'  corresponding batch labels for the downsampled PCs
downsample_pcs <- function(pcs, frac_cells) {

  # Check that there is a minimum number of cells
  num_cells <- nrow(pcs)
  if (!num_cells > 1) {
    stop("There are not enough cells to perform downsampling.")
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
