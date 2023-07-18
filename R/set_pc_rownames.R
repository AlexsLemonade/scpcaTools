#' Label rownames of provided PCs and remove cells with NA batch label from PCs
#'
#' @param pcs The PCs with unlabeled rownames
#' @param batches Vector of cell-wise batch information whose length is the number of
#'   rows in `pcs`
#'
#' @return PCs with labeled rownames and NA batch cells removed
set_pc_rownames <- function(pcs, batches) {
  # Check dimensions
  if (nrow(pcs) != length(batches)) {
    stop("Incompatable PC and batch information dimensions after removing NAs.")
  }

  # Remove NA batch cells
  retain_indices <- which(!is.na(batches))
  batches <- batches[retain_indices]
  pcs <- pcs[retain_indices, ]

  # Set PC rownames to be the batches
  rownames(pcs) <- batches

  return(pcs)
}
