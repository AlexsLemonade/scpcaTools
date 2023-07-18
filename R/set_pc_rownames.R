#' Label rownames of provided PCs and remove NA batch cells from PCs
#'
#' @param integrated_pcs The PCs with unlabeled rownames
#' @param batches Vector of cell-wise batch information whose length is the number of
#'   rows in `integrated_pcs`
#'
#' @return PCs with labeled rownames and NA batch cells removed
set_pc_rownames <- function(integrated_pcs, batches) {

  retain_indices <- which(!is.na(batches))
  batches <- batches[retain_indices]
  integrated_pcs <- integrated_pcs[retain_indices,]

  # Check dimensions
  if (nrow(integrated_pcs) != length(batches)) {
    stop("Incompatable PC and batch information dimensions after removing NAs.")
  }

  # Set PC rownames to be the batches
  rownames(integrated_pcs) <- batches

  return(integrated_pcs)
}
