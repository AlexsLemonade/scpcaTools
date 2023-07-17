#' Remove NA batch cells from PCs for use in integration metric calculations
#'
#' @param integrated_pcs The PCs with NA batch cells to be removed
#' @param batches Vector of cell-wise batch information whose length is the number of
#'   rows in `integrated_pcs`
#'
#' @return List with two items: `pcs`: PCs with NA batch cells removed and
#'   rownames assigned as batch; `indices`: Which rows were _retained_
remove_batch_nas_from_pcs <- function(integrated_pcs, batches) {

  retain_indices <- which(!is.na(batches))
  batches <- batches[retain_indices]
  integrated_pcs <- integrated_pcs[retain_indices,]

  # Check dimensions
  if (nrow(integrated_pcs) != length(batches)) {
    stop("Incompatable PC and batch information dimensions after removing NAs.")
  }

  # Set PC rownames to be the batches
  rownames(integrated_pcs) <- batches

  return(
    list(
      pcs = integrated_pcs,
      indices = retain_indices
    )
  )
}
