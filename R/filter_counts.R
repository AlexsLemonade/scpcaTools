#' Filter counts matrix using DropletUtils::emptyDrops
#'
#' @param sce SingleCellExperiment with unfiltered gene x cell counts matrix.
#' @param fdr_cutoff FDR cutoff to use for DropletUtils::emptyDrops.
#'   Default is 0.01.
#' @param ... Any arguments to be passed into DropletUtils::emptyDrops.
#'
#' @return SingleCellExperiment with filtered gene x cell matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' filter_counts(sce = sce_object)
#' }
filter_counts <- function(sce, fdr_cutoff = 0.01, ...) {

  if(!is(sce,"SingleCellExperiment")){
    stop("Input must be a SingleCellExperiment object.")
  }

  # calculate probability of being an empty droplet
  empty_df <- DropletUtils::emptyDrops(sce, ...)
  if(any(empty_df$FDR > fdr_cutoff & empty_df$Limited)){
    warning(glue::glue("`niters` may be set too low for emptyDrops filtering.",
                    " Current value is {empty_df@metadata$niters}."))
  }
  cells <- rownames(empty_df)[which(empty_df$FDR <= fdr_cutoff)]

  # subset original counts matrix by cells that pass filter
  sce <- sce[, cells]
  return(sce)
}
