#' Filter counts matrix using DropletUtils::emptyDrops
#'
#' @param sce SingleCellExperiment with unfiltered gene x cell counts matrix
#'
#' @return SingleCellExperiment with filtered gene x cell matrix
#' @export
#'
#' @examples
#' \dontrun{
#' filter_counts(sce = sce_object)
#' }
filter_counts <- function(sce) {

  if(!is(sce,"SingleCellExperiment")){
    stop("Input must be a SingleCellExperiment object.")
  }

  # grab counts from single cell experiment
  counts <- counts(sce)

  # calculate probability of being an empty droplet
  empty_df <- emptyDrops(counts)
  cells <- rownames(empty_df[which(empty_df$Limited == "TRUE" & empty_df$FDR <= 0.01),])

  # subset original counts matrix by cells that pass filter
  counts <- counts[, cells]
  sce <- SingleCellExperiment(list(counts = counts))
  return(sce)
}
