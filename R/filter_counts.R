#' Filter counts matrix using DropletUtils::emptyDropsCellRanger
#'
#' This function will filter a SingleCellExperiment object using DropletUtils::emptyDropsCellRanger() by default,
#'   or DropletUtils::emptyDrops(), as well as any associated alternative experiments. If mean expression and percent detected
#'   were previously calculated in the columns `mean` and `detected`, respectively, these
#'   will be removed from both the main and alternative experiments.
#'
#' @param sce SingleCellExperiment with unfiltered gene x cell counts matrix.
#' @param cr_like Logical indicating whether or not to use DropletUtils::emptyDropsCellRanger.
#'   Default is set to TRUE.
#' @param fdr_cutoff FDR cutoff to use for DropletUtils::emptyDropsCellRanger or DropletUtils::emptyDrops.
#'   Default is 0.01.
#' @param seed An optional random seed for reproducibility.
#' @param umi_cutoff The minimum UMI count for cells to pass filtering, only used if emptyDropsCellRanger and emptyDrops fails.
#'   Default is 100.
#' @param ... Any arguments to be passed into DropletUtils::emptyDropsCellRanger or DropletUtils::emptyDrops.
#'
#' @return SingleCellExperiment with filtered gene x cell matrix.
#'
#' @import SingleCellExperiment
#'
#' @export
#'
#' @examples
#' \dontrun{
#' filter_counts(sce = sce_object)
#' }
filter_counts <- function(sce, cr_like = TRUE, fdr_cutoff = 0.01, seed = NULL, umi_cutoff= 100, ...) {

  set.seed(seed)

  if(!is(sce,"SingleCellExperiment")){
    stop("Input must be a SingleCellExperiment object.")
  }

  if(!is.logical(cr_like)){
    stop("cr_like must be set as TRUE or FALSE")
  }

  if(!is.numeric(umi_cutoff) & umi_cutoff >= 0){
    stop("umi_cutoff must be an integer greater than or equal to 0")
  }

  # calculate probability of being an empty droplet
  if(cr_like){
    empty_df <- tryCatch(
      DropletUtils::emptyDropsCellRanger(m = round(counts(sce))),
      error = function(x){NULL}
    )
    metadata(sce)$filtering_method <- "emptyDropsCellRanger"
  } else {
    empty_df <- tryCatch(
      DropletUtils::emptyDrops(counts(sce), ...),
      error = function(x){NULL}
    )
    metadata(sce)$filtering_method <- "emptyDrops"
  }

  # if emptyDrops doesn't fail, filter by FDR cutoff
  if(!is.null(empty_df)){
    if(any(empty_df$FDR > fdr_cutoff & empty_df$Limited, na.rm = TRUE)){
      warning(glue::glue("`niters` may be set too low for emptyDrops filtering.",
                         " Current value is {empty_df@metadata$niters}."))
    }
    cells <- rownames(empty_df)[which(empty_df$FDR <= fdr_cutoff)]
  } else {
    # if emptyDrops fails, filter using a hard UMI threshold
    cells <- sce$sum > umi_cutoff
    # replace filtering method with UMI cutoff as method
    metadata(sce)$filtering_method <- "UMI cutoff"
    metadata(sce)$umi_cutoff <- umi_cutoff
  }

  # subset original counts matrix by cells that pass filter
  sce <- sce[, cells]

  # remove feature stats that are no longer valid
  drop_cols <- colnames(rowData(sce)) %in% c('mean', 'detected')
  rowData(sce) <- rowData(sce)[!drop_cols]
  for (alt in altExpNames(sce)) { # alt expts too
    drop_cols = colnames(rowData(altExp(sce, alt))) %in% c('mean', 'detected')
    rowData(altExp(sce, alt)) <- rowData(altExp(sce, alt))[!drop_cols]
  }
  return(sce)
}
