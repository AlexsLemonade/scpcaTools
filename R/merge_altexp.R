#' Add an alternative experiment to a SingleCellExperiment by "left join"
#'
#' This function takes SingleCellExperiment (SCE) and adds an Alternative Experiment (AltExp) to it.
#' If the two experiments do not have the same set of cells, this will use the first
#' SCE as the base set, adding columns (cells) to the AltExp so that the matrices match on that dimension.
#' Columns (cells) that are not found in the base experiment are discarded.
#'
#' @param sce a SingleCellExperiment object
#' @param alt_exp a second SummarizedExperiment to add as an alternative experiment
#'   Most likely, this will also be a SingleCellExperiment, but it does not have to be.
#'   The AltExp in the output will be a SingleCellExperiment
#' @param alt_name The name to use for the alternative experiment slot
#'
#' @return A SingleCellExperiment with an AltExp slot containing another SingleCellExperiment.
#' @export
#'
#' @import SingleCellExperiment
#' @importFrom S4Vectors metadata
#'
merge_altexp <- function(sce, alt_exp, alt_name){
  if(!is(sce, "SingleCellExperiment")){
    stop("sce must be a SingleCellExperiment object")
  }
  if(!is(alt_exp, "SummarizedExperiment")){
    stop("alt_exp must be a SummarizedExperiment (or SingleCellExperiment) object")
  }
  sce_cells <- colnames(sce)
  alt_cells <- colnames(alt_exp)
  alt_rows <- rownames(alt_exp)
  # perform "left join": create blank columns for cells without alt counts
  alt_missing_cells <- sce_cells[!sce_cells %in% alt_cells]
  empty_counts <- Matrix::Matrix(0,
                                 nrow = length(alt_rows),
                                 ncol = length(alt_missing_cells),
                                 dimnames = list(alt_rows, alt_missing_cells))
  alt_counts_all <- cbind(counts(alt_exp), empty_counts)
  # sort to original order (dropping alt-only cells)
  alt_counts_all <- alt_counts_all[, sce_cells]
  altExp(sce, alt_name) <- SingleCellExperiment(assays = list(counts = alt_counts_all),
                                                metadata = metadata(alt_exp))
  return(sce)
}
