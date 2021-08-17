#' Add an Alternative Experiment to a SingleCellExperiment by "left join"
#'
#' @param sce
#' @param feature_counts
#' @param feature_name
#' @param feature_metadata
#'
#' @return
#' @export
#'
#' @examples
add_altexp <- function(sce, feature_counts, feature_name, feature_metadata = NULL){
  cell_barcodes <- colnames(sce)
  feature_cell_barcodes <- colnames(feature_counts)
  features <- rownames(feature_counts)
  # perform "left join": create blank columns for cells without feature counts
  feature_missing_cells <- cell_barcodes[!cell_barcodes %in% feature_cell_barcodes]
  empty_counts <- Matrix::Matrix(0,
                                 nrow = length(features),
                                 ncol = length(feature_missing_cells),
                                 dimnames = list(features, feature_missing_cells))
  feature_counts_all <- cbind(feature_counts, empty_counts)
  feature_counts_all <- feature_counts_all[, cell_barcodes]
  altExp(sce, feature_name) <- SingleCellExperiment(assays = list(counts = feature_counts_all),
                                                    metadata = feature_metadata)
  return(sce)
}
