#' Title
#'
#' @param quant_dir
#'
#' @return
#' @export
#'
#' @examples
read_cellranger <- function(quant_dir) {

  sce <- DropletUtils::read10xCounts(file.path(quant_dir,
                                              "outs", "filtered_feature_bc_matrix.h5"),
                                    sample.names = basename(quant_dir),
                                    col.names = TRUE)

  # for consistency with other quantifiers:
  # change the column names just the barcode value, which is the first part of the barcode name
  # drop colData
  colnames(sce) <- stringr::str_extract(colnames(sce), "^([ACGT]+)")
  SummarizedExperiment::colData(sce) <- NULL
  return(sce)
}
