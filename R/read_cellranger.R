#' Read in counts data processed with Cellranger
#'
#' @param quant_dir Path to directory where output files are located.
#'
#' @return SingleCellExperiment of gene x cell counts matrix
#' @export
#'
read_cellranger <- function(quant_dir) {

  cellranger_file <- file.path(quant_dir, "outs", "filtered_feature_bc_matrix.h5")
  if(!file.exists(cellranger_file)){
    stop("Missing filtered_feature_bc_matrix.h5 file from cellranger output")
  }

  sce <- DropletUtils::read10xCounts(cellranger_file,
                       sample.names = basename(quant_dir),
                       col.names = TRUE)

  # for consistency with other quantifiers:
  # change the column names just the barcode value, which is the first part of the barcode name
  # drop colData
  colnames(sce) <- str_extract(colnames(sce), "^([ACGT]+)")
  SummarizedExperiment::colData(sce) <- NULL
  return(sce)
}
