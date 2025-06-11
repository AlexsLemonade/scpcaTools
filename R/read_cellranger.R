#' Read in counts data processed with Cell Ranger
#'
#' @param cellranger_h5_file Path to H5AD file output from Cell Ranger.
#'
#' @return SingleCellExperiment of gene x cell counts matrix
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Import data from cellranger output file
#' read_cellranger(cellranger_h5_file)
#' }
#'
read_cellranger <- function(cellranger_h5_file) {

  if (!file.exists(cellranger_h5_file)) {
    stop(glue::glue("{cellranger_h5_file} does not exist"))
  }

  sce <- DropletUtils::read10xCounts(
    cellranger_h5_file,
    col.names = TRUE
  )

  # for consistency with other quantifiers:
  # change the column names just the barcode value, which is the first part of the barcode name
  # drop colData and metadata
  colnames(sce) <- str_extract(colnames(sce), "^([ACGT]+)")
  SummarizedExperiment::colData(sce) <- NULL
  SummarizedExperiment::metadata(sce) <- list()
  return(sce)
}
