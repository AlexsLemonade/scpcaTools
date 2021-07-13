#' Convert colData from SingleCellExperiment to a data.frame
#'
#' @param sce SingleCellExperiment object.
#'
#' @return data.frame containing colData from the SingleCellExperiment object
#'   and a column, "cell_id", corresponding to the cell_barcodes or columns of the
#'   SingleCellExperiment object.
#'
#' @import SingleCellExperiment
#'
#' @export
#'
coldata_to_df <- function(sce) {

  # make sure that input is a SingleCellExperiment
  if(!is(sce, "SingleCellExperiment")){
    stop("Input must be a SingleCellExperiment object.")
  }

  # make sure that input has colData
  if(ncol(colData(sce)) == 0){
    warning("SingleCellExperiment has empty colData slot.")
  }

  # convert colData to data.frame with column containing cell barcodes
  df <- as.data.frame(colData(sce)) %>%
    tibble::rownames_to_column(var = "cell_id")
  return(df)
}

#' Convert rowData from SingleCellExperiment to a data.frame
#'
#' @param sce SingleCellExperiment object.
#'
#' @return data.frame containing rowData from the SingleCellExperiment object
#'   and a column, "gene_id", corresponding to the gene names or rows of the
#'   SingleCellExperiment object.
#' @export
#'
rowdata_to_df <- function(sce) {

  # make sure that input is a SingleCellExperiment
  if(!is(sce, "SingleCellExperiment")){
    stop("Input must be a SingleCellExperiment object.")
  }

  # make sure that input has rowData
  if(any(is.na(rowData(sce)))){
    stop("SingleCellExperiment has empty rowData slot")
  }

  # convert rowData to data.frame column containing gene names
  df <- as.data.frame(rowData(sce)) %>%
    tibble::rownames_to_column(var = "gene_id")
  return(df)
}
