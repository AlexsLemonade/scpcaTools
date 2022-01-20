#' Add cellhash sample information to altExperiment rowData
#'
#' @param sce SingleCellExperiment object.
#' @param barcode_table A data_frame of barcode_id and sample_id.
#' @param altexp_id The name of the alternative experiment that contains the cellhash data.
#'   Default value: "cellhash"
#'
#'
#' @return SingleCellExperiment with colData slot containing the demultiplexing calls
#'
#' @import SingleCellExperiment
#'
#'
#' @export
#'
#' @examples
#' \dontrun{

#' }
add_barcode_table <- function(sce, barcode_table, altexp_id = "cellhash"){
  # check that input is a SingleCellExperiment
  if(!is(sce, "SingleCellExperiment")){
    stop("sce must be a SingleCellExperiment object")
  }
  if(!altexp_id %in% altExpNames(sce)){
    stop(glue::glue("altexp_id `{altexp_id}` not found as an experiment in sce"))
  }
  # check that input is a SingleCellExperiment
  if(!(is(barcode_table, "data.frame") &&
     all(c("barcode_id", "sample_id") %in% barcode_table))){
    stop("barcode_table must be a data frame with columns `barcode_id` and `sample_id`")
  }

  # get barcodes in table
  altexp_rows <- rowData(altExp(sce, altexp_id)) |>
    as.data.frame() |>
    tibble::rowid_to_column("barcode_id")

  barcode_rowdata <- altexp_barcodes |>
    dplyr::left_join(barcode_table) |>
    dplyr::select("barcode_id", "sample_id")

  rowData(altExp(sce))$barcode_id <- barcode_rowdata$barcode_id
  rowData(altExp(sce))$sample_id <- barcode_rowdata$sample_id

  return(sce)
}


#' Add sample of origin calls using DropletUtils::hashedDrops()
#'
#' @param sce SingleCellExperiment object.
#' @param altexp_id The name of the alternative experiment that contains the cellhash data
#' @param ... Any additional arguments to be passed to DropletUtils::hashedDrops()
#'
#' @return SingleCellExperiment with colData slot containing the demultiplexing calls
#'
#' @import SingleCellExperiment
#'
#'
#' @export
#'
#' @examples
#' \dontrun{

#' }
add_hashedDrops_calls <- function(sce,
                                  altexp_id = "cellhash",
                                  ...){

  # check that input is a SingleCellExperiment
  if(!is(sce, "SingleCellExperiment")){
    stop("sce must be a SingleCellExperiment object")
  }
}
