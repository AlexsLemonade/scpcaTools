#' Add cellhash sample information to altExperiment rowData
#'
#' @param sce SingleCellExperiment object.
#' @param hashsample_table A data_frame of barcode_id and sample_id.
#' @param altexp_id The name of the alternative experiment that contains the cellhash data.
#'   Default value: "cellhash"
#'
#'
#' @return SingleCellExperiment with rowData for the altExp containing barcode and sample ids
#'
#' @import SingleCellExperiment
#' @importFrom rlang .data
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # add sample information for a cellhash table
#' add_hashsample_table(sce = sce,
#'                      hashsample_table = barcode_ids)
#' }
add_hashsample_table <- function(sce, hashsample_table, altexp_id = "cellhash"){
  # check that input is a SingleCellExperiment
  if(!is(sce, "SingleCellExperiment")){
    stop("sce must be a SingleCellExperiment object")
  }
  if(!altexp_id %in% altExpNames(sce)){
    stop(glue::glue("altexp_id `{altexp_id}` not found as an experiment in sce"))
  }
  # check that input is a SingleCellExperiment
  if(!(is(hashsample_table, "data.frame") &&
     all(c("barcode_id", "sample_id") %in% colnames(hashsample_table)))){
    stop("hashsample_table must be a data frame with columns `barcode_id` and `sample_id`")
  }

  # get barcodes in table
  altexp_rowdata <- rowData(altExp(sce, altexp_id)) |>
    as.data.frame() |>
    tibble::rownames_to_column("barcode_id")

  barcode_rowdata <- altexp_rowdata |>
    dplyr::left_join(hashsample_table, by = "barcode_id") |>
    dplyr::select("barcode_id", "sample_id")

  # test that each barcode has a corresponding sample
  missing_barcodes <- barcode_rowdata |>
    dplyr::filter(is.na(.data$sample_id)) |>
    dplyr::pull(.data$barcode_id)
  if (length(missing_barcodes)){
    warning(paste0("The following `barcode_id`s do not have corresponding `sample_id`s:",
                   paste0(missing_barcodes, collapse = ", ")))
  }

  rowData(altExp(sce))$barcode_id <- barcode_rowdata$barcode_id
  rowData(altExp(sce))$sample_id <- barcode_rowdata$sample_id

  return(sce)
}

