#' Add cellhash sample information to altExperiment rowData
#'
#' @param sce SingleCellExperiment object.
#' @param hashsample_table A data_frame of barcode_id and sample_id.
#' @param altexp_id The name of the alternative experiment that contains the cellhash data.
#'   Default value: "cellhash"
#' @param removed_unlabeled Should data for barcodes with no corresponding sample be removed
#'   from the altExp. Default `FALSE`.
#' @param replace_rownames Replace rownames for the altExp with sample_ids.
#'   Requires `remove_unlabeled` to be TRUE. Default `FALSE`.
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
#' add_cellhash_ids(sce = sce,
#'                  hashsample_table = barcode_ids)
#' }
add_cellhash_ids <- function(sce, hashsample_table,
                             altexp_id = "cellhash",
                             remove_unlabeled = FALSE,
                             replace_rownames = FALSE){
  # check that input is a SingleCellExperiment
  if(!is(sce, "SingleCellExperiment")){
    stop("sce must be a SingleCellExperiment object")
  }
  if(!altexp_id %in% altExpNames(sce)){
    stop(glue::glue("altexp_id `{altexp_id}` not found as an experiment in sce"))
  }
  # check that hashsample_table is a data frame containing a sample_id and barcode_id column
  if(!(is(hashsample_table, "data.frame") &
     all(c("barcode_id", "sample_id") %in% colnames(hashsample_table)))){
    stop("hashsample_table must be a data frame with columns `barcode_id` and `sample_id`")
  }
  # check requirements for replace_rownames
  if(replace_rownames & !remove_unlabeled){
    stop("replace_rownames = TRUE requires remove_unlabeled to be TRUE")
  }

  # get rowData from sce and add barcode_id column if needed
  altexp_rowdata <- rowData(altExp(sce, altexp_id)) |>
    as.data.frame()
  if (!"barcode_id" %in% colnames(altexp_rowdata)){
    altexp_rowdata <- altexp_rowdata |>
      tibble::rownames_to_column("barcode_id")
  }


  barcode_rowdata <- altexp_rowdata |>
    dplyr::left_join(hashsample_table, by = "barcode_id") |>
    dplyr::select("barcode_id", "sample_id") |>
    dplyr::distinct()

  if (any(duplicated(barcode_rowdata$barcode_id))){
    stop("`barcode_ids` in `hashsample_table` can not be duplicated.")
  }

  # test that each barcode has a corresponding sample
  missing_barcodes <- barcode_rowdata |>
    dplyr::filter(is.na(.data$sample_id)) |>
    dplyr::pull(.data$barcode_id)

  if (remove_unlabeled){
    missing_rows <- rownames(altExp(sce, altexp_id)) %in% missing_barcodes
    altExp(sce, altexp_id) <- altExp(sce, altexp_id)[!missing_rows,]
    barcode_rowdata <- barcode_rowdata[!missing_rows,]
  } else if (length(missing_barcodes) > 0){
    warning(paste0("The following `barcode_id`s do not have corresponding `sample_id`s:",
                   paste0(missing_barcodes, collapse = ", ")))
  }

  rowData(altExp(sce, altexp_id))$barcode_id <- barcode_rowdata$barcode_id
  rowData(altExp(sce, altexp_id))$sample_id <- barcode_rowdata$sample_id

  if (replace_rownames){
    rownames(altExp(sce, altexp_id)) <- rowData(altExp(sce, altexp_id))$sample_id
  }

  return(sce)
}


#' Add cellhash demultiplexing results using DropletUtils::hashedDrops
#'
#' @param sce SingleCellExperiment object.
#' @param altexp_id The name of the alternative experiment that contains the cellhash data.
#'   Default value: "cellhash"
#' @param ... Other arguments to pass to DropletUtils::hashedDrops
#'
#'
#' @return SingleCellExperiment with colData columns containing the demultiplexing calls
#'
#' @import SingleCellExperiment
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # add cell calls from DropletUtils::hashedDrops()
#'   add_demux_hashedDrops(sce = sce)
#' }
add_demux_hashedDrops <- function(sce, altexp_id = "cellhash", ...){
  # check that input is a SingleCellExperiment
  if(!is(sce, "SingleCellExperiment")){
    stop("sce must be a SingleCellExperiment object")
  }
  if(!altexp_id %in% altExpNames(sce)){
    stop(glue::glue("altexp_id `{altexp_id}` not found as an experiment in sce"))
  }
  # get ids to use: if sample ids are present, use those,
  #  otherwise barcodes or rownames, in that order
  sample_ids <- rowData(altExp(sce, altexp_id))$sample_id
  if(is.null(sample_ids)){
    sample_ids <- rowData(altExp(sce, altexp_id))$barcode_id
  }
  if (is.null(sample_ids)){
    warning("No sample ids are present for demux results, using row names")
    sample_ids <- rownames(altExp(sce, altexp_id))
  }

  # calculate cellhash results
  hash_result <- DropletUtils::hashedDrops(altExp(sce, altexp_id), ...)
  stopifnot(all.equal(rownames(hash_result), colnames(sce)))

  # extract results with sample ids for main rowData table
  hashedDrops_bestid = sample_ids[hash_result$Best]
  hashedDrops_id = ifelse(hash_result$Confident,
                          hashedDrops_bestid,
                          NA_character_)


  ### add all table results to altExp rowData

  # rename results with prefix
  colnames(hash_result) <- paste0("hashedDrops_", colnames(hash_result))
  # add hashedDrops results to altExp
  colData(altExp(sce, altexp_id))[, colnames(hash_result)] <- hash_result

  # add labeled results
  altExp(sce, altexp_id)$hashedDrops_id <- hashedDrops_id
  altExp(sce, altexp_id)$hashedDrops_bestid <- hashedDrops_bestid
  # add id to main sce table
  sce$hashedDrops_id <- hashedDrops_id

  return(sce)
}
