#' Merge a list of SCEs as preparation for formal integration
#'
#' This function takes a named (ideally by a form of library ID) list of
#'  SingleCellExperiment (SCE) objects and merges them into one SCE object.
#'  All names in the colData slot must match across provided SCE objects.
#'  The resulting colData slot will contain information in the given batch column
#'  to differentiate originating SCE objects. Only features (genes) that are present
#'  in all provided SCE objects will be retained. To preserve original feature-level
#'  information, all final rowData slot column names will be appended with the
#'  given SCE's name, as `<column_name>-<sce_name>`. Currently any present altExps
#'  are not retained.
#'
#' @param sce_list A named list of SingleCellExperiment objects
#' @param batch_column A character value giving the resulting colData column name
#'  to differentiate originating SingleCellExperiment objects Default value is `batch`.
#' @param retain_coldata_cols A vector of colData columns which should be retained
#'  in the main experiment of the final merged SCE object.
#' @param preserve_rowdata_cols A vector of column names that appear in originating
#'  SingleCellExperiment objects' rowData slots which should not be renamed with
#'  the given SingleCellExperiment object name. These are generally columns which are
#'  not specific to the given library's preparation or statistics. For example,
#'  such an array might contain items like "Gene", "ensembl_ids", "gene_symbol", etc.
#'
#' @return A SingleCellExperiment object containing all SingleCellExperiment objects
#'   present in the inputted list
#' @export
#'
#' @import SingleCellExperiment
merge_sce_list <- function(sce_list = list(),
                           batch_column = "batch",
                           retain_coldata_cols = c("sum",
                                                   "detected",
                                                   "total",
                                                   "subsets_mito_sum",
                                                   "subsets_mito_detected",
                                                   "subsets_mito_percent",
                                                   "miQC_pass"),
                           preserve_rowdata_cols = NULL) {

  # Check `sce_list`----------------------
  if (is.null(names(sce_list))) {
    stop("Individual SingleCellExperiment objects in `sce_list` must be named.")
  }
  if (length(sce_list) < 2) {
    warning("There are fewer than two SCE objects in the provided `sce_list` so there is nothing to merge.")
    # Early return:
    return(sce_list)
  }

  # Check input vectors ----------------
  if(length(retain_coldata_cols) == 0) {
    warning("All colData will be removed from the the merged SCE.
            Please check that `retain_coldata_cols` was correctly specified.")
  }

  # Subset SCEs to shared features and ensure appropriate naming ------------------

  # First, obtain intersection among all SCE objects
  shared_features <- sce_list |>
    purrr::map(rownames) |>
    purrr::reduce(intersect)

  # Prepare SCEs
  sce_list <- purrr::imap(sce_list,
                          prepare_sce_for_merge,
                          batch_column,
                          shared_features,
                          retain_coldata_cols,
                          preserve_rowdata_cols)


  # Create the merged SCE from the processed list ------------------
  merged_sce <- do.call(cbind, sce_list)


  # Return merged SCE object ----------------------------
  return(merged_sce)

}



#' Helper function to prepare an SCE object for merging
#'
#' @param sce The SCE object to be prepared
#' @param sce_name The name of the SCE object
#' @param batch_column The name of the batch column will which be added to the
#'   colData slot
#' @param shared_features A vector of features (genes) that all SCEs to be merged
#'   have in common
#' @param retain_coldata_cols A vector of columns to retain in the colData slot
#' @param preserve_rowdata_cols A vector of rowData columns which should not be
#'   renamed
#'
#' @return An updated SCE that is prepared for merging
prepare_sce_for_merge <- function(sce,
                                  sce_name,
                                  batch_column,
                                  shared_features,
                                  retain_coldata_cols,
                                  preserve_rowdata_cols) {

  # Current functionality does not retain any present altExps
  sce <- removeAltExps(sce)

  # Subset to shared features
  sce <- sce[shared_features,]

  # Add `sce_name` to colData row names so cell barcodes can be mapped to originating SCE
  rownames(colData(sce)) <- glue::glue("{rownames(colData(sce))}-{sce_name}")

  ##### rowData #####
  # Add `sce_name` ID to rowData column names except for those
  #  present in `preserve_rowdata_cols`
  original_colnames <- colnames(rowData(sce))

  colnames(rowData(sce)) <- ifelse(
    original_colnames %in% preserve_rowdata_cols,
    original_colnames,
    glue::glue("{original_colnames}-{sce_name}")
  )

  ##### colData #####
  # Retain only the columns present in `retain_coldata_cols`
  coldata_names <- names(colData(sce))
  if (!(all(retain_coldata_cols %in% coldata_names))) {
    stop("Error: Provided columns to retain are not present
            in all SCE objects.")
  }
  colData(sce) <- colData(sce) |>
    as.data.frame() |>
    dplyr::select( dplyr::all_of(retain_coldata_cols) ) |>
    S4Vectors::DataFrame()

  # Add batch column
  sce[[batch_column]] <- sce_name


  # return the processed SCE
  return(sce)
}
