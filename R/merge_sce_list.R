#' Merge a list of SCEs as preparation for formal integration
#'
#' This function takes an optionally-named (if named, ideally by a form of
#'  library ID) list of SingleCellExperiment (SCE) objects and merges them into
#'  one SCE object. The resulting colData slot will contain information in the
#'  given batch column to differentiate originating SCE objects based on the SCE
#'  names, or index if no names are provided. If not already present, cell barcodes
#'  will also be created as a colData column. Only features (genes) that are
#'  present in all provided SCE objects will be retained. To preserve original
#'  feature-level information, all final rowData slot column names will be
#'  appended with the given SCE's name, as `<column_name>-<sce_name>` sxcept for
#'  columns indicated to preserve. Currently any present altExps are not retained.
#'
#' @param sce_list A list of SingleCellExperiment objects. The list may optionally
#'  be named with batch information. If no names are provided, names will be generated
#'  based on the SCE's index.
#' @param batch_column A character value giving the resulting colData column name
#'  to differentiate originating SingleCellExperiment objects. Often these values
#'  are library IDs. Default value is `library_id`.
#' @param barcode_column A character value giving the colData column name that
#'  stores cell barcodes. If this column does not yet exist it will be created.
#'  Default value is  `barcode`.
#' @param retain_coldata_cols A vector of colData columns which should be retained
#'  in the the final merged SCE object.
#' @param preserve_rowdata_cols A vector of column names that appear in originating
#'  SCE objects' rowData slots which should not be renamed with
#'  the given SCE object name or index is name is not given. These are generally
#'  columns which are not specific to the given library's preparation or statistics.
#'  For example, such a vector might contain items like "Gene", "ensembl_ids", etc.
#'
#' @return A SingleCellExperiment object containing all SingleCellExperiment objects
#'   present in the inputted list
#' @export
#'
#' @import SingleCellExperiment
merge_sce_list <- function(sce_list = list(),
                           batch_column = "library_id",
                           barcode_column = "barcode",
                           retain_coldata_cols = c("sum",
                                                   "detected",
                                                   "total",
                                                   "subsets_mito_sum",
                                                   "subsets_mito_detected",
                                                   "subsets_mito_percent",
                                                   "miQC_pass",
                                                   "prob_compromised"),
                           preserve_rowdata_cols = NULL) {

  # Check `sce_list`----------------------
  if (is.null(names(sce_list))) {
    warning(
      glue::glue(
        "Individual SCE objects in `sce_list` are not named, so batches will be
        named based on their list index in the merged SCE object.")
    )
    names(sce_list) <- 1:length(sce_list)
  }

  if (length(sce_list) < 2) {
    warning("There are fewer than two SCE objects in the provided `sce_list` so there is nothing to merge.")
    # Early return:
    return(sce_list)
  }

  # Check `retain_coldata_cols` ----------------
  if (length(retain_coldata_cols) == 0) {
    warning("All colData will be removed from the the merged SCE.
     Please check that `retain_coldata_cols` was correctly specified.")
  }

  # Subset SCEs to shared features and ensure appropriate naming ------------------

  # First, obtain intersection among all SCE objects
  shared_features <- sce_list |>
    purrr::map(rownames) |>
    purrr::reduce(intersect)

  if (length(shared_features) == 0) {
    stop("There are no shared features among provided SCE objects.
         They cannot be merged.")
  }

  # Second, determine all the column names that are present in any SCE so it can
  #  be created in any missing SCEs with `NA` values
  all_colnames <- purrr::map(sce_list, ~names(colData(.))) |>
    unlist() |>
    unname() |>
    unique()

  # Check that the `retain_coldata_cols` are present in at least one SCE, and
  #  error if the column exists nowhere.
  if (!(any(retain_coldata_cols %in% all_colnames))) {
    stop("Error: The provided `retain_coldata_cols` are not present in any SCEs.")
  }

  # Prepare SCEs
  sce_list <- purrr::imap(sce_list,
                          prepare_sce_for_merge,
                          batch_column,
                          barcode_column,
                          shared_features,
                          # ensure that barcodes are retained
                          c(barcode_column, retain_coldata_cols),
                          preserve_rowdata_cols,
                          all_colnames)


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
#' @param barcode_column The name of the column to be created (or not if it
#'   already exists) in the colData slot to store cell barcodes
#' @param shared_features A vector of features (genes) that all SCEs to be merged
#'   have in common
#' @param retain_coldata_cols A vector of columns to retain in the colData slot
#' @param preserve_rowdata_cols A vector of rowData columns which should not be
#'   renamed
#' @param expected_coldata_names A vector of column names that are expected to be
#'   present in the SCE. Any missing columns will be created with `NA` values.
#'
#' @return An updated SCE that is prepared for merging
prepare_sce_for_merge <- function(sce,
                                  sce_name,
                                  batch_column,
                                  barcode_column,
                                  shared_features,
                                  retain_coldata_cols,
                                  preserve_rowdata_cols,
                                  expected_coldata_names) {

  # Current functionality does not retain any present altExps
  sce <- removeAltExps(sce)

  # Subset to shared features
  sce <- sce[shared_features,]

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
  observed_coldata_names <- names(colData(sce))

  # If the barcode column is not already present, create it
  if (!(barcode_column %in% observed_coldata_names)) {
    colData(sce)[,barcode_column] <- rownames(colData(sce))
    observed_coldata_names <- c(observed_coldata_names, barcode_column)
  }

  # Ensure barcode_column is expected and retained
  if (!(barcode_column %in% expected_coldata_names)) {
    expected_coldata_names <- c(expected_coldata_names, barcode_column)
  }
  if (!(barcode_column %in% retain_coldata_cols)) {
    retain_coldata_cols <- c(retain_coldata_cols, barcode_column)
  }

  # Ensure all columns are present in all SCEs by adding `NA` columns as needed
  missing_columns <- setdiff(expected_coldata_names, observed_coldata_names)
  for (missing_col in missing_columns) {
    # Create the missing column only if it should be retained
    if (missing_col %in% retain_coldata_cols) {
      colData(sce)[[missing_col]] <- NA
    }
  }

  # Retain only the columns present in `retain_coldata_cols`
  colData(sce) <- colData(sce)[, retain_coldata_cols, drop=FALSE]

  # Add `sce_name` to colnames so cell ids can be mapped to originating SCE
  colnames(sce) <- glue::glue("{colnames(sce)}-{sce_name}")

  # Add batch column
  sce[[batch_column]] <- sce_name


  # return the processed SCE
  return(sce)
}
