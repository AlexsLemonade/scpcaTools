#' Merge a list of SCEs as preparation for formal integration
#'
#' This function takes an optionally-named (if named, ideally by a form of
#'  library ID) list of SingleCellExperiment (SCE) objects and merges them into
#'  one SCE object. At least some genes must be present in all SCEs in order to
#'  merge them. Currently any present altExps are not retained.
#'
#'  Original SCE contents are modified or retained as follows:
#'  - The resulting colData slot will include a new column specified by
#'    `batch_column` (default "library_id") that either holds the originating SCE
#'    object's name (referred to as `sce_name` here), of if it is unnamed then
#'    its index in the provided `sce_list`.
#'  - The resulting colData slot will include another new column `cell_id_column`
#'    (default "cell_id") that is will contain the SCE's column names (i.e. colData
#'    rownames). Often, but not always, the rowname holds a unique cell barcode.
#'    This column serves to match cells back to both originating batch and cell name.
#'  - The resulting colData rownames will be updated to match `cell_id`.
#'  - The resulting rowData slot column names will be appended with the given
#'    SCE's name, as `{sce_name}-{column_name}` except for columns whose names
#'    are indicated to preserve with the `preserve_rowdata_cols` argument.
#'
#'
#' @param sce_list A list of SingleCellExperiment objects. The list may optionally
#'  be named with batch information. If no names are provided, names will be generated
#'  based on the SCE's index.
#' @param batch_column A character value giving the resulting colData column name
#'  to differentiate originating SingleCellExperiment objects. Often these values
#'  are unique library IDs. Default value is `library_id`.
#' @param retain_coldata_cols A vector of colData columns which should be retained
#'  in the the final merged SCE object.
#' @param preserve_rowdata_cols A vector of column names that appear in originating
#'  SCE objects' rowData slots which should not be renamed with
#'  the given SCE object name or index is name is not given. These are generally
#'  columns which are not specific to the given library's preparation or statistics.
#'  For example, such a vector might contain items like "Gene", "ensembl_ids", etc.
#' @param cell_id_column A character value giving the resulting colData colum name
#'  to hold unique cell IDs formatted as their batch and original rowname. Default
#'  value is `cell_id`.
#'
#' @return A SingleCellExperiment object containing all SingleCellExperiment objects
#'   present in the inputted list
#' @export
#'
#' @import SingleCellExperiment
merge_sce_list <- function(sce_list = list(),
                           batch_column = "library_id",
                           retain_coldata_cols = c("sum",
                                                   "barcode",
                                                   "detected",
                                                   "total",
                                                   "subsets_mito_sum",
                                                   "subsets_mito_detected",
                                                   "subsets_mito_percent",
                                                   "miQC_pass",
                                                   "prob_compromised"),
                           preserve_rowdata_cols = NULL,
                           cell_id_column = "cell_id") {

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
    warning("All pre-existing colData will be removed from the the merged SCE.
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
  sce_list <- sce_list |>
    purrr::imap(prepare_sce_for_merge,
                batch_column = batch_column,
                cell_id_column = cell_id_column,
                shared_features = shared_features,
                retain_coldata_cols = retain_coldata_cols,
                preserve_rowdata_cols = preserve_rowdata_cols)


  # Create the merged SCE from the processed list ------------------
  merged_sce <- do.call(cbind, sce_list)


  # Return merged SCE object ----------------------------
  return(merged_sce)

}



#' Helper function to prepare an SCE object for merging
#'
#' @param sce The SCE object to be prepared
#' @param sce_name The name of the SCE object
#' @param batch_column The name of the batch column which will be added to the
#'   colData slot
#' @param cell_id_column The name of the cell_id column which will be added to the
#'   colData slot
#' @param shared_features A vector of features (genes) that all SCEs to be merged
#'   have in common
#' @param retain_coldata_cols A vector of columns to retain in the colData slot.
#'   If columns are missing from the data, they will be filled with `NA` values.
#'   The exceptions to this are `barcode_column` and `batch_column` which will be
#'   populated with respective information.
#' @param preserve_rowdata_cols A vector of rowData columns which should not be
#'   renamed
#'
#' @return An updated SCE that is prepared for merging
prepare_sce_for_merge <- function(sce,
                                  sce_name,
                                  batch_column,
                                  cell_id_column,
                                  shared_features,
                                  retain_coldata_cols,
                                  preserve_rowdata_cols) {

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
    glue::glue("{sce_name}-{original_colnames}")
  )

  ##### colData #####
  observed_coldata_names <- names(colData(sce))

  # Ensure all columns are present in all SCEs by adding `NA` columns as needed
  missing_columns <- setdiff(retain_coldata_cols, observed_coldata_names)
  for (missing_col in missing_columns) {
    # Create the missing column only if it should be retained, and it is NOT the
    # barcode column
    if (missing_col %in% retain_coldata_cols) {
      colData(sce)[[missing_col]] <- NA
    }
  }

  # Retain only the columns present in `retain_coldata_cols`
  # Use drop=FALSE to ensure result is a DataFrame
  colData(sce) <- colData(sce)[, retain_coldata_cols, drop = FALSE]

  # Add batch column
  sce[[batch_column]] <- sce_name

  # Add cell_id column
  sce[[cell_id_column]] <- colnames(sce)

  # Add `sce_name` to colnames so cell ids can be mapped to originating SCE
  colnames(sce) <- glue::glue("{sce_name}-{colnames(sce)}")

  # return the processed SCE
  return(sce)
}
