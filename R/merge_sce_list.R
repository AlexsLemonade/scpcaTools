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
#'    object's name (referred to as `sce_name` here), or if it is unnamed then
#'    its index in the provided `sce_list`.
#'  - The resulting colData slot will include another new column `cell_id_column`
#'    (default "cell_id") that will contain the SCE's original column names (i.e.
#'    original colData rownames). Often, but not always, this rowname holds a
#'    unique cell barcode.
#'  - The resulting colData rownames will be be prefixed with `{sce_name-}`.
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
#' @param cell_id_column A character value giving the resulting colData column name
#'  to hold unique cell IDs formatted as their original row name. Default
#'  value is `cell_id`.
#'
#' @return A SingleCellExperiment object containing all SingleCellExperiment objects
#'   present in the inputted list
#' @export
#'
#' @import SingleCellExperiment
merge_sce_list <- function(sce_list = list(),
                           batch_column = "library_id",
                           retain_coldata_cols = c(
                             "sum",
                             "detected",
                             "total",
                             "subsets_mito_sum",
                             "subsets_mito_detected",
                             "subsets_mito_percent",
                             "miQC_pass",
                             "prob_compromised",
                             "barcode"
                           ),
                           preserve_rowdata_cols = NULL,
                           cell_id_column = "cell_id") {
  # Check `sce_list`----------------------
  if (is.null(names(sce_list))) {
    warning(
      glue::glue(
        "Individual SCE objects in `sce_list` are not named, so batches will be
        named based on their list index in the merged SCE object."
      )
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
  all_colnames <- purrr::map(sce_list, ~ names(colData(.))) |>
    unlist() |>
    unname() |>
    unique()

  # Check that the `retain_coldata_cols` are present in at least one SCE, and
  #  error if the column exists nowhere.
  if (!(any(retain_coldata_cols %in% all_colnames))) {
    warning("The provided `retain_coldata_cols` are not present in any SCEs.")
  }

  # check that library id and sample id are present in metadata
id_checks <- sce_list |>
  purrr::map(\(sce){
    all(c("library_id", "sample_id") %in% names(metadata(sce)))
  }) |>
  unlist()
  
if (!all(id_checks)){
  stop("The metadata for each SCE object must contain `library_id` and `sample_id`.")
}



  # Prepare SCEs
  sce_list <- sce_list |>
    purrr::imap(prepare_sce_for_merge,
      batch_column = batch_column,
      cell_id_column = cell_id_column,
      shared_features = shared_features,
      retain_coldata_cols = retain_coldata_cols,
      preserve_rowdata_cols = preserve_rowdata_cols
    )

  # Create the merged SCE from the processed list ------------------
  merged_sce <- do.call(cbind, sce_list)

  # combine library metadata into a single data frame
  metadata_list <- metadata(merged_sce)

  # create a single list of library metadata
  library_metadata <- metadata_list[names(metadata_list) == "library_metadata"]

  # get a single list of all library ids
  library_ids <- metadata(merged_sce)[names(metadata_list) == "library_id"] |>
    unlist() |>
    unique()

  # get a single list of all sample ids
  sample_ids <- metadata(merged_sce)[names(metadata_list) == "sample_id"] |>
    unlist() |>
    unique()

  # combine into one metadata list
  metadata_list <- list(
    library_id = library_ids,
    sample_id = sample_ids,
    library_metadata = library_metadata
  )

  # if object has sample metadata then combine into a single data frame
  if("sample_metadata" %in% names(metadata_list)){
    sample_metadata <- metadata_list[names(metadata_list) == "sample_metadata"] |>
      purrr::map(as.data.frame) |>
      dplyr::bind_rows() |>
      unique()

    metadata_list[["sample_metadata"]] <- sample_metadata
  }

  metadata(merged_sce) <- metadata_list

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
  sce <- sce[shared_features, ]

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
    # Create the missing column only if it should be retained
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

  # separate library and sample metadata
  metadata_list <- metadata(sce)
  library_metadata <- metadata_list[names(metadata_list) != "sample_metadata"]
  sample_metadata <- metadata_list$sample_metadata

  # combine into one list
  metadata_list <- list(
    library_id = metadata(sce)[["library_id"]],
    sample_id = metadata(sce)[["sample_id"]],
    library_metadata = library_metadata,
    sample_metadata = sample_metadata
  )

  # replace existing metadata
  metadata(sce) <- metadata_list

  # return the processed SCE
  return(sce)
}
