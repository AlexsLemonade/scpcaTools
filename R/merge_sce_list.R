#' Merge a list of SCEs into one SCE object
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
#' @param include_altexp Boolean for whether or not any present alternative experiments
#'  should be included in the final merged object. Default is TRUE.
#'
#' @return A SingleCellExperiment object containing all SingleCellExperiment objects
#'   present in the inputted list
#' @export
#'
#' @import SingleCellExperiment
merge_sce_list <- function(
    sce_list = list(),
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
    cell_id_column = "cell_id",
    include_altexp = TRUE) {

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

  # If we are including altExps, extract them now for separate merging, to be
  #  added back into the final merged SCE at the end
  if (include_altexp) {
    altexp_list <- sce_list |>
      purrr::map(altExp)
  }
  sce_list <- sce_list |>
    purrr::map(removeAltExps)


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
  all_colnames <- sce_list |>
    purrr::map(
      \(sce) names(colData(sce))
    ) |>
    unlist() |>
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

  if (!all(id_checks)) {
    stop("The metadata for each SCE object must contain `library_id` and `sample_id`.")
  }

  # Prepare main experiment of SCEs for merging
  sce_list <- sce_list |>
    purrr::imap(
      prepare_sce_for_merge,
      batch_column = batch_column,
      cell_id_column = cell_id_column,
      shared_features = shared_features,
      retain_coldata_cols = retain_coldata_cols,
      preserve_rowdata_cols = preserve_rowdata_cols
    )

  # get a list of metadata from the list of sce objects
  # each library becomes an element within the metadata components
  metadata_list <- sce_list |>
    purrr::map(metadata) |>
    purrr::transpose()

  # flatten and reduce the library ids and sample ids
  metadata_list$library_id <- metadata_list$library_id |>
    unlist() |>
    unique()

  metadata_list$sample_id <- metadata_list$sample_id |>
    unlist() |>
    unique()

  # if object has sample metadata then combine into a single data frame
  if ("sample_metadata" %in% names(metadata_list)) {
    sample_metadata <- metadata_list$sample_metadata |>
      dplyr::bind_rows() |>
      unique()

    # check that all sample ids are found in the new sample metadata and warn if not
    if (!all(metadata_list$sample_id %in% sample_metadata$sample_id)) {
      warning("Not all sample ids are present in metadata(merged_sce)$sample_metadata.")
    }

    # replace sample metadata in metadata list
    metadata_list$sample_metadata <- sample_metadata
  }

  # Create the merged SCE from the processed list and replace existing metadata list with merged metadata
  merged_sce <- do.call(cbind, sce_list)
  metadata(merged_sce) <- metadata_list

  # If we are including altExps, process them and add to the merged sce
  if (include_altexp) {

    # Find all shared features
    altexp_features <- sce_list |>
      purrr::map(
        \(sce) rownames(altExp(sce))
      ) |>
      purrr::reduce(union)

    # Add NA values for features where needed, and otherwise prepare for merging
    altexp_list <- altexp_list |>
      purrr::imap(
        prepare_altexps_for_merge,
        altexp_features
      )

    # merge altExps
    merged_altexps <- do.call(cbind, altexp_list)

    # add the merged altExp to the merged_sce
    # arbitrarily select the first altExp name to use here, since they are all the same
    altExp(merged_sce, altExpNames(sce_list[[1]])) <- merged_altexps
  }

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
prepare_sce_for_merge <- function(
  sce,
  sce_name,
  batch_column,
  cell_id_column,
  shared_features,
  retain_coldata_cols,
  preserve_rowdata_cols) {

  # Subset to shared features
  sce <- sce[shared_features, ]

  ##### rowData #####
  # Add `sce_name` ID to rowData column names except for those
  #  present in `preserve_rowdata_cols`
  original_colnames <- colnames(rowData(sce))

  # only rename rowData if there is any rowData to begin with
  if (length(original_colnames) > 0) {
    colnames(rowData(sce)) <- ifelse(
      original_colnames %in% preserve_rowdata_cols,
      original_colnames,
      glue::glue("{sce_name}-{original_colnames}")
    )
  }

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

  # get metadata list for updating it
  metadata_list <- metadata(sce)

  # first check that this library hasn't already been merged
  if ("library_metadata" %in% names(metadata_list)) {
    stop("This SCE object appears to be a merged object. We do not support merging objects with objects that have already been merged.")
  }

  # create library and sample metadata.
  # library metadata will hold all the previous metadata fields, to avoid conflicts
  library_metadata <- metadata_list[names(metadata_list) != "sample_metadata"]
  sample_metadata <- metadata_list$sample_metadata

  # combine into one list
  metadata_list <- list(
    library_id = metadata(sce)[["library_id"]],
    sample_id = metadata(sce)[["sample_id"]],
    library_metadata = library_metadata, # this will be all previous metadata for the given library
    sample_metadata = sample_metadata # this will be the same as the previous sample_metadata column
  )

  # replace existing metadata
  metadata(sce) <- metadata_list

  # return the processed SCE
  return(sce)
}



#' Prepare altExps for merge by ensuring that all altExps have the same features.
#'
#' For any features in `altexp_features` that are missing from the SCE's altExp,
#'  add those features in with `NA` counts, and update the rowData slot to match.
#'  This involves creating a new SCE to replace the current altExp.
#'
#' @param sce The SCE object whose altExp should be prepared
#' @param sce_name The name for this SCE
#' @param altexp_features Vector of features that should be present in the altExp
#'
#' @return An updated SCE object with all altExp features present
prepare_altexps_for_merge <- function(
    sce,
    sce_name,
    altexp_features) {

  sce_altexp <- altExp(sce)

  # Determine which features are missing
  missing_features <- setdiff(altexp_features, rownames(sce_altexp))
  n_missing <- length(missing_features)

  # Update altExp if any features are missing
  if (n_missing > 0) {

    # First, establish new rowData with all NA values
    # this data.frame only contains rows for missing features
    new_rowdata <- matrix(
      NA,
      nrow = n_missing,
      ncol = ncol(rowData(sce_altexp))
    ) |>
      as.data.frame() |>
      purrr::set_names(colnames(rowData(sce_altexp)))

    # create new rownames to be added in when we re-DataFrame this
    new_rownames <- c(
      rownames(sce_altexp), # existing rownames
      missing_features # new rownames
    )

    # Next, establish new assay matrices
    assay_names <- assayNames(sce_altexp)
    new_assays <- assay_names |>
      purrr::map(
        update_altexp_assay,
        sce_altexp,
        altexp_features
      )
    names(new_assays) <- assay_names

    # Create a new altexp starting with new assay matrices
    new_altexp <- SingleCellExperiment(assays = new_assays)

    # Add the new rowData rows into the SCE, with updated column names
    rowData(new_altexp) <- sce_altexp |>
      rowData() |>
      as.data.frame() |>
      rbind(new_rowdata) |>
      # re DataFrame with _all_ feature rownames
      DataFrame(row.names = new_rownames)

    # Add colData and metadata
    colData(new_altexp) <- colData(altExp(sce))
    metadata(new_altexp) <- metadata(altExp(sce))

    # Replace the old altExp, while preparing it for merge
    altExp(sce) <- prepare_sce_for_merge(
      new_altexp,
      sce_name,
      batch_column = "library_id",
      cell_id_column = "cell_id",
      shared_features = altexp_features,
      retain_coldata_cols = colnames(colData(new_altexp)),
      preserve_rowdata_cols = NA # rename them all
    )
  }

  # Return the SCE with updated altExp
  return(sce)
}



#' Create a new assay sparse matrix with `NA` values for missing features
#'
#' @param assay_name The name of the assay to update
#' @param sce_altexp Alternative experiment SCE to modify
#' @param altexp_features All features that should end up in the matrix
#'
#' @return Updated sparse matrix
update_altexp_assay <- function(
    assay_name,
    sce_altexp,
    altexp_features) {

  # pull out existing matrix
  assay_sparse_matrix <- assay(sce_altexp, assay_name)

  # define the new full matrix
  new_matrix <- matrix(
    nrow = length(altexp_features),
    ncol = ncol(assay_sparse_matrix),
    dimnames = list(altexp_features, colnames(assay_sparse_matrix))
  )
  # fill in existing features
  new_matrix[rownames(assay_sparse_matrix),] <- matrix(assay_sparse_matrix)

  # make it sparse again
  new_matrix <- as(new_matrix, "CsparseMatrix")

  return(new_matrix)

}
