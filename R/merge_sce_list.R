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
      "barcodes"
    ),
    preserve_rowdata_cols = NULL,
    cell_id_column = "cell_id",
    include_altexp = TRUE) {

  ## Checks --------------------------
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

  # Check `retain_coldata_cols`
  if (length(retain_coldata_cols) == 0) {
    warning("All pre-existing colData will be removed from the the merged SCE.
     Please check that `retain_coldata_cols` was correctly specified.")
  }

  # Check altExp compatibility, if we are including them
  # altExps for a given name must all have the same features and the same assays
  if (include_altexp) {

    # This is a list of lists of altexp information for later use:
    # (altexp_name = list(  features = c(features), assays = c(assays) ))
    altexp_attributes <- check_altexps(sce_list)

  }


  ## Subset SCEs to shared features and ensure appropriate naming ------------------

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

  ## Prepare main experiment of SCEs for merging --------------------
  sce_list <- sce_list |>
    purrr::imap(
      prepare_sce_for_merge,
      batch_column = batch_column,
      cell_id_column = cell_id_column,
      shared_features = shared_features,
      retain_coldata_cols = retain_coldata_cols,
      preserve_rowdata_cols = preserve_rowdata_cols
    )


  ## Handle metadata ---------------------------------------------
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


  ## Handle altExps ------------------------------------------------------

  # If we are including altExps, process them and save to list to add to merged SCE
  merged_altexps <- list()
  if (include_altexp) {

    # First we need to determine the final column names of the merged_sce (not yet made)
    #  for use in altExp code. Later we'll apply this order to the merged_sce itself.
    # These values are cell ids: `{sce_name}-{barcode}`
    all_merged_barcodes <- sce_list |>
      purrr::map(colnames) |>
      purrr::reduce(union)

    # Find all altExp names present in the SCE objects.
    # We will prepare a merged altExp for each of these.
    altexp_names <- sce_list |>
      purrr::map(
        \(sce) altExpNames(sce)
      ) |>
      purrr::reduce(union)

    for (altexp_name in altexp_names) {
      # Determine which SCEs contain this altExp, and create list of those altExps
      altexp_list <- sce_list |>
        purrr::keep(\(sce) altexp_name %in% altExpNames(sce)) |>
        purrr::map(altExp)

      # Create and save the merged altExp for this altexp_name
      merged_altexps[[altexp_name]] <- create_merged_altexp(
        altexp_list,
        all_merged_barcodes
      )

    }
  }

  # Remove altExps from SCEs prior to main experiment merge
  # If none are present, this code has no effect.
  sce_list <- sce_list |>
    purrr::map(removeAltExps)

  # Create the merged SCE from the processed list
  merged_sce <- do.call(cbind, sce_list)

  # Replace existing metadata list with merged metadata
  metadata(merged_sce) <- metadata_list

  # Add the merged altE into the main merged_sce
  if (include_altexp) {

    # Ensure compatible column names
    # (this is probably not necessary but doesn't hurt...)
    merged_sce <- merged_sce[,all_merged_barcodes]

    # Add the merged altexps into the merged sce
    for (altexp_name in names(merged_altexps)) {
      altExp(merged_sce, altexp_name) <- merged_altexps[[altexp_name]]
    }
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
  missing_columns <- setdiff(
    retain_coldata_cols,
    observed_coldata_names
  )
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
    sample_metadata = sample_metadata # this will be the same as the previous sample_metadata
  )

  # replace existing metadata
  metadata(sce) <- metadata_list

  # return the processed SCE
  return(sce)
}



#' Prepare altExps for merge and create a list of merged altExps for each altExp name
#'
#'
#' @param altexp_list List of altexps to merge
#' @param all_merged_barcodes Vector of column names (`{sce_name}-{barcode}`) to include
#'   in the final merged altExp. This vector includes _all_ SCEs, not only those
#'   with this altExp name.
#'
#' @return A list of merged altExps to include the final merged SCE object
create_merged_altexp <- function(
    altexp_list,
    all_merged_barcodes) {

  # Create vector of all features
  # this order will be used for the final assay matrix/ces
  altexp_features <- altexp_list |>
    purrr::map(rownames) |>
    purrr::reduce(union)

  # Determine which assays are present for this altexp_name. We'll need a matrix
  # for each of these in the final merged object
  altexp_assay_names <- altexp_list |>
    purrr::map(assayNames) |>
    purrr::reduce(union)

  # Create new merged assay matrices
  new_assays <- altexp_assay_names |>
    purrr::map(
      build_new_altexp_assay,
      altexp_list,
      altexp_features,
      all_merged_barcodes
    )
  names(new_assays) <- altexp_assay_names


  # Create merged altExp
  merged_altexp <- SingleCellExperiment(assays = new_assays)

  # TODO: Add rowData and colData to merged_altexp

  return(merged_altexp)

}




#' Build a new sparse matrix for merging altExps
#'
#' @param assay_name Name of assay of interest (e.g., "counts")
#' @param altexp_list List of altExps which should be included in the new matrix
#' @param all_merged_features Vector of matrix row names, corresponding to the full
#'   set of features for this altExp
#' @param all_merged_barcodes Vector of matrix column names, corresponding to all cells
#'   which will be in the final merged altExp
#'
#' @return Sparse matrix
build_new_altexp_assay <- function(
    assay_name,
    altexp_list,
    all_merged_features,
    all_merged_barcodes) {

  # Establish new matrix with all NA values
  new_matrix <- matrix(
    data = NA,
    nrow = length(all_merged_features),
    ncol = length(all_merged_barcodes),
    dimnames = list(
      all_merged_features,
      all_merged_barcodes
    )
  )

  # Substitute existing assays into the matrix, if they exist
  for (altexp in altexp_list) {
    # Note that column names were already formatted as `{sce_name}-{barcode}` by
    #  the main SCE merging code
    if (assay_name %in% assayNames(altexp)) {
      # as.matrix() is needed here
      new_matrix[rownames(altexp), colnames(altexp)] <- as.matrix( assay(altexp, assay_name) )
    }

  }
  # sparsify
  new_matrix <- as(new_matrix, "CsparseMatrix")

  return(new_matrix)
}



#' Helper function to check altExp compatibility
#'
#' @param sce_list List of SCEs with altExps to check
#'
#' @return List of named lists with altExp information for use when preparing to merge,
#'   with each sublist formatted as:
#'   altexp_name = list(features = c(features), assays = c(assays))
check_altexps <- function(sce_list) {

  # Attribute list to save for later use
  altexp_attributes <- list()

  # Find all altExp names present in the SCE objects.
  altexp_names <- sce_list |>
    purrr::map(
      \(sce) altExpNames(sce)
    ) |>
    purrr::reduce(union)

  # For each in altexp_names (if present), do they have the same features?
  # If not, error out
  for (altexp_name in altexp_names) {

    # all altExps for this name
    altexp_list <- sce_list |>
      purrr::keep(\(sce) altexp_name %in% altExpNames(sce)) |>
      purrr::map(altExp, altexp_name)

    # find their union of features
    all_features <- altexp_list |>
      purrr::map(rownames) |>
      purrr::reduce(union) |>
      sort()

    # create logical vector for presence of all features
    features_present <- altexp_list |>
      purrr::map_lgl(
        \(alt_sce) identical(all_features, sort(rownames(alt_sce)))
      )

    if (!all(features_present)) {
      stop(
        glue::glue("The {altexp_name} alternative experiments do not share the same set of features.")
      )
    }

    # check for same assays
    all_assays <- altexp_list |>
      purrr::map(assayNames) |>
      purrr::reduce(union)|>
      sort()

    # create logical vector for presence of all assays
    assays_present <- altexp_list |>
      purrr::map_lgl(
        \(alt_sce) identical(all_assays, sort(assayNames(alt_sce)))
      )

    # TODO: we may want to drop assays that aren't present in all altexps, rather than dying.
    if (!all(assays_present)) {
      stop(
        glue::glue("The {altexp_name} alternative experiments do not share the same set of assays.")
      )
    }

    # Save to altexp_attributes for later use
    altexp_attributes[[altexp_name]] <- list(
      "features" = all_features,
      "assays"   = all_assays
    )

  }
  return(altexp_attributes)
}
