#' Merge a list of SCEs into one SCE object
#'
#' This function takes an optionally-named (if named, ideally by a form of
#'  library ID) list of SingleCellExperiment (SCE) objects and merges them into
#'  one SCE object. At least some genes must be present in all SCEs in order to
#'  merge them. By default, alternative experiments (altExps) are retained in the
#'  final merged object, but each altExp of a given name is required to have identical
#'  features.
#'
#'  Original SCE contents are modified or retained as follows:
#'  - The resulting colData slot will include a new column specified by
#'    `batch_column` (default "library_id") that either holds the originating SCE
#'    object's name (referred to as `sce_name` here), or if it is unnamed then
#'    its index in the provided `sce_list`.
#'  - The resulting colData slot will include another new column `cell_id_column`
#'    (default "cell_id") that will contain the SCE's original column names (i.e.
#'    original colData row names). Often, but not always, this row name holds a
#'    unique cell barcode.
#'  - Of the original colData columns, only column names provided in the argument
#'    `retain_coldata_cols` will be retained.
#'  - The resulting colData rownames will be be prefixed with `{sce_name-}`.
#'  - The resulting rowData slot column names will be appended with the given
#'    SCE's name, as `{sce_name}-{column_name}` except for columns whose names
#'    are indicated to preserve with the `preserve_rowdata_cols` argument.
#'
#'  SCE altExp contents are modified or retained as follows:
#'  - As with the main experiment, the additional columns `batch_column` and
#'  `cell_id_column` will be added to the colData slot.
#'  - Of the original altExp colData columns, only column names provided in
#'    the argument `retain_altexp_coldata_cols`, as specified for each named
#'    altExp will be retained.
#'  - The resulting rowData slot column names will be appended with the given
#'    SCE's name, as `{sce_name}-{column_name}` except for the column `target_type`
#'    which is commonly present in CITE-seq alternative experiments.
#'
#'
#' @param sce_list A list of SingleCellExperiment objects. The list may optionally
#'  be named with batch information. If no names are provided, names will be generated
#'  based on the SCE's index.
#' @param batch_column A character value giving the resulting colData column name
#'  to differentiate originating SingleCellExperiment objects. Often these values
#'  are unique library IDs. Default value is `"library_id"`.
#' @param cell_id_column A character value giving the resulting colData column name
#'  to hold unique cell IDs formatted as their original row name. Default
#'  value is `"cell_id"`.
#' @param retain_coldata_cols A vector of colData columns which should be retained
#'  in the the final merged SCE object. If columns are missing from any SCE to be merged,
#'  they will be created and populated with `NA` values. A vector of default columns
#'  to retain is given in the function definition.
#' @param include_altexp Boolean for whether altExps, if present, should be
#' included in the final merged object. Default is `TRUE`.
#' @param preserve_rowdata_cols A vector of column names that appear in originating
#'  SCE objects' rowData slots which should not be renamed with
#'  the given SCE object name or index is name is not given. These are generally
#'  columns which are not specific to the given library's preparation or statistics.
#'  For example, such a vector might contain items like "Gene", "ensembl_ids", etc.
#'  Default value is `NULL`.
#' @param retain_altexp_coldata_cols Named list containing vectors of column names
#'  that should be retained in altExp colData. Elements are named by the altExp
#'  for which the columns should be retained. If any given altExp name is not
#'  present in any SCE, it will be ignored. If columns are missing from any given
#'  altExp to be merged, they will be created and populated with `NA` values.
#'  Default value is `NULL`.
#'
#' @return A SingleCellExperiment object containing all SingleCellExperiment objects
#'   present in the inputted list
#' @export
#'
#' @examples
#' \dontrun{
#' # Merge list of SCEs, specifying a different batch column name
#' merge_sce_list(
#'   sce_list = list("sce1" = sce1, "sce2" = sce2),
#'   batch_column = "batch"
#' )
#'
#'
#' #' # Merge list of SCEs but do include any alternative experiments in the merged object
#' merge_sce_list(
#'   sce_list = list("sce1" = sce1, "sce2" = sce2),
#'   include_altexp = FALSE
#' )
#'
#'
#' # Merge list of SCEs and include alternative experiments in the merged object.
#' # The provided list of SCEs may contain alternative experiments named `"adt"`
#' #  and/or `"other_altexp"` which, if present, are expected to respectively
#' #  have the columns shown in the `retain_altexp_coldata_cols`.
#' merge_sce_list(
#'   sce_list = list("sce1" = sce1, "sce2" = sce2),
#'   # columns to retain in the given alternative experiment, if it is present
#'   retain_altexp_coldata_cols = list(
#'     "adt" = c("discard", "high.controls"),
#'     "other_altexp" = c("first_column", "second_column")
#'   )
#' )
#' }
#'
#' @import SingleCellExperiment
merge_sce_list <- function(
    sce_list = list(),
    batch_column = "library_id",
    cell_id_column = "cell_id",
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
    include_altexp = TRUE,
    preserve_rowdata_cols = NULL,
    retain_altexp_coldata_cols = NULL) {
  if (is.null(names(sce_list))) {
    warning(
      glue::glue(
        "Individual SCE objects in `sce_list` are not named, so batches will be
        named based on their list index in the merged SCE object."
      )
    )
    names(sce_list) <- seq_along(sce_list)
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

  # Check metadata - library_id and sample_id must be present
  check_metadata(sce_list)

  # Check altExp compatibility, if we are including them
  # altExps for a given name must all have the same features and the same assays
  if (include_altexp) {
    # This is a list of lists of altexp information for later use:
    # (altexp_name = list(  features = c(features), assays = c(assays) ))
    altexp_attributes <- get_altexp_attributes(sce_list)

    # Define more main SCE colData columns to keep, based on altExp names
    coldata_suffixes <- c("sum", "detected", "percent")
    altexp_columns <- glue::glue("altexps_{names(altexp_attributes)}") |>
      purrr::map(
        \(prefix) {
          glue::glue("{prefix}_{coldata_suffixes}")
        }
      ) |>
      unlist()

    # Update retain_coldata_cols
    retain_coldata_cols <- c(retain_coldata_cols, altexp_columns)

    # Check each altExp metadata, for SCEs with the altExp
    names(altexp_attributes) |>
      purrr::walk(
        \(altexp_name) {
          sce_list |>
            purrr::keep(\(sce) altexp_name %in% altExpNames(sce)) |>
            purrr::map(altExp, altexp_name) |>
            check_metadata()
        }
      )
  } else {
    # Remove altexps if we are not including them
    sce_list <- sce_list |>
      purrr::map(removeAltExps)
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


  ## Handle altExps ------------------------------------------------------
  # If we are including altExps, process them and save to list to add to merged SCE
  if (include_altexp) {
    for (altexp_name in names(altexp_attributes)) {
      expected_assays <- altexp_attributes[[altexp_name]][["assays"]]
      expected_features <- altexp_attributes[[altexp_name]][["features"]]
      altexp_coldata_cols <- purrr::pluck(retain_altexp_coldata_cols, altexp_name)

      sce_list <- sce_list |>
        purrr::imap(
          prepare_altexp_for_merge,
          altexp_name,
          expected_assays,
          expected_features,
          batch_column = batch_column,
          cell_id_column = cell_id_column,
          retain_altexp_coldata_cols = altexp_coldata_cols
        )
    }
  }

  # Create the merged SCE from the processed list
  merged_sce <- do.call(combineCols, unname(sce_list))

  # Update metadata in merged objects, using the unmerged sce_list
  metadata(merged_sce) <- sce_list |>
    purrr::map(metadata) |>
    prepare_merged_metadata()

  if (include_altexp) {
    for (altexp_name in names(altexp_attributes)) {
      has_altexp_name <- sce_list |>
        purrr::map_lgl(\(sce) (altexp_name %in% altExpNames(sce)))

      # For any SCEs without this altExp, create the library_id and sample_id metadata
      additional_metadata <- sce_list |>
        purrr::discard(has_altexp_name) |>
        purrr::map(extract_metadata_for_altexp)

      # Update metadata in altExps that were originally present
      altexp_metadata_list <- sce_list |>
        purrr::keep(has_altexp_name) |>
        purrr::map(altExp, altexp_name) |>
        purrr::map(metadata) |>
        # Tack on the metadata we created for libraries without altExps
        c(additional_metadata)

      # Ensure correct order
      altexp_metadata_list <- altexp_metadata_list[names(sce_list)]

      metadata(altExp(merged_sce, altexp_name)) <- prepare_merged_metadata(altexp_metadata_list)
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
#' @param is_altexp Boolean if we are preparing an altExp or not.
#'   Default is `FALSE`. If FALSE, the SCE colnames will not be modified.
#'
#' @return An updated SCE that is prepared for merging
prepare_sce_for_merge <- function(
    sce,
    sce_name,
    batch_column,
    cell_id_column,
    shared_features,
    retain_coldata_cols,
    preserve_rowdata_cols,
    is_altexp = FALSE) {
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

  # Only modify column names if we are not working with an altExp.
  # This avoids having double `{sce_name}-{sce_name}` prefixes
  if (!is_altexp) {
    # Add `sce_name` to colnames so cell ids can be mapped to originating SCE
    colnames(sce) <- glue::glue("{sce_name}-{colnames(sce)}")
  }

  # return the processed SCE
  return(sce)
}



#' Prepare an alternative experiment for merging
#'
#' If the given `altexp_name` name is missing from the SCE, a new one with all
#'   NA assays will be created.
#' @param sce SCE with altExp to prepare
#' @param sce_name SCE name
#' @param altexp_name Name of altExp of interest
#' @param expected_assays Vector of assays that should be in each altExp
#' @param expected_features Vector of features that should be in each altExp
#' @param batch_column The name of the batch column which will be added to the
#'   colData slot.
#' @param cell_id_column The name of the cell_id column which will be added to the
#'   colData slot.
#' @param retain_altexp_coldata_cols Named list of columns that should be retained
#'  in alternative experiment colData. Each name should correspond to the alternative
#'  experiment in which it should be retained.
#' @param preserve_rowdata_cols altExp rowData columns which should not be renamed
#'
#' @return An SCE with an updated altExp
prepare_altexp_for_merge <- function(
    sce,
    sce_name,
    altexp_name,
    expected_assays,
    expected_features,
    batch_column,
    cell_id_column,
    retain_altexp_coldata_cols,
    preserve_rowdata_cols = c("target_type")) {
  if (!altexp_name %in% altExpNames(sce)) {
    return(sce)
  }

  # Now, prepare this altexp for merge
  altExp(sce, altexp_name) <- prepare_sce_for_merge(
    altExp(sce, altexp_name),
    sce_name,
    batch_column = batch_column,
    cell_id_column = cell_id_column,
    shared_features = expected_features,
    retain_coldata_cols = retain_altexp_coldata_cols,
    preserve_rowdata_cols = preserve_rowdata_cols,
    is_altexp = TRUE
  )

  return(sce)
}


#' Helper function to check altExp compatibility
#'
#' @param sce_list List of SCEs with altExps to check
#'
#' @return List of named lists with altExp information for use when preparing to merge,
#'   with each sublist formatted as:
#'   altexp_name = list(features = c(features), assays = c(assays))
get_altexp_attributes <- function(sce_list) {
  # Attribute list to save for later use
  altexp_attributes <- list()

  # Find all altExp names present in the SCE objects.
  altexp_names <- sce_list |>
    purrr::map(altExpNames) |>
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
      purrr::reduce(union)

    # create logical vector for presence of all features
    features_present <- altexp_list |>
      purrr::map_lgl(
        \(alt_sce) setequal(rownames(alt_sce), all_features)
      )

    if (!all(features_present)) {
      stop(
        glue::glue("The {altexp_name} alternative experiments do not share the same set of features.")
      )
    }

    # check for same assays
    all_assays <- altexp_list |>
      purrr::map(assayNames) |>
      purrr::reduce(union)

    # Save to altexp_attributes for later use
    altexp_attributes[[altexp_name]] <- list(
      "features" = all_features,
      "assays"   = all_assays
    )
  }
  return(altexp_attributes)
}





#' Prepare an updated list of metadata for the merged SCE
#'
#' The metadata will be reformatted to contain the following fields:
#' - `library_id`:       vector of all library ids in the merged object
#' - `sample_id`:        vector of all sample ids in the merged object
#' - `library_metadata`: pre-merge metadata list for each library, named by SCE name
#' - `sample_metadata`:  a data frame with all distinct rows from the pre-merge sample_metadata data frames,
#'   with all values coerced to character
#'
#' @param metadata_list List of metadata to update
#'
#' @return Updated metadata list to include in the merged SCE
prepare_merged_metadata <- function(metadata_list) {
  # Define vectors of library and sample ids
  metadata_library_ids <- metadata_list |>
    purrr::map_chr("library_id") |>
    unname()

  metadata_sample_ids <- metadata_list |>
    # can't use map_chr in case we have multiplexed
    purrr::map("sample_id") |>
    unlist() |>
    unname()

  # Grab names to check contents
  transposed_names <- metadata_list |>
    purrr::map(names) |>
    purrr::reduce(union)

  # first check that this library hasn't already been merged
  if ("library_metadata" %in% transposed_names) {
    stop(paste(
      "This SCE object appears to be a merged object",
      "We do not support merging objects with objects that have already been merged."
    ))
  }

  # Create new version of metadata_list without any sample_metadata fields
  # Note that this will work even if `sample_metadata` is not present
  library_metadata <- metadata_list |>
    purrr::map(
      \(meta) meta[which(names(meta) != "sample_metadata")]
    )

  # combine into final metadata list
  final_metadata_list <- list(
    # vector of all library ids
    library_id = metadata_library_ids,
    # vector of all sample ids
    sample_id = metadata_sample_ids,
    # list of existing metadata for each library, with sample_metadata removed
    library_metadata = library_metadata
  )

  # if object has sample metadata then combine into a single data frame,
  #  and add this to the final_metadata_list
  if ("sample_metadata" %in% transposed_names) {
    sample_metadata <- metadata_list |>
      purrr::map("sample_metadata") |>
      purrr::map(\(df) {
        # make sure all column types are compatible first
        df |>
          dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
      }) |>
      dplyr::bind_rows() |>
      dplyr::distinct()

    # separate into individual sample ids for the check below
    # to help checking cellhash samples w/ multiple sample_ids per library_id
    all_sample_ids <- sample_metadata$sample_id |>
      purrr::map(\(x) stringr::str_split_1(x, pattern = ",")) |>
      unlist()

    # check that all sample ids are found in the new sample metadata and warn if not
    if (!all(metadata_sample_ids %in% all_sample_ids)) {
      warning("Not all sample ids are present in metadata(merged_sce)$sample_metadata.")
    }

    # add to final_metadata_list
    final_metadata_list$sample_metadata <- sample_metadata
  }

  return(final_metadata_list)
}



#' Helper function to check that expected metadata fields are present in a given
#' list of SCEs.
#'
#' The default expected fields are `library_id` and `sample_id`. If either is missing from
#'   any SCE, an error is thrown. This function does not return anything.
#'
#' @param sce_list List of SCEs to check
#' @param expected_fields a vector of metadata fields that should be present
check_metadata <- function(sce_list, expected_fields = c("library_id", "sample_id")) {
  metadata_checks <- sce_list |>
    purrr::map_lgl(\(sce) {
      all(expected_fields %in% names(metadata(sce)))
    })

  if (!all(metadata_checks)) {
    stop(glue::glue("
       The metadata for each SCE object must contain {stringr::str_flatten_comma(expected_fields, ', and ')}.
       If `include_altexp` is TRUE, these fields must also be present in all altExps.
     "))
  }
}

#' Helper function to extract main SCE metadata for inclusion in an altExp
#'
#' @param sce SCE object to extract metadata from
#'
#' @return List with fields `library_id` and `sample_id`
extract_metadata_for_altexp <- function(sce) {
  list(
    library_id = metadata(sce)$library_id,
    sample_id = metadata(sce)$sample_id
  )
}
