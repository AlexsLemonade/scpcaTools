#' Merge a list of SCEs as preparation for formal integration
#'
#' This function takes a named (ideally by a form of library ID) list of
#'  SingleCellExperiment (SCE) objects and merges them into one SCE object.
#'  All names in the colData slot must match across provided SCE objects.
#'  The resulting colData slot will contain information in the given batch column
#'  to differentiate originating SCE objects. Only features (genes) that are present
#'  in all provided SCE objects will be retained. To preserve original feature-level
#'  information, all final rowData slot column names will be appended with the
#'  given SCE's name, as `<column_name>-<sce_name>`.
#'
#' @param sce_list A named list of SingleCellExperiment objects
#' @param batch_column A character value giving the resulting colData column name
#'  to differentiate originating SingleCellExperiment objects Default value is `batch`.
#' @param preserve_rowdata_columns An array of columns that appear in originating
#'  SingleCellExperiment object's rowData slot which should not be renamed with
#'  the given SingleCellExperiment objet name. These are generally columns which are
#'  not specific to the given library's preparation or statistics. For example,
#'  such an array might contain items like "Gene", "ensembl_ids", "gene_symbol", etc.
#'
#' @return A SingleCellExperiment containing all SingleCellExperiment objects present
#'   in the inputted list
#' @export
#'
#' @import SingleCellExperiment
#'
merge_sce_list <- function(sce_list = list(),
                           batch_column = "batch",
                           preserve_rowdata_columns = NULL) {
  # Ensure `sce_list` is named (according to library IDs) ----------------------
  if (is.null(names(sce_list))) {
    stop("Individual SingleCellExperiment objects in `sce_list` must be named.")
  }
  # Ensure `sce_list` has >=2 items --------------------------------------------
  if (length(sce_list) < 2) {
    stop("The `sce_list` must contain 2 or more SingleCellExperiment objects to merge.")
  }

  # Check that colnames of colData match across all SCE objects ---------------

  # Create of colData column names for each SCE
  sce_colnames_list <- purrr::map(sce_list, ~colnames(colData(.)))

  # Compare each list item back to the first set of colnames, and throw an error
  #  if column names do not match
  check_name_diffs <- function(test_list) {
    # Compare test_list back to the first namelist.
    name_diffs <- setdiff(sce_colnames_list[[1]], test_list)
    if (length(name_diffs) != 0) {
      stop("Error: All SingleCellExperiment objects in the provided list must have
             the same colData slot column names.")
    }
  }
  purrr::walk(
    sce_colnames_list[2:length(sce_colnames_list)],
    check_name_diffs
  )

  # TODO: Alternatively, we can also add a "dummy" column here with all NAs for any "loner columns"
  #  We'd throw a warning when this happens.

  # Subset SCEs to shared features and ensure appropriate naming ------------------

  # First, obtain intersection among all SCE objects
  shared_features <- sce_list |>
    purrr::map(rownames) |>
    purrr::reduce(intersect)

  # Now, prepare SCEs by subsetting to shared features, updating column names, and
  #   adding the batch column
  prepare_sce_for_merge <- function(sce, sce_name) {

    # TODO: Are we stilll doing this?
    # remove any alternative experiments before merging
    #alt_names <- altExpNames(sce)
    #if(!is.null(alt_names)){
    #  altExp(sce) <- NULL
    #}

    # Subset to shared features
    sce <- sce[shared_features,]

    # Add `sce_name` to colData row names so cell barcodes can be mapped to originating SCE
    colnames(sce) <- glue::glue("{colnames(sce)}-{sce_name}")

    # Add in `batch_column` with `sce_name` to track originating SCE
    sce[[batch_column]] <- sce_name

    # Add `sce_name` ID to rowData column names except for those present in `preserve_rowdata_columns`
    original_colnames <- colnames(rowData(sce))
    colnames(rowData(sce)) <- ifelse(
      original_colnames %in% preserve_rowdata_columns,
      original_colnames,
      glue::glue("{original_colnames}-{sce_name}")
    )

    # return the processed SCE
    return(sce)
  }

  sce_list <- purrr::map2(sce_list,
                          names(sce_list),
                          prepare_sce_for_merge)


  # Create the merged SCE from the processed list ------------------
  merged_sce <- do.call(cbind, sce_list)


  # TODO: ARE WE STILL DOING THIS?
  # Retain only colData names that are the `batch_column` or columns added by
  #  scuttle::addPerCellQC() during the scpca-downstream-analyses processing
  #mito_names <- names(colData(merged_sce)) %>%
  #  stringr::str_subset("^subsets_mito")
#
#  retain_cols <-  c(batch_column,
#                    "celltype",
#                    # scuttle::addPerCellQC() columns
#                    "sum", "detected", mito_names)
#  retain_pos <- which(names(colData(merged_sce)) %in% retain_cols)
#  colData(merged_sce) <- colData(merged_sce)[, retain_pos]


  # Return merged SCE object ----------------------------
  return(merged_sce)

}
