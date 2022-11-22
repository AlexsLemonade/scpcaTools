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
#' @param retain_coldata_columns_main A vector of colData columns which should be retained
#'  in the main experiment of the final merged SCE object.
#' @param retain_coldata_columns_altexps A vector of colData columns which should be retained
#'  in the altExps experiments of the final merged SCE object.
#' @param retain_altexps A logical indicating whether any present altExp slots should
#'   be retained in the merged SCE. Default is `TRUE`.
#' @param preserve_rowdata_columns A vector of column names that appear in originating
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
                           retain_altexps = TRUE,
                           retain_coldata_columns_main = c("sum",
                                                           "detected",
                                                           "total",
                                                           "subsets_mito_sum",
                                                           "subsets_mito_detected",
                                                           "subsets_mito_percent",
                                                           "miQC_pass"),
                           retain_coldata_columns_altexps = c(),
                           preserve_rowdata_columns = NULL) {
  # Ensure `sce_list` is named (according to library IDs) ----------------------
  if (is.null(names(sce_list))) {
    stop("Individual SingleCellExperiment objects in `sce_list` must be named.")
  }
  # Ensure `sce_list` has >=2 items --------------------------------------------
  if (length(sce_list) < 2) {
    stop("The `sce_list` must contain 2 or more SingleCellExperiment objects to merge.")
  }

  # Subset SCEs to shared features and ensure appropriate naming ------------------

  # First, obtain intersection among all SCE objects
  shared_features <- sce_list |>
    purrr::map(rownames) |>
    purrr::reduce(intersect)


  process_columns <- function(sce, sce_name, retain_columns) {
    # helper function to rename and subset columns in an SCE object when
    # formatting for merging

    ##### rowData #####
    # Add `sce_name` ID to rowData column names except for those
    #  present in `preserve_rowdata_columns`
    original_colnames <- colnames(rowData(sce))

    colnames(rowData(sce)) <- ifelse(
      original_colnames %in% preserve_rowdata_columns,
      original_colnames,
      glue::glue("{original_colnames}-{sce_name}")
    )

    ##### colData #####
    # Retain only the columns present in `retain_columns`
    coldata_names <- names(colData(sce))
    if (!(all(retain_columns %in% coldata_names))) {
      stop("Error: Provided columns in `retain_coldata_columns` are not present
            in all SCE objects.")
    }
    colData(sce) <- colData(sce) |>
      as.data.frame() |>
      dplyr::select( dplyr::all_of(retain_columns) ) |>
      DataFrame()

    # Add batch column
    sce[[batch_column]] <- sce_name


    return(sce)
  }



  prepare_sce_for_merge <- function(sce, sce_name) {
    # helper function for preparing SCE for merging

    # If specified, remove altExps
    #if (retain_altexps == FALSE) {
      sce <- removeAltExps(sce)
    #}

    # Subset to shared features
    sce <- sce[shared_features,]

    # Add `sce_name` to colData row names so cell barcodes can be mapped to originating SCE
    # Note that this will also affect altExps
    colnames(sce) <- glue::glue("{colnames(sce)}-{sce_name}")

    # Rename rowData and subset colData columns for both main and any present altExps
    sce <- process_columns(sce, sce_name, retain_coldata_columns_main)
    for (altexp_name in altExpNames(sce)) {
      altExp(sce, altexp_name) <- process_columns( altExp(sce, altexp_name),
                                                   sce_name,
                                                   retain_coldata_columns_altexps)
    }

    # return the processed SCE
    return(sce)
  }


  sce_list <- purrr::map2(sce_list,
                          names(sce_list),
                          prepare_sce_for_merge)


  # Create the merged SCE from the processed list ------------------
  # TODO: Fails when altExps are present
  merged_sce <- do.call(cbind, sce_list)


  # Return merged SCE object ----------------------------
  return(merged_sce)

}
