# This function can perform the crux integration based on the given integration method


#' Integrate a combined set of SingleCellExperiment objects using a specified
#'  integration method. The final SCE object will contain an additional reducedDim
#'  entitled `{integration_method}_PCA`.
#'
#' @param combined_sce A combined SCE object as prepared by `scpcaTools::merge_sce_list()`.
#' @param integration_method The integration method to use. One of `fastMNN` or
#'  `harmony` (case-insensitive)
#' @param batch_column The column in the combined SCE object indicating batches
#'  being integrated. Default is "sample".
#' @param retain_uncorrected_assays A character vector of uncorrected assays to retain in
#'   the integrated SCE object. By default, both "counts" and "logcounts" are retained.
#' @param return_corrected_expression A boolean indicating whether corrected expression
#'   values determined by the given integration method should be included in the
#'   integrated SCE object. Note that `harmony` does not calculate corrected expression,
#'   so this argument is ignored for this integration method. Default is `TRUE`.
#' @param harmony_do_PCA A `harmony`-specific argument (ignored otherwise) indicating whether
#'   integration should perform PCA. If `FALSE` (default), `harmony` will use the existing PCs to perform
#'   integration. If `TRUE`, `harmony` will will start from normalized expression when integrating and
#'   re-calculate PCs.
#' @param harmony_covariate_cols A `harmony`-specific argument (ignored otherwise)
#'   providing other columns in the combined SCE besides `batch_column` that `harmony`
#'   should consider as covariates during integration.
#' @param seed Random seed to set for integration. This is only set if the value is not `NULL`.
#' @param ... Any additional parameters to be passed to the given integration method
#'
#' @return An update combined SCE object containing
#' @export
integrate_sces <- function(combined_sce,
                           integration_method,
                           batch_column = "sample",
                           retain_uncorrected_assays = c("counts", "logcounts"),
                           return_corrected_expression = TRUE,
                           harmony_do_PCA = FALSE,
                           harmony_covariate_cols = c(),
                           seed = NULL,
                           ...) {

  # Set seed
  if (!(is.null(seed))) {
    set.seed(seed)
  }

  # Check integration_method
  integration_method <- tolower(integration_method)
  allowed_methods <- c("fastMNN", "harmony")
  if (!(integration_method %in% tolower(allowed_methods))) {
    stop("Integration method must be one of `fastMNN` or `harmony` (case-insensitive.")
  }

  # Check batch_column
  if (!(batch_column %in% names(colData(combined_sce)))) {
    stop("The provided `batch_column` column must be in the `combined_sce` colData.")
  }

  # Perform integration with specified method
  if (integration_method == "fastmnn") {
    integrated_sce <- integrate_fastmnn(combined_sce,
                                        batch_column,
                                        ...)
    # Add corrected expression to combined_sce if specified
    if (return_corrected_expression) {
      assay(combined_sce, "fastMNN_corrected") <- assay(integrated_sce, "reconstructed")
    }
    # define PCs
    integrated_pcs <- reducedDim(integrated_sce, "corrected")
  }
  if (integration_method == "harmony") {
    # here the result is the PCs:
    integrated_pcs <- integrate_harmony(combined_sce,
                                        batch_column,
                                        harmony_do_PCA,
                                        harmony_covariate_cols,
                                        ...)
  }


  # Filter assays
  assays_to_remove <- setdiff(assayNames(combined_sce),
                              retain_uncorrected_assays)
  for (assay_name in assays_to_remove) {
    assay(combined_sce, assay_name) <- NULL
  }

  # Add PCs from integration into combined_sce
  reducedDim(integrated_sce, glue::glue("{integration_method}_PCA")) <- integrated_pcs

  # Return the combined_sce with integration results included
  return(combined_sce)
}



#' Integrate a list of SingleCellExperiment (SCE) objects using the fastMNN
#' approach from the `batchelor` package.
#' Source: http://www.bioconductor.org/packages/release/bioc/html/batchelor.html
#'
#' @param combined_sce A combined SCE object as prepared by `scpcaTools::merge_sce_list()`.
#' @param batch_column The variable in `combined_sce` indicating batches.
#' @param ... Additional arguments to pass into `batchelor::fastMNN()`
#'
#' @return The integrated SCE object
#'
#' @import SingleCellExperiment
integrate_fastmnn <- function(combined_sce,
                              batch_column,
                              ...) {

  # Check the combined_sce for logcounts
  if (!("logcounts" %in% names(assays(combined_sce)))) {
    stop("The `combined_sce` object requires a `logcounts` assay for fastMNN integration.")
  }

  # Perform integration
  integrated_sce <- batchelor::fastMNN(combined_sce,
                                       batch = colData(combined_sce)[,batch_column],
                                       ...)

  return(integrated_sce)
}




#' Integrate a list of SingleCellExperiment (SCE) objects using the harmony
#' approach from the `harmony` package.
#' Source: https://cran.r-project.org/web/packages/harmony/index.html
#'
#' @param combined_sce A combined SCE object as prepared by `scpcaTools::merge_sce_list()`.
#' @param batch_column The variable in `combined_sce` indicating batches.
#' @param harmony_do_PCA Whether integration should perform PCA (`TRUE`) and start with logcounts, or
#'  use existing PCs to start (`FALSE`).
#' @param harmony_covariate_cols A vector of other columns in the combined SCE besides
#'  `batch_column` should be considered as covariates during integration.
#' @param ... Additional arguments to pass into `harmony::harmonyMatrix()`
#'
#' @return The integrated SCE object
#'
#' @import SingleCellExperiment
integrate_harmony <- function(combined_sce,
                              batch_column,
                              harmony_do_PCA,
                              harmony_covariate_cols,
                              ...) {

  # Check the combined_sce, depending on do_PCA, and define the harmony input_matrix accordingly
  if (harmony_do_PCA) {
    if (!("logcounts" %in% names(assays(combined_sce)))) {
      stop("The combined_sce object requires a `logcounts` assay for harmony integration when `harmony_do_PCA` is `TRUE`.")
    }
    input_matrix <- logcounts(combined_sce)
  } else {
    if (!("PCA" %in% reducedDimNames(combined_sce))) {
      stop("The combined_sce object requires a `PCA` reduced dimension for harmony integration when `harmony_do_PCA` is `FALSE`.")
    }
    input_matrix <- reducedDim(combined_sce, "PCA")
  }

  # Setup harmony metadata

  # First, ensure that harmony_covariate_cols are present
  if (!(all(harmony_covariate_cols %in% names(colData(combined_sce))))) {
    stop("The combined_sce object does not contain all columns specified in `harmony_covariate_cols`.")
  }
  # Second, we need to remove the column `cell_id` if it already exists in the SCE
  #  since it will conflict with harmony.
  if ("cell_id" %in% names(colData(combined_sce))) {
    colData(combined_sce)$cell_id <- NULL
  }

  # Now we can prepare the metadata such that `cell_id` is the colData rownames per
  #  harmony's expectation
  covariate_cols <- unique(c(batch_column, harmony_covariate_cols))
  harmony_metadata <- tibble::as_tibble(colData(combined_sce),
                                        rownames = "cell_id") %>%
    dplyr::select(cell_id,
                  batch_column,
                  dplyr::all_of(covariate_cols))



  # Perform integration
  harmony_results <- harmony::HarmonyMatrix(input_matrix,
                                            meta_data = harmony_metadata,
                                            do_pca    = harmony_do_PCA,
                                            vars_use  = covariate_cols,
                                            ...)
  # return resulting PCs
  return(harmony_results)
}
