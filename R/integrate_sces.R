# This function can perform the crux integration based on the given integration method


#' Integrate a merged set of SingleCellExperiment objects using a specified
#'  integration method. The final SCE object will contain an additional reducedDim
#'  entitled `\{integration_method\}_PCA`. All original SCE assays are retained.
#'  `fastMNN` integration uses all default setting, but additional parameters can be supplied.
#'  `harmony` integration uses existing PCs, but additional parameters can be supplied.
#'
#' @param merged_sce A merged SCE object as prepared by `scpcaTools::merge_sce_list()`.
#' @param integration_method The integration method to use. One of `fastMNN` or
#'  `harmony` (case-sensitive)
#' @param batch_column The column in the merged SCE object indicating batches
#' @param covariate_cols A vector of additional columns in the merged SCE to
#'   consider as covariates during integration. Currently, this is used only by
#'   `harmony`.
#' @param batch_lambda The ridge regression penalty to use when correcting using the `batch_column`.
#'   This is only used by `harmony`. Default is 1.
#' @param covariate_lambda A vector of ridge regression penalties to use with each covariate.
#'   A lambda value must be provided for each covariate column.
#'   This is only used by `harmony`.
#' @param return_corrected_expression A boolean indicating whether corrected expression
#'   values determined by the given integration method should be included in the
#'   integrated SCE object. Note that `harmony` does not calculate corrected expression,
#'   so this argument is ignored for this integration method. Default is `FALSE`.
#' @param seed Random seed to set for integration. This is only set if the value is not `NULL`.
#' @param ... Any additional parameters to be passed to the given integration method
#'
#' @return An updated merged SCE object containing a new reduced dimension named
#'   `\{integration_method\}_PCA` containing the corrected PCs. In the case of `fastMNN`
#'   integration, the SCE also includes a corrected expression assay (`fastMNN_corrected`)
#'
#' @import SingleCellExperiment
#'
#' @export
integrate_sces <- function(merged_sce,
                           integration_method = c("fastMNN", "harmony"),
                           batch_column = "sample",
                           covariate_cols = c(),
                           batch_lambda = 1,
                           covariate_lambda = c(),
                           return_corrected_expression = FALSE,
                           seed = NULL,
                           ...) {
  # Set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # make sure that input is a SingleCellExperiment
  if (!is(merged_sce, "SingleCellExperiment")) {
    stop("The `merged_sce` must be a SingleCellExperiment object created with `scpcaTools::merge_sce_list().`")
  }


  # Check integration_method
  integration_method <- match.arg(integration_method)

  # Check batch_column
  if (!(batch_column %in% names(colData(merged_sce)))) {
    stop("The provided `batch_column` column must be in the `merged_sce` colData.")
  }
  # Ensure there are no NAs
  batches <- colData(merged_sce)[[batch_column]]
  if (any(is.na(batches))) {
    stop("`NA` values are not allowed in the `batch_column`.")
  }
  # Ensure >=2 batches
  if (length(unique(batches)) < 2) {
    stop("At least two batches are required for integration.")
  }


  # Perform integration with specified method
  if (integration_method == "fastMNN") {
    if (length(covariate_cols) > 0) {
      warning("fastMNN cannot use additional covariates, so `covariate_cols` will be ignored.")
    }

    integrated_sce <- integrate_fastmnn(
      merged_sce,
      batch_column,
      ...
    )
    # Add corrected expression to merged_sce if specified
    if (return_corrected_expression) {
      assay(merged_sce, "fastMNN_corrected") <- assay(integrated_sce, "reconstructed")
    }
    # define PCs
    integrated_pcs <- reducedDim(integrated_sce, "corrected")
  } else if (integration_method == "harmony") {
    # warn that this is not possible
    if (return_corrected_expression) {
      warning("`harmony` does not calculate corrected expression values, so none can be returned.")
    }

    # here the result is the PCs:
    integrated_pcs <- integrate_harmony(
      merged_sce,
      batch_column,
      covariate_cols,
      batch_lambda,
      covariate_lambda,
      ...
    )
  }

  # Add PCs from integration into merged_sce (name is lowercase)
  reducedDim(merged_sce, glue::glue("{integration_method}_PCA")) <- integrated_pcs

  # Return the merged_sce with integration results included
  return(merged_sce)
}



#' Integrate a list of SingleCellExperiment (SCE) objects using the fastMNN
#' approach from the `batchelor` package.
#' Source: http://www.bioconductor.org/packages/release/bioc/html/batchelor.html
#'
#' @param merged_sce A merged SCE object as prepared by `scpcaTools::merge_sce_list()`.
#' @param batch_column The variable in `merged_sce` indicating batches.
#' @param ... Additional arguments to pass into `batchelor::fastMNN()`
#'
#' @return An integrated SCE object produced by `fastMNN`
#'
#' @import SingleCellExperiment
integrate_fastmnn <- function(merged_sce,
                              batch_column,
                              ...) {
  # Check that batchelor is installed
  if (!requireNamespace("batchelor", quietly = TRUE)) {
    stop("The `batchelor` package must be installed to use fastMNN.")
  }

  # Check the merged_sce for logcounts
  if (!("logcounts" %in% names(assays(merged_sce)))) {
    stop("The `merged_sce` object requires a `logcounts` assay for fastMNN integration.")
  }

  # Perform integration
  integrated_sce <- batchelor::fastMNN(merged_sce,
    batch = colData(merged_sce)[, batch_column],
    ...
  )

  return(integrated_sce)
}




#' Integrate a list of SingleCellExperiment (SCE) objects using the harmony
#' approach from the `harmony` package.
#' Source: https://cran.r-project.org/web/packages/harmony/index.html
#'
#' @param merged_sce A merged SCE object as prepared by `scpcaTools::merge_sce_list()`.
#' @param batch_column The column in the merged SCE object indicating batches
#'   being integrated.
#' @param covariate_cols A vector of other columns to consider as
#'   covariates during integration.
#' @param batch_lambda The ridge regression penalty to use when correcting using the `batch_column`.
#'   Default is 1.
#' @param covariate_lambda A vector of ridge regression penalties to use with each covariate.
#'   A lambda value must be provided for each covariate column.
#' @param ... Additional arguments to pass into `harmony::RunHarmony()`
#'
#' @return Integrated PCs as calculated by `harmony`
#'
#' @import SingleCellExperiment
integrate_harmony <- function(merged_sce,
                              batch_column,
                              covariate_cols = c(),
                              batch_lambda = 1,
                              covariate_lambda = c(),
                              ...) {
  # Check that harmony is installed
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("The `harmony` package must be installed to use this integration method.")
  }

  # Check the merged_sce
  if (!("PCA" %in% reducedDimNames(merged_sce))) {
    stop("The merged_sce object requires a `PCA` reduced dimension for harmony integration.")
  }

  # Check covariate columns
  if (!(all(covariate_cols %in% names(colData(merged_sce))))) {
    stop("The provided covariate columns are not all present in the `merged_sce` colData.")
  }

  # if no covariate lambda provided, set to 1
  if (length(covariate_lambda) == 0) {
    covariate_lambda <- rep(1, length(covariate_cols))
  }

  # check that the lambda values match up
  if (length(covariate_lambda) != length(covariate_cols)) {
    stop("The number of covariate columns must be equal to the number of covariate lambda values.")
  }

  # Setup harmony metadata
  covariate_cols <- c(batch_column, covariate_cols)

  # get a vector of all lambdas
  lambda_vec <- c(batch_lambda, covariate_lambda)

  # Perform integration
  harmony_object <- harmony::RunHarmony(
    merged_sce,
    group.by.vars = covariate_cols,
    reduction.use = "PCA",
    lambda = lambda_vec,
    verbose = FALSE,
    ...
  )

  # pull out corrected pcs
  harmony_results <- reducedDim(harmony_object, "HARMONY")

  # Ensure PCs have rownames
  if (is.null(rownames(harmony_results))) {
    rownames(harmony_results) <- rownames(reducedDim(merged_sce, "PCA"))
  }

  # return resulting PCs
  return(harmony_results)
}
