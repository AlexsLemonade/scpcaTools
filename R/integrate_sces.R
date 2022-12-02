# This function can perform the crux integration based on the given integration method


#' Integrate a combined set of SingleCellExperiment objects using a specified
#'  integration method. The final SCE object will contain an additional reducedDim
#'  entitled `{integration_method}_PCA`. All original SCE assays are retained.
#'  `fastMNN` integration uses all default setting, but additional parameters can be supplied.
#'  `harmony` integration uses existing PCs, but additional parameters can be supplied.
#'
#' @param combined_sce A combined SCE object as prepared by `scpcaTools::merge_sce_list()`.
#' @param integration_method The integration method to use. One of `fastMNN` or
#'  `harmony` (case-sensitive)
#' @param covariate_cols A vector of covariates (e.g. the batches) to consider during integration.
#'   For both `fastMNN` and `harmony`, this should be SCE column indicating batches.
#'   For `harmony` specifically, additionally SCE columns may be supplied.
#'   Default is `c("sample")`.
#' @param return_corrected_expression A boolean indicating whether corrected expression
#'   values determined by the given integration method should be included in the
#'   integrated SCE object. Note that `harmony` does not calculate corrected expression,
#'   so this argument is ignored for this integration method. Default is `TRUE`.
#' @param seed Random seed to set for integration. This is only set if the value is not `NULL`.
#' @param ... Any additional parameters to be passed to the given integration method
#'
#' @return An update combined SCE object containing a new reduced dimension named
#'   `{integration_method}_PCA` containing the corrected PCs. In the case of `fastMNN`
#'   integration, the SCE also includes a corrected expression assay (`fastMNN_corrected`)
#'
#' @import SingleCellExperiment
#'
#' @export
integrate_sces <- function(combined_sce,
                           integration_method = c("fastMNN", "harmony"),
                           covariate_cols = "sample",
                           return_corrected_expression = TRUE,
                           seed = NULL,
                           ...) {

  # Set seed
  set.seed(seed)

  # make sure that input is a SingleCellExperiment
  if(!is(combined_sce, "SingleCellExperiment")){
    stop("The `combined_sce` must be a SingleCellExperiment object created with `scpcaTools::merge_sce_list().`")
  }

  # Check integration_method
  integration_method <- match.arg(integration_method)

  # Check covariate_cols
  if (length(covariate_cols) < 1) {
    stop("Error: You must specify covariate (e.g. batch) columns in `covariate_cols`.")
  }
  if (length(covariate_cols) > 1 & integration_method == "fastMNN") {
    stop("Error: fastMNN can only handle a single covariate (batch).")
  }
  if (!(all(covariate_cols %in% names(colData(combined_sce))))) {
    stop("The provided covariate columns are not all present in the `combined_sce` colData.")
  }

  # Perform integration with specified method
  if (integration_method == "fastMNN") {

    integrated_sce <- integrate_fastmnn(combined_sce,
                                        covariate_cols,
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
                                        covariate_cols,
                                        ...)
  }

  # Add PCs from integration into combined_sce (name is lowercase)
  reducedDim(combined_sce, glue::glue("{integration_method}_PCA")) <- integrated_pcs

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
#' @param covariate_cols A vector of other columns, including the batch, to
#'   consider as covariates during integration.
#' @param ... Additional arguments to pass into `harmony::harmonyMatrix()`
#'
#' @return The integrated SCE object
#'
#' @import SingleCellExperiment
integrate_harmony <- function(combined_sce,
                              covariate_cols,
                              ...) {

  # Check the combined_sce
  if (!("PCA" %in% reducedDimNames(combined_sce))) {
    stop("The combined_sce object requires a `PCA` reduced dimension for harmony integration.")
  }

  # Setup harmony metadata
  harmony_metadata <- tibble::as_tibble(colData(combined_sce)) %>%
    dplyr::select(dplyr::all_of(covariate_cols))

  # Perform integration
  harmony_results <- harmony::HarmonyMatrix(reducedDim(combined_sce, "PCA"),
                                            meta_data = harmony_metadata,
                                            vars_use  = covariate_cols,
                                            ...)
  # Ensure PCs have rownames
  if (is.null(rownames(harmony_results))) {
    rownames(harmony_results) <- rownames(reducedDim(combined_sce, "PCA"))
  }


  # return resulting PCs
  return(harmony_results)
}
