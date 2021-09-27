#' Calculate QC metrics for all reads and mitochondrial subset for each cell in a SingleCellExperiment
#'
#' @param sce SingleCellExperiment object.
#' @param mito Character vector of mitochondrial gene names in the same format as rownames of SingleCellExperiment object.
#' @param miQC Logical indicating whether or not to calculate the posterior probability of a cell being compromised using
#'   the linear mixture model in miQC. Default is FALSE.
#' @param ... Any additional arguments to be passed to scuttle::addPerCellQCMetrics.
#'
#' @return SingleCellExperiment with colData slot containing calculated QC metrics.
#'  Columns include sum, detected, mito_sum, mito_detected, mito_percent, and total.
#'
#' @import SummarizedExperiment
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # add per cell QC metrics using only genes detected > 5% of cells
#' add_cell_mito_qc(sce = sce,
#'                  mito = mito_genes,
#'                  threshold = 5)
#' }
add_cell_mito_qc <- function(sce, mito, miQC = FALSE, ...){

  # check that input is a SingleCellExperiment
  if(!is(sce, "SingleCellExperiment")){
    stop("sce must be a SingleCellExperiment object")
  }

  # check that mito is not empty, otherwise resulting colData will be innacurate
  if(length(mito) == 0){
    stop("Mitochondrial gene list not used, cannot calculate mitochondrial metrics.")
  }

  # check that mito genes are present in sce, otherwise colData will have 0's for mito columns
  if(length(intersect(mito, rownames(sce))) == 0){
    warning("sce does not contain genes corresponding to the list of mito gene names")
  }

  # check that miQC is logical
  if(!is.logical(miQC)){
    stop("miQC must be set as TRUE or FALSE")
  }

  # add per cell QC with mitochondrial subset
  sce <- scuttle::addPerCellQCMetrics(
    sce,
    subsets = list(mito = mito[mito %in% rownames(sce)]),
    ...
  )


  if(miQC){
    # generate linear mixture model of probability of cells being compromised
    model <- miQC::mixtureModel(sce)

    # use filter cells, but keeping all cells, to add a column to colData with prob_compromised
    sce <- miQC::filterCells(sce, model, posterior_cutoff = 1, verbose = FALSE)
    metadata(sce)$miQC_model <- model
  }

  return(sce)
}
