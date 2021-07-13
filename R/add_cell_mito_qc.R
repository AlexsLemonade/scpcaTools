#' Calculate QC metrics for all reads and mitochondrial subset for each cell in a SingleCellExperiment
#'
#' @param sce SingleCellExperiment object.
#' @param mito Character vector of mitochondrial gene names in the same format as rownames of SingleCellExperiment object.
#' @param ... Any additional arguments to be passed to scater::addPerCellQC.
#'
#' @return SingleCellExperiment with colData slot containing calculated QC metrics.
#'  Columns include sum, detected, subsets_mito_sum, subsets_mito_detected, subsets_mito_percent, and total.
#' @export
#'
#' @examples
#' \dontrun{
#' # add per cell QC metrics using only genes detected > 5% of cells
#' add_cell_mito_qc(sce = sce,
#'                  mito = mito_genes,
#'                  threshold = 5)
#' }
add_cell_mito_qc <- function(sce, mito, ...){

  # check that input is a SingleCellExperiment
  if(!is(sce, "SingleCellExperiment")){
    stop("sce must be a SingleCellExperiment object")
  }

  # check that mito is not empty, otherwise resulting colData will be innacurate
  if(length(mito) == 0){
    warning("Mitochondrial gene list not used, resulting mitochondrial metrics will not be accurate.")
  }

  # check that mito genes are present in sce
  if(length(intersect(mito, rownames(sce))) == 0){
    warning("sce does not contain genes corresponding to the list of mito gene names.")
  }

  # add per cell QC with mitochondrial subset
  scater::addPerCellQC(sce,
                       subsets = list(mito = mito[mito %in% rownames(sce)]),
                       ...)
}
