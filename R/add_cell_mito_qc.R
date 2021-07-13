#' Title
#'
#' @param sce
#' @param mito
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
add_cell_mito_qc <- function(sce, mito, threshold){

  if(!is(sce, "SingleCellExperiment")){
    stop("sce must be a SingleCellExperiment object")
  }

  if(is.na(mito)){
    stop("Must input list of mitochondrial genes using mito argument.")
  }
  if(is.na(intersect(mito, rownames(sce)))){
    stop("sce does not contain genes corresponding to the list of mito gene names.")
  }


  scater::addPerCellQC(sce,
                       subsets = list(mito = mito[mito %in% rownames(sce)]),
                       threshold = threshold)
}
