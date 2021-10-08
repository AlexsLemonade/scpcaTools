#' Add miQC model and probability compromised to a SingleCellExperiment object
#'
#'
#' @param sce SingleCellExperiment object.
#'
#' @return SingleCellExperiment with prob_compromised column added to colData
#'   If the model fits properly, there will also be a miQC_model slot added to the metadata.
#'
#' @import SummarizedExperiment
#'
#' @export
#'
add_miQC <- function(sce){
  # check that input is a SingleCellExperiment
  if(!is(sce, "SingleCellExperiment")){
    stop("sce must be a SingleCellExperiment object")
  }
  # check that sce has subsets_mito_percent
  if(!("subsets_mito_percent" %in% colnames(colData(sce)))){
    stop("sce must have subsets_mito_percent in the column data. Use scuttle::addPerCellQCMetrics or similar to add it.")
  }
  # check if prob_compromised exists
  if(!is.null(sce$prob_compromised)){
    warning("prob_compromised was already calculated and will be replaced.")
  }

  # generate linear mixture model of probability of cells being compromised
  model <- tryCatch(
    miQC::mixtureModel(sce),
    error = function(x){NA}
  )
  if (class(model) == "flexmix"  && length(model@components) < 2){
    # first model fit gave only one component, give it one more chance.
    model <- miQC::mixtureModel(sce)
  }
  if (class(model) != "flexmix" || length(model@components) < 2){
    # no good fit, fill prob_compromised with NA
    sce$prob_compromised <- NA_real_
  }else{
    sce <- miQC::filterCells(sce, model = model, posterior_cutoff = 1, verbose = FALSE)
    metadata(sce)$miQC_model <- model
  }
  return(sce)
}
