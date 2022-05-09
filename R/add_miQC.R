#' Add miQC model and probability compromised to a SingleCellExperiment object
#'
#'
#' @param sce SingleCellExperiment object.
#' @param posterior_cutoff Optional posterior cutoff used for the miQC filtering calculation (default 0.75)
#' @param seed An optional random seed for reproducibility.
#'
#' @return SingleCellExperiment with prob_compromised column added to colData
#'   If the model fits properly, there will also be a miQC_model slot added to the metadata.
#'
#' @import SummarizedExperiment
#'
#' @export
#'
add_miQC <- function(sce, posterior_cutoff = 0.75, seed = NULL){
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
  # check if miQC_pass exists
  if(!is.null(sce$miQC_pass)){
    warning("miQC_pass was already calculated and will be replaced.")
  }

  # set seed
  set.seed(seed)

  # generate linear mixture model of probability of cells being compromised
  model <- NULL
  pass_cells <- NULL
  model_attempt <- 0

  # This can fail in a few ways, so we will wrap the next steps in a while/try loop
  while(model_attempt < 3 &&
        (!is(model, "flexmix") || length(model@components) < 2 ) &&
        is.null(pass_cells)){
    model_attempt <- model_attempt + 1
    try({
      model <- miQC::mixtureModel(sce)
      # filter step can fail, so we will try this too, keeping the passing cell ids
      pass_cells <- colnames(
        miQC::filterCells(sce,
                          model = model,
                          posterior_cutoff = posterior_cutoff,
                          verbose = FALSE)
      )
    }, silent = TRUE)
  }

  if (!is(model, "flexmix") || length(model@components) < 2){
    # no good fit, fill prob_compromised with NA
    sce$prob_compromised <- NA_real_
  }else{
    # add the prob_compromised value for all cells
    sce <- miQC::filterCells(sce,
                             model = model,
                             posterior_cutoff = 1,
                             enforce_left_cutoff = FALSE,
                             verbose = FALSE)
    # test whether cells passed filtering
    if(is.null(pass_cells)){
      sce$miQC_pass <- NA
    }else{
      sce$miQC_pass <- colnames(sce) %in% pass_cells
    }
    metadata(sce)$miQC_model <- model
  }
  return(sce)
}
