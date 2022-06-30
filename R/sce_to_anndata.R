
sce_to_anndata <- function(sce){

  if(!is(sce,"SingleCellExperiment")){
    stop("Input must be a SingleCellExperiment object.")
  }

  anndata <- basilisk::basiliskRun(fun = function(sce){

    adata <- zellkonverter::SCE2AnnData(sce)

    zellkonverter::AnnData2SCE(adata)

  }, env = zellkonverterAnnDataEnv, sce = sce)

  return(anndata)
}
