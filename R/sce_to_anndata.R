
#' Convert SingleCellExperiment objects to AnnData file stored as HDF5 file
#'
#' @param sce SingleCellExperiment object to be converted to AnnData
#' @param anndata_file Name of
#'
#' @return Converted AnnData object
#'
#' @import SingleCellExperiment
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sce_to_anndata(sce = sce_object,
#'                anndata_file = "test_anndata.h5")
#' }
sce_to_anndata <- function(sce, anndata_file){

  if (!requireNamespace("zellkonverter", quietly = TRUE)) {
    warning("The zellkonverter package must be installed to convert objects to AnnData. Returning unmodified sce object")
    return(sce)
  }

  if(!is(sce,"SingleCellExperiment")){
    stop("Input must be a SingleCellExperiment object.")
  }

  # check that filename is in the proper format for writing h5
  if(!(stringr::str_ends(anndata_file, ".hdf5|.h5"))){
    stop("anndata file must end in either `.hdf5` or .`h5`")
  }

  # remove miQC model from metadata
  if(!is.null(metadata(sce)$miQC_model)){
    metadata(sce)$miQC_model <- NA
    warning("miQC model cannot be converted between SCE and AnnData.")
  }

  zellkonverter::writeH5AD(sce, file = anndata_file)

}
