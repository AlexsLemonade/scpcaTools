
#' Convert SingleCellExperiment objects to AnnData file stored as HDF5 file
#'
#' @param sce SingleCellExperiment object to be converted to AnnData as an HDF5 file
#' @param anndata_file Path to output AnnData file. Must be an `.h5` or `.hdf5`
#'
#' @return SingleCellExperiment object equivalent to the AnnData object found in the
#'   HDF5 file produced.
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
    stop("The zellkonverter package must be installed to convert objects to AnnData. No output file written.")
  }

  if(!is(sce,"SingleCellExperiment")){
    stop("Input must be a SingleCellExperiment object.")
  }

  # check that filename is in the proper format for writing h5
  if(!(stringr::str_ends(anndata_file, ".hdf5|.h5"))){
    stop("`anndata_file` must end in either '.hdf5' or '.h5'")
  }

  # remove miQC model from metadata
  if(!is.null(metadata(sce)$miQC_model)){
    metadata(sce)$miQC_model <- NA
    warning("miQC model cannot be converted between SCE and AnnData.")
  }

  # export SCE object as AnnData to HDF5 file
  zellkonverter::writeH5AD(sce, file = anndata_file)

  # read HDF5 file back in to create SCE
  converted_sce <- zellkonverter::readH5AD(anndata_file)

  invisible(converted_sce)

}
