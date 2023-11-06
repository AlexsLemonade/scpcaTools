#' Convert SingleCellExperiment objects to AnnData file stored as HDF5 file
#'
#' @param sce SingleCellExperiment object to be converted to AnnData as an HDF5 file
#' @param anndata_file Path to output AnnData file. Must be an `.h5` or `.hdf5`
#' @param x_assay_name Name of assay in SCE object to save as X in AnnData. Default is "counts".
#'
#' @return original SingleCellExperiment object used as input (invisibly)
#' **Note that any columns present in the `rowData` of an SCE object that contains
#' duplicated information, e.g. duplicate gene identifiers, are converted to
#' categorical data by the `anndata` package.
#'
#' @import SingleCellExperiment
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sce_to_anndata(
#'   sce = sce_object,
#'   anndata_file = "test_anndata.h5"
#' )
#' }
sce_to_anndata <- function(sce, anndata_file, x_assay_name = "counts") {
  if (!requireNamespace("zellkonverter", quietly = TRUE)) {
    stop("The zellkonverter package must be installed to convert objects to AnnData. No output file written.")
  }

  if (!is(sce, "SingleCellExperiment")) {
    stop("Input must be a SingleCellExperiment object.")
  }

  if(ncol(sce) < 2){
    stop("Input SingleCellExperiment must contain at least 2 cells.")
  }

  if(nrow(sce) < 2){
    stop("Input SingleCellExperiment must contain at least 2 genes or features.")
  }

  # check that filename is in the proper format for writing h5
  if (!(stringr::str_ends(anndata_file, ".hdf5|.h5"))) {
    stop("`anndata_file` must end in either '.hdf5' or '.h5'")
  }

  # make sure assay is found in sce object
  if (!x_assay_name %in% assayNames(sce)) {
    stop("`x_assay_name` is not an assay in `sce`")
  }

  # assign SCE to new variable to avoid modifying input SCE
  sce_to_convert <- sce

  # remove any objects or dataframes
  metadata_to_keep <- metadata(sce_to_convert) |>
    purrr::discard(is.object) |>
    purrr::discard(is.list)

  # print out warning that removed objects won't be converted
  removed_metadata <- setdiff(names(metadata(sce_to_convert)),  names(metadata_to_keep))
  if(length(removed_metadata) > 0){
    glue::glue("{removed_metadata} cannot be converted between SCE and AnnData.") |>
      purrr::walk(message)
  }

  # reset metadata
  metadata(sce_to_convert) <- metadata_to_keep

  # zellkonverter (or pandas) does not like NA values to be stored in character vectors
  colData(sce_to_convert) <- colData(sce_to_convert) |>
    as.data.frame() |>
    dplyr::mutate(
      dplyr::across(dplyr::where(\(x) all(is.na(x))), as.logical)
    ) |>
    DataFrame(row.names = colnames(sce_to_convert))

  rowData(sce_to_convert) <- rowData(sce_to_convert) |>
    as.data.frame() |>
    dplyr::mutate(
      dplyr::across(dplyr::where(\(x) all(is.na(x))), as.logical)
    ) |>
    DataFrame(row.names = rownames(sce_to_convert))

  # export SCE object as AnnData to HDF5 file
  zellkonverter::writeH5AD(sce_to_convert,
    file = anndata_file,
    X_name = x_assay_name
  )
  invisible(sce)
}
