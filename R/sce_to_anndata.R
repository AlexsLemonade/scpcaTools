#' Convert SingleCellExperiment objects to AnnData file stored as HDF5 file
#'
#' @param sce SingleCellExperiment object to be converted to AnnData as an HDF5 file
#' @param anndata_file Path to output AnnData file. Must end in `.h5`, `.hdf5`, or `.h5ad`
#' @param x_assay_name Name of assay in SCE object to save as X in AnnData. Default is "counts".
#' @param compression Type of compression to use when exporting hdf5 files. Options are
#'   "none", "gzip", or "lzf". Default is "gzip".
#' @param ... Any arguments to be passed into zellkonverter::writeH5AD
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
sce_to_anndata <- function(
    sce,
    anndata_file,
    x_assay_name = "counts",
    compression = c("gzip", "none", "lzf"),
    ...) {
  if (!requireNamespace("zellkonverter", quietly = TRUE)) {
    stop("The zellkonverter package must be installed to convert objects to AnnData. No output file written.")
  }

  if (!is(sce, "SingleCellExperiment")) {
    stop("Input must be a SingleCellExperiment object.")
  }

  if (ncol(sce) < 2) {
    stop("Input SingleCellExperiment must contain at least 2 cells.")
  }

  if (nrow(sce) < 2) {
    stop("Input SingleCellExperiment must contain at least 2 genes or features.")
  }

  # check that filename is in the proper format for writing h5
  if (!(stringr::str_ends(anndata_file, ".hdf5|.h5|.h5ad"))) {
    stop("`anndata_file` must end in either '.hdf5', '.h5', '.h5ad'")
  }

  # make sure assay is found in sce object
  if (!x_assay_name %in% assayNames(sce)) {
    stop("`x_assay_name` is not an assay in `sce`")
  }

  compression <- match.arg(compression)

  # assign SCE to new variable to avoid modifying input SCE
  sce_to_convert <- sce

  # keep only vectors and data frames; no objects or DataFrames
  metadata_to_keep <- metadata(sce_to_convert) |>
    purrr::keep(\(x) {
      is.atomic(x) || is.data.frame(x)
    })


  # prepare data frames for anndata conversion:
  # - remove any list columns
  # - ensure any NA values in character columns are formatted as character
  # - ensure any
  metadata_to_keep <- metadata_to_keep |>
    purrr::map(
      \(meta) {
        if (is.data.frame(meta)) {
          meta <- meta |>
            dplyr::select(dplyr::where(\(col) !is.list(col))) |>
            dplyr::mutate(
              dplyr::across(dplyr::where(is.character), \(col) tidyr::replace_na(col, "NA")),
              dplyr::across(dplyr::where(\(col) all(is.na(col))), as.logical)
            ) |>
            # ensure it's not a tibble; zellkonverter drops these since no python type
            as.data.frame()
        } else {
          meta
        }
      }
    )


  # print out warning that removed objects won't be converted
  removed_metadata <- setdiff(names(metadata(sce_to_convert)), names(metadata_to_keep))
  if (length(removed_metadata) > 0) {
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
    X_name = x_assay_name,
    compression = compression,
    ...
  )
  invisible(sce)
}
