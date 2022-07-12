
#' Convert SingleCellExperiment object to Seurat object
#'
#' @param sce SingleCellExperiment object
#'
#' @return Seurat object
#'
#' @import SingleCellExperiment
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sce_to_seurat(sce = sce_object)
#' }
sce_to_seurat <- function(sce){

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The Seurat package must be installed to create a Seurat object. No output returned.")
  }

  if(!is(sce,"SingleCellExperiment")){
    stop("Input must be a SingleCellExperiment object.")
  }

  # Seurat counts need to be integers
  sce_counts <- round(counts(sce))
  # Seurat will not like zero count cells
  sce_sum <- Matrix::colSums(sce_counts)
  seurat_cells <- names(sce_sum)[sce_sum > 0]
  sce_counts <- sce_counts[, seurat_cells]

  # convert metadata to Seurat compatible formats
  coldata <- as.data.frame(colData(sce))
  rowdata <- as.data.frame(rowData(sce))


  # create seurat object
  seurat_obj <- Seurat::CreateSeuratObject(counts = sce_counts,
                                           meta.data = coldata,
                                           var.features = rowdata)

  # add metadata separately adding it while creating leaves the misc slot empty without warning
  seurat_obj@misc = metadata(sce)

  # grab names of altExp, if any
  alt_names <- altExpNames(sce)

 # add each altExp from the SCE to the Seurat object
  for (name in alt_names){

    # we want to make sure that we are only adding counts for the cells
    # that currently exist in the object
    seurat_obj_cells <- colnames(seurat_obj)

    # round counts and calculate total counts
    alt_counts <- round(counts(altExp(sce, name)))
    alt_sum <- Matrix::colSums(alt_counts)

    # identify cells with > 0 counts and intersect with cells already in seurat
    alt_seurat_cells <- names(alt_sum)[alt_sum > 0]
    cells_to_keep <- intersect(seurat_obj_cells, seurat_obj_cells)
    alt_counts <- alt_counts[, cells_to_keep]

    # subset seurat object to only have cells that will be in alt_counts before adding new assay
    seurat_obj <- seurat_obj[, cells_to_keep]

    # add new assay along with associated row data, if any
    rowdata <- as.data.frame(rowData(altExp(sce, name)))
    seurat_obj[[name]] <- Seurat::CreateAssayObject(counts = alt_counts)
    seurat_obj[[name]]@var.features <- rowdata

    # check that altExp data is present as a new assay, since Seurat sometimes fails without warning
    if(!is.null(seurat_obj[[name]])){
      warning("Unable to convert altExp data found in SCE to Seurat.")
    }

  }

  return(seurat_obj)

}
