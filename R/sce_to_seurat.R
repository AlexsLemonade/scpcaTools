
#' Convert SingleCellExperiment object to Seurat object
#'
#' @param sce SingleCellExperiment object
#' @param assay_name The assay name (default "counts") to include in
#'   the Seurat object. This name will be applied as the assay name in
#'   the Seurat object. If the default "counts" assay is used, then the
#'   assay name will instead be "RNA," consistent with Seurat defaults.
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
sce_to_seurat <- function(sce,
                          assay_name = "counts") {

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The Seurat package must be installed to create a Seurat object. No output returned.")
  }

  if(!is(sce,"SingleCellExperiment")){
    stop("Input must be a SingleCellExperiment object.")
  }

  if (!assay_name %in% assayNames(sce)) {
    stop("Provided assay name is not present in the SingleCellExperiment object.")
  }

  # remove miQC model from metadata
  if(!is.null(metadata(sce)$miQC_model)){
    metadata(sce)$miQC_model <- NULL
    message("miQC model will not be included in Seurat object.")
  }


  # Seurat counts need to be integers
  # extract the given `assay_name`
  sce_counts <- round(assay(sce, assay_name))

  # Seurat will not like zero count calls
  sce_sum <- Matrix::colSums(sce_counts)
  seurat_cells <- names(sce_sum)[sce_sum > 0]
  sce_counts <- sce_counts[, seurat_cells]

  # convert metadata to Seurat compatible formats
  coldata <- as.data.frame(colData(sce))
  rowdata <- as.data.frame(rowData(sce))


  # create seurat object, with new assay name
  # use "RNA" as assay name if `counts` is being used from SCE
  if (assay_name == "counts") {
    assay_name <- "RNA"
  }
  seurat_obj <- Seurat::CreateSeuratObject(counts = sce_counts,
                                           meta.data = coldata,
                                           assay = assay_name)

  # add rowdata and metadata separately adding it while creating leaves the slots empty without warning
  seurat_obj[[assay_name]]@var.features <- rowdata
  seurat_obj@misc <- metadata(sce)

  # grab names of altExp, if any
  alt_names <- altExpNames(sce)

 # add each altExp from the SCE to the Seurat object
  for (name in alt_names){

    # we want to make sure that we are only adding counts for the cells
    # that currently exist in the object
    seurat_obj_cells <- colnames(seurat_obj)

    # round counts and calculate total counts
    alt_counts <- round(counts(altExp(sce, name)))

    # subset altExp counts and seurat object to shared cells
    cells_to_keep <- intersect(seurat_obj_cells, colnames(alt_counts))
    alt_counts <- alt_counts[, cells_to_keep]

    # subset seurat object to only have cells that will be in alt_counts before adding new assay
    seurat_obj <- seurat_obj[, cells_to_keep]

    # add new assay along with associated row data, if any
    seurat_obj[[name]] <- Seurat::CreateAssayObject(counts = alt_counts)
    rowdata <- as.data.frame(rowData(altExp(sce, name)))
    seurat_obj[[name]]@var.features <- rowdata

    # check that altExp data is present as a new assay, since Seurat sometimes fails without warning
    if(is.null(seurat_obj[[name]])){
      warning("Unable to convert altExp data found in SCE to Seurat.")
    }

  }

  return(seurat_obj)

}
