#' Add sample metadata to colData of SCE object
#'
#' @param sce SingleCellExperiment object to add sample metadata to colData.
#'   Must contain the `batch_column` as a column in `colData`.
#' @param batch_column A character value giving the colData column name
#'  to differentiate originating SingleCellExperiment objects. Often these values
#'  are unique library IDs. Default value is `library_id`.
#'
#' @return SingleCellExperiment object with the contents of `metadata(sce)$sample_metadata`
#'   added as new columns to the `colData(sce)`.
#' @export
#'
metadata_to_coldata <- function(sce,
                                batch_column = "library_id") {

  # make sure input is sce
  if (!is(sce, "SingleCellExperiment")) {
    stop("`sce` must be a SingleCellExperiment object.")
  }

  # check that sample_metadata present in metadata of sce object
  metadata_list <- metadata(sce)
  if(!"sample_metadata" %in% names(metadata_list)){
    stop("No `sample_metadata` in SCE object to add.")
  }

  # check that batch id column is present
  if(!batch_column %in% colnames(colData(sce))) {
    stop(glue::glue("{batch_column} is not in `colData(sce)`."))
  }

  # pull out sample metadata
  sample_metadata_df <- metadata_list$sample_metadata
  # check that batch column exists
  if(!batch_column %in% colnames(sample_metadata_df)){
    stop(glue::glue("{batch_column} is not a column in `metadata(sce)$sample_metadata"))
  }

  # check that batches match in colData and sample metadata
  coldata_batches <- sce[[batch_column]] |>
    unique()
  metadata_batches <- sample_metadata_df[[batch_column]] |>
    unique()
  if(!all(metadata_batches %in% coldata_batches)){
    stop("The batches in `metadata(sce)$sample_metadata` are not found in the sce object.")
  }

  # join coldata with sample metadata
  coldata_df <- as.data.frame(colData(sce)) |>
    dplyr::left_join(sample_metadata_df, by = batch_column)

  # replace existing coldata
  colData(sce) <- DataFrame(coldata_df,
                            row.names = rownames(coldata_df))

  # return modified sce with sample metadata
  return(sce)

}
