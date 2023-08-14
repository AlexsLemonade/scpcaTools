#' Add data frame with sample metadata to an SCE object
#'
#' @param sce SingleCellExperiment object to add metadata too
#' @param metadata_df A dataframe with metadata to add to the sce.
#'
#' @return SingleCellExperiment object with a `metadata_df` added to the `sample_metadata`
#'   slot of the metadata list
#' @export
#'
add_sample_metadata <- function(sce,
                                metadata_df){

  # make sure input is sce
  if (!is(sce, "SingleCellExperiment")) {
    stop("`sce` must be a SingleCellExperiment object.")
  }

  # make sure metadata is a dataframe
  if (!is.data.frame(metadata_df)){
    stop("`metadata_df` must be a data.frame")
  }

  # add sample metadata
  metadata(sce)$sample_metadata <- metadata_df

  return(sce)

}
