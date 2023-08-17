#' Add data frame with sample metadata to an SCE object
#'
#' @param sce SingleCellExperiment object to add metadata too
#' @param metadata_df A dataframe with metadata to add to the sce. Must contain a
#'   column named `sample_id`.
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

  # sample id column should be present in metadata
  if(!"sample_id" %in% colnames(metadata_df)){
    stop("No column named `sample_id` in `metadata_df`")
  }

  # check that sample ids in the object are found in the metadata data frame
  if(!all(metadata(sce)$sample_id %in% metadata_df$sample_id)){
    stop("Sample ids in SCE object are not all present in `metadata_df`")
  }

  # filter to relevant sample ids
  metadata_df <- metadata_df |>
    dplyr::filter(.data$sample_id %in% metadata(sce)$sample_id)

  # add sample metadata
  metadata(sce)$sample_metadata <- metadata_df

  return(sce)

}
