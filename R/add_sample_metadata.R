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
    stop("`sample_id_column` not found in `metadata_df`")
  }

  # filter to relevant sample ids
  metadata_df <- metadata_df |>
    dplyr::filter(.data$sample_id %in% metadata(sce)$sample_id)

  # check that sample id in object corresponds to the sample id column in metadata being added
  if(!all.equal(sort(metadata_df$sample_id), sort(metadata(sce)$sample_id))){
    stop("Sample ids in SCE object do not match the contents of `sample_id_column` in `metadata_df`")
  }

  # add sample metadata
  metadata(sce)$sample_metadata <- metadata_df

  return(sce)

}
