#' Add data frame with sample metadata to an SCE object
#'
#' @param sce SingleCellExperiment object to add metadata too
#' @param metadata_df A dataframe with metadata to add to the sce.
#' @param sample_id_column Name of column in `metadata_df` that contains the sample ids.
#'   The contents of this column should be equivalent to the `sample_id` in `metadata(sce)`
#'
#' @return SingleCellExperiment object with a `metadata_df` added to the `sample_metadata`
#'   slot of the metadata list
#' @export
#'
add_sample_metadata <- function(sce,
                                metadata_df,
                                sample_id_column){

  # make sure input is sce
  if (!is(sce, "SingleCellExperiment")) {
    stop("`sce` must be a SingleCellExperiment object.")
  }

  # make sure metadata is a dataframe
  if (!is.data.frame(metadata_df)){
    stop("`metadata_df` must be a data.frame")
  }

  # sample id column should be present in metadata
  if(!sample_id_column %in% colnames(metadata_df)){
    stop("`sample_id_column` not found in `metadata_df`")
  }

  # check that sample id in object corresponds to the sample id column in metadata being added
  sce_sample_ids <- metadata(sce)$sample_id |>
    stringr::str_split_1(pattern = ",")
  metadata_sample_ids <- metadata_df[[sample_id_column]]

  if(!all.equal(sort(sce_sample_ids), sort(metadata_sample_ids))){
    stop("Sample ids in SCE object do not match the contents of `sample_id_column` in `metadata_df`")
  }

  # add sample metadata
  metadata(sce)$sample_metadata <- metadata_df

  return(sce)

}
