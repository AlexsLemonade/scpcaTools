#' Add sample metadata to colData of SCE object
#'
#' @param sce SingleCellExperiment object to add sample metadata to colData.
#'   Must contain the `batch_column` as a column in `colData`.
#' @param join_columns A character value giving the names of the columns in colData
#'  to use for joining with the `sample_metadata`. Default is `library_id`
#'
#' @return SingleCellExperiment object with the contents of `metadata(sce)$sample_metadata`
#'   added as new columns to the `colData(sce)`.
#' @import SingleCellExperiment
#'
#' @export

metadata_to_coldata <- function(sce,
                                join_columns = "library_id") {

  # make sure input is sce
  if (!is(sce, "SingleCellExperiment")) {
    stop("`sce` must be a SingleCellExperiment object.")
  }

  # check that sample_metadata present in metadata of sce object
  # if it's not present, just skip everything and return the original object
  metadata_list <- metadata(sce)
  if(!"sample_metadata" %in% names(metadata_list)){
    warning("No `sample_metadata` in SCE object to add.")
    return(sce)
  }

    # check that batch id column is present
    if(!all(join_columns %in% colnames(colData(sce)))) {
      stop("One or more of the specified `join_columns` are not in `colData(sce)`.")
    }

    # pull out sample metadata
    sample_metadata_df <- metadata_list$sample_metadata
    # check that batch column exists
    if(!all(join_columns %in% colnames(sample_metadata_df))){
      stop("One or more of the specified `join_columns` are not a column in `metadata(sce)$sample_metadata.")
    }

    # check that join columns match in colData and sample metadata
    mismatching_columns <- purrr::map(
      join_columns,
      \(column){
        if(!all(sce[[column]] %in% sample_metadata_df[[column]])){
          return(column)
        }
      }) |> unlist()

    if(!is.null(mismatching_columns)){
      warning(
        glue::glue("Contents of {mismatching_columns} do not match the `metadata(sce)$sample_metadata.`")
        )
    }

    # join coldata with sample metadata
    coldata_df <- as.data.frame(colData(sce)) |>
      dplyr::left_join(sample_metadata_df, by = join_columns)

    # replace existing coldata
    colData(sce) <- DataFrame(coldata_df,
                              row.names = rownames(coldata_df))

    }

  # return modified sce with sample metadata
  return(sce)

}
