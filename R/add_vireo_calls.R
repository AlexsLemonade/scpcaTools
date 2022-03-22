#' Add genetic demultiplexing results to a SCE object from Vireo
#'
#' @param sce SingleCellExperiment object.
#' @param vireo_df A data frame of vireo results
#'
#'
#' @return SingleCellExperiment with a `colData` column `hashedDrops_sampleid` containing the confident demultiplexing calls.
#'   Other results from `DropletUtils::hashedDrops()` are included in the `rowData` for the cellhash `altExp`,
#'   with the prefix `hashedDrops_`. See [DropletUtils::hashedDrops()] for the contents of these fields.
#'
#' @import SingleCellExperiment
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # add cell calls from DropletUtils::hashedDrops()
#'   add_demux_hashedDrops(sce = sce)
#' }
add_demux_vireo <- function(sce, vireo_df){
  # check that input is a SingleCellExperiment
  if(!is(sce, "SingleCellExperiment")){
    stop("sce must be a SingleCellExperiment object")
  }

  # check required columns of vireo_df
  vireo_expected_cols = c("cell", "donor_id")
  if(!all(vireo_expected_cols %in% colnames(vireo_df))){
    stop("Columns in vireo_df do not include 'cell' and/or 'donor_id', as required.")
  }

  # rename results with prefix & normalize
  vireo_df <- vireo_df |>
    dplyr::rename_with(~ paste0("vireo_", .x)) |>
    # add sampleid with finalized single results
    dplyr::mutate(vireo_sampleid = ifelse(vireo_donor_id %in% c("doublet","unassigned"),
                                          NA_character_,
                                          vireo_donor_id)) |>
    dplyr::relocate(vireo_sampleid) # move vireo_sampleid first

  if(!all(vireo_df$vireo_cell %in% colnames(sce))){
    warning("Cell id(s) from vireo do not match cells in sce object.")
  }

  # merge vireo results with SCE column names
  vireo_merge <- data.frame(cell = colnames(sce)) |>
    dplyr::left_join(vireo_df, by = c("cell" = "vireo_cell")) |>
    tibble::column_to_rownames("cell")

  # add vireo results to SCE colData
  colData(sce)[, colnames(vireo_merge)] <- vireo_merge

  return(sce)
}
