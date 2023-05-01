#' Add genetic demultiplexing results to a SCE object from Vireo
#'
#' @param sce SingleCellExperiment object.
#' @param vireo_table A data frame of vireo results.
#'   The expected format is that of the vireo `donor_ids.tsv` file, minimally containing columns
#'   'cell' and 'donor_id' (but likely to contain other statistics as well).
#'
#'
#' @return SingleCellExperiment with a `colData` column `vireo_sampleid` containing the confident demultiplexing calls.
#'   Other results imported from the vireo data frame are included in `colData` as well,
#'   with the prefix `vireo_`. See https://vireosnp.readthedocs.io/ for more information.
#'
#' @import SingleCellExperiment
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # add cell calls from a vireo data frame
#'   add_demux_vireo(sce, vireo_table = vireo_df)
#' }
add_demux_vireo <- function(sce, vireo_table){
  # check that input is a SingleCellExperiment
  if(!is(sce, "SingleCellExperiment")){
    stop("`sce` must be a SingleCellExperiment object")
  }

  # check required columns of vireo_table
  vireo_expected_cols = c("cell", "donor_id")
  if(!(is(vireo_table, "data.frame") &
       all(vireo_expected_cols %in% colnames(vireo_table)))){
    stop("`vireo_table` must be a data frame with columns `cell` and `donor_id`")
  }

  # rename results with prefix & normalize
  vireo_table <- vireo_table |>
    dplyr::rename_with(~ paste0("vireo_", .x)) |>
    # add sampleid with finalized single results
    dplyr::mutate(vireo_sampleid = ifelse(.data$vireo_donor_id %in% c("doublet", "unassigned"),
                                          NA_character_,
                                          .data$vireo_donor_id)) |>
    dplyr::relocate("vireo_sampleid") # move vireo_sampleid first

  if(!all(vireo_table$vireo_cell %in% colnames(sce))){
    warning("Cell id(s) from `vireo_table` do not match cells in `sce` object.")
  }
  # if sample_id is present in the SCE metadata, check that samples are as expected
  if(!is.null(metadata(sce)$sample_id)){
    vireo_samples <- unique(vireo_table$vireo_sampleid[!is.na(vireo_table$vireo_sampleid)])
    if(!all(vireo_samples %in% metadata(sce)$sample_id)){
      warning("Sample IDs in `vireo_table` do not match those in the `sce` metadata.")
    }
  }

  # merge vireo results with SCE column names
  vireo_merge <- data.frame(cell = colnames(sce)) |>
    dplyr::left_join(vireo_table, by = c("cell" = "vireo_cell")) |>
    tibble::column_to_rownames("cell")

  # add vireo results to SCE colData
  colData(sce)[, colnames(vireo_merge)] <- vireo_merge

  return(sce)
}
