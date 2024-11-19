# nolint start: cyclocomp_linter.
#' Read in counts data processed with Alevin or Alevin-fry
#'
#' @param quant_dir Path to directory where output files are located.
#' @param fry_mode Logical indicating if Alevin-fry was used for quantification.
#'   Implies the input data is in matrix market format.
#'   Default is FALSE.
#' @param include_unspliced Whether or not to include the unspliced reads in the counts matrix.
#'   If TRUE, the main "counts" assay will contain unspliced reads and spliced reads and an additional "spliced"
#'   assay will contain spliced reads only. If TRUE, requires that data has been aligned to a reference containing
#'   spliced and unspliced reads.
#'   Default is TRUE.
#' @param feature_data Logical indicating if the data being read in contains feature data.
#'   Default is FALSE. If feature data is being read in there are no spliced or unspliced reads.
#'   All reads will be counted and stored in the `counts` assay.
#' @param round_counts Logical indicating in the count matrix should be rounded to integers on import.
#'   Only used if `fry_mode` is FALSE. Default is TRUE.
#' @param library_id Optional library identifier
#' @param sample_id Optional sample identifier.
#'   If multiplexed samples are included in a library, this may be a vector.
#' @param project_id Optional project identifier.
#' @param tech_version Technology or kit used to process library (i.e. 10Xv3, 10Xv3.1).
#' @param assay_ontology_term_id Optional Experimental Factor Ontology term associated with the provided `tech_version`
#' @param seq_unit Optional sequencing unit associated with the library (i.e. cell, nucleus)
#'
#' @return SingleCellExperiment of unfiltered gene x cell counts matrix.
#' @export
#'
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#'
#' @examples
#' \dontrun{
#'
#' # Import output files processed with either Alevin or Alevin-fry with alignment to
#' # cDNA only, including only spliced cDNA in final counts matrix
#' read_alevin(quant_dir)
#'
#' # Import output files processed with either Alevin or Alevin-fry with alignment to
#' # cDNA + introns and including all unspliced cDNA in final counts matrix
#' read_alevin(quant_dir,
#'   include_unspliced = TRUE
#' )
#'
#' # Import output files processed with alevin-fry USA mode
#' # including all unspliced cDNA in final counts matrix
#' read_alevin(quant_dir,
#'   fry_mode = TRUE,
#'   include_unspliced = TRUE
#' )
#' }
read_alevin <- function(
    quant_dir,
    fry_mode = FALSE,
    include_unspliced = TRUE,
    feature_data = FALSE,
    round_counts = TRUE,
    library_id = NULL,
    sample_id = NULL,
    project_id = NULL,
    tech_version = NULL,
    assay_ontology_term_id = NULL,
    seq_unit = NULL) {
  # checks for *_mode
  if (!is.logical(fry_mode)) {
    stop("fry_mode must be set as TRUE or FALSE")
  }
  if (!is.logical(include_unspliced)) {
    stop("include_unspliced must be set as TRUE or FALSE")
  }
  if (!is.logical(round_counts)) {
    stop("round_counts must be set as TRUE or FALSE")
  }

  # make sure that include unspliced and feature data are not both set as TRUE
  if (include_unspliced && feature_data) {
    stop("Feature data does not have unspliced reads, cannot use `include_unspliced=TRUE`")
  }

  # check that the expected quant directory exists
  if (!dir.exists(file.path(quant_dir, "alevin"))) {
    stop("Missing alevin directory with output files")
  }

  # read metadata
  meta <- read_alevin_metadata(
    quant_dir,
    tech_version,
    assay_ontology_term_id = assay_ontology_term_id,
    seq_unit = seq_unit,
    library_id = library_id,
    sample_id = sample_id,
    project_id = project_id
  )

  # if alevin-fry USA and MTX format directly create SCE object with fishpond
  if (fry_mode) {
    # only check for usa mode for non-feature data
    if (!feature_data && !meta$usa_mode) {
      stop("Output files not in USA mode")
    }

    # define assays to include in SCE object based on include_unspliced
    if (include_unspliced) {
      assay_formats <- list("counts" = c("S", "A", "U"), "spliced" = c("S", "A"))
      meta$transcript_type <- c("total", "spliced")
    } else if (!feature_data) {
      assay_formats <- list("counts" = c("S", "A"))
      meta$transcript_type <- "spliced"
    } else {
      # if using feature data, use default output format for loadFry of scRNA
      assay_formats <- "scrna"
      meta$transcript_type <- "feature_counts"
    }

    # must be both alevin-fry and usa mode to use fishpond
    sce <- fishpond::loadFry(
      fryDir = quant_dir,
      outputFormat = assay_formats
    )
    if (round_counts) {
      assays(sce) <- lapply(assays(sce), round)
    }
  } else {
    # read in any non-USA formatted alevin-fry data or Alevin data
    counts <- read_tximport(quant_dir)

    # set transcript type based on including unspliced or not
    if (include_unspliced) {
      meta$transcript_type <- c("total", "spliced")
    } else if (!feature_data) {
      meta$transcript_type <- "spliced"
    } else {
      meta$transcript_type <- "feature_counts"
    }

    # generate the SCE object containing either counts and spliced assays or just counts assay
    sce <- build_sce(
      counts,
      include_unspliced,
      round_counts
    )
  }

  # remove any empty droplets that may have snuck in with rounding
  sce <- sce[, which(colSums(counts(sce)) > 0)]

  # add the metadata to the SCE
  meta$include_unspliced <- include_unspliced
  metadata(sce) <- meta

  return(sce)
}
# nolint end

#' Read in counts data processed with Alevin or alevin-fry in tximport-compatible formats
#'
#' @param quant_dir Path to alevin output directory.
#'
#' @return unfiltered & uncollapsed gene x cell counts matrix
#'
read_tximport <- function(quant_dir) {
  # check that all files exist in quant_dir
  # use tximport for all non-usa mode
  alevin_files <- c("quants_mat_cols.txt", "quants_mat_rows.txt", "quants_mat.gz")

  missing <- !file.exists(file.path(quant_dir, "alevin", alevin_files))
  if (any(missing)) {
    missing_files <- paste(alevin_files[missing], collapse = ", ")
    stop(paste0("Missing Alevin output file(s): ", missing_files))
  }

  txi <- suppressMessages(tximport::tximport(
    file.path(quant_dir, "alevin", "quants_mat.gz"),
    type = "alevin"
  ))
  counts <- as(txi$counts, "CsparseMatrix")
}



#' Read alevin metadata from json files
#'
#' @param quant_dir Path alevin output directory.
#' @param tech_version Technology or kit used to process library (i.e. 10Xv3, 10Xv3.1).
#' @param assay_ontology_term_id Optional Experimental Factor Ontology term associated with the provided `tech_version`
#' @param seq_unit Optional sequencing unit associated with the library (i.e. cell, nucleus)
#' @param library_id Optional library identifier
#' @param sample_id Optional sample identifier.
#'   If multiplexed samples are included in a library, this may be a vector.
#' @param project_id Optional project identifier.
#'
#' @return A list containing alevin run metadata,
#'   with NULL values for missing elements.
#'
#' @noRd
read_alevin_metadata <- function(quant_dir,
                                 tech_version,
                                 assay_ontology_term_id = NULL,
                                 seq_unit = NULL,
                                 library_id = NULL,
                                 sample_id = NULL,
                                 project_id = NULL) {
  cmd_info_path <- file.path(quant_dir, "cmd_info.json")
  permit_json_path <- file.path(quant_dir, "generate_permit_list.json")
  # Unused file, but leaving for future reference
  # collate_json_path <- file.path(quant_dir, "collate.json")
  quant_json_path <- file.path(quant_dir, "quant.json")
  aux_meta_path <- file.path(quant_dir, "aux_info", "meta_info.json")

  if (!file.exists(quant_json_path)) {
    # file for alevin-fry < 0.4.1
    quant_json_path <- file.path(quant_dir, "meta_info.json")
  }

  # get cmd_info and aux_info/meta_info.json, which should always be present
  if (file.exists(cmd_info_path)) {
    cmd_info <- jsonlite::read_json(cmd_info_path)
  } else {
    stop("cmd_info.json is missing")
  }
  if (file.exists(aux_meta_path)) {
    aux_meta <- jsonlite::read_json(aux_meta_path)
  } else {
    stop("meta_info.json in aux_info folder is missing")
  }

  # Read other info files if they exist. Otherwise, create dummy values
  if (file.exists(permit_json_path)) {
    permit_info <- jsonlite::read_json(permit_json_path)
  } else {
    permit_info <- list()
  }
  if (file.exists(quant_json_path)) {
    quant_info <- jsonlite::read_json(quant_json_path)
  } else {
    quant_info <- list()
  }

  # Create a metadata list
  meta <- list(
    library_id = library_id,
    sample_id = sample_id,
    project_id = project_id,
    salmon_version = cmd_info$salmon_version,
    reference_index = cmd_info[["index"]],
    total_reads = aux_meta[["num_processed"]],
    mapped_reads = aux_meta[["num_mapped"]]
  )
  # using $ notation  for `salmon_version` to get partial matching due to salmon 1.5.2 bug
  # see https://github.com/COMBINE-lab/salmon/issues/691

  # if we have permit_info data, we used alevin-fry, otherwise alevin
  if (length(permit_info) == 0) {
    meta$mapping_tool <- "alevin"
  } else {
    meta$mapping_tool <- "alevin-fry"
  }

  # Add other metadata
  # assume all alevin-fry tool versions are the same
  meta$alevinfry_version <- permit_info[["version_str"]]
  meta$af_permit_type <- permit_info[["permit-list-type"]]
  meta$af_resolution <- quant_info[["resolution_strategy"]]
  meta$af_tx2gene <- cmd_info[["tgMap"]]
  meta$usa_mode <- quant_info[["usa_mode"]]
  meta$af_num_cells <- quant_info[["num_quantified_cells"]]
  meta$tech_version <- tech_version
  meta$assay_ontology_term_id <- assay_ontology_term_id
  meta$seq_unit <- seq_unit


  return(meta)
}
