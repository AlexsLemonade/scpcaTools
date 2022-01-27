#' Read in counts data processed with Alevin or Alevin-fry
#'
#' @param quant_dir Path to directory where output files are located.
#' @param mtx_format Logical indicating if input data is in matrix market format.
#'   Default is FALSE.
#' @param intron_mode Logical indicating if the files included alignment to intronic regions.
#'   Default is FALSE.
#' @param usa_mode Logical indicating if Alevin-fry was used, if USA mode was invoked.
#'   Implies the input data is in matrix market format.
#'   Default is FALSE.
#' @param which_counts Which type of counts should be included,
#'   only counts aligned to spliced cDNA ("spliced") or all spliced and unspliced cDNA ("unspliced").
#'   Applies if `intron_mode` or `usa_mode` is TRUE.
#'   Default is "spliced".
#' @param tech_version Technology or kit used to process library (i.e. 10Xv3, 10Xv3.1).
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
#'             intron_mode = TRUE,
#'             which_counts = "unspliced")
#'
#' # Import output files processed with alevin-fry USA mode
#' # including all unspliced cDNA in final counts matrix
#' read_alevin(quant_dir,
#'             usa_mode = TRUE,
#'             which_counts = "unspliced")
#'
#'}
read_alevin <- function(quant_dir,
                        mtx_format = FALSE,
                        intron_mode = FALSE,
                        usa_mode = FALSE,
                        which_counts = c("spliced", "unspliced"),
                        tech_version = NULL){

  which_counts <- match.arg(which_counts)

  # checks for *_mode
  if(!is.logical(mtx_format)){
    stop("mtx_format must be set as TRUE or FALSE")
  }
  if(!is.logical(intron_mode)){
    stop("intron_mode must be set as TRUE or FALSE")
  }
  if(!is.logical(usa_mode)){
    stop("usa_mode must be set as TRUE or FALSE")
  }

  # check that usa_mode and intron_mode are not used together
  if(usa_mode & intron_mode){
    stop("Can only read counts using either usa mode or intron mode.")
  }

  # check that the expected quant directory exists
  if(!dir.exists(file.path(quant_dir, "alevin"))){
    stop("Missing alevin directory with output files")
  }

  # read metadata
  meta <- read_alevin_metadata(quant_dir, tech_version)

  # Read the count data
  if(mtx_format | usa_mode) {
    if(usa_mode & meta$usa_mode != TRUE){
      stop("Output files not in USA mode")
    }
    counts <- read_alevin_mtx(quant_dir)
  } else {
    counts <- read_tximport(quant_dir)
  }
  if (intron_mode | usa_mode) {
    counts <- collapse_intron_counts(counts, which_counts)
    meta$transcript_type <- which_counts
  }

  # make the SCE object
  sce <- SingleCellExperiment(assays = list(counts = counts),
                              metadata = meta)
  return(sce)
}

#' Read in counts data processed with Alevin-fry in with mtx output.
#'
#' @param quant_dir Path to alevin output directory.
#'
#' @return unfiltered and uncollapsed gene x cell counts matrix
#'
read_alevin_mtx <- function(quant_dir){

  # check that all files exist in quant_dir
  alevin_files <- c("quants_mat_cols.txt", "quants_mat_rows.txt", "quants_mat.mtx")

  # check that all files exist in quant directory
  if(!dir.exists(file.path(quant_dir, "alevin"))){
    stop("Missing alevin directory with output files")
  }

  missing <- !file.exists(file.path(quant_dir, "alevin", alevin_files))
  if(any(missing)) {
    missing_files <- paste(alevin_files[missing], collapse = ", ")
    stop(paste0("Missing Alevin output file(s): ", missing_files))
  }

  # read in .mtx files
  counts <- Matrix::readMM(file = file.path(quant_dir, "alevin", "quants_mat.mtx"))|>
    Matrix::t() |>
    as("dgCMatrix")
  cols <- readLines(file.path(quant_dir, "alevin", "quants_mat_cols.txt"))
  rows <- readLines(file.path(quant_dir, "alevin", "quants_mat_rows.txt"))
  dimnames(counts) <- list(cols, rows)
  return(counts)
}

#' Read in counts data processed with Alevin or alevin-fry in tximport-compatible formats
#'
#' @param quant_dir Path to alevin output directory.
#'
#' @return unfiltered & uncollapsed gene x cell counts matrix
#'
read_tximport <- function(quant_dir){

  # check that all files exist in quant_dir
  # use tximport for all non-usa mode
  alevin_files <- c("quants_mat_cols.txt", "quants_mat_rows.txt", "quants_mat.gz")

  missing <- !file.exists(file.path(quant_dir, "alevin", alevin_files))
  if(any(missing)) {
    missing_files <- paste(alevin_files[missing], collapse = ", ")
    stop(paste0("Missing Alevin output file(s): ", missing_files))
  }

  txi <- suppressMessages(tximport::tximport(
    file.path(quant_dir, "alevin", "quants_mat.gz"),
    type = "alevin"
  ))
  counts <- txi$counts
}

#' Read alevin metadata from json files
#'
#' @param quant_dir Path alevin output directory.
#' @param tech_version Technology or kit used to process library (i.e. 10Xv3, 10Xv3.1).
#'
#' @return A list containing alevin run metadata,
#'   with NULL values for missing elements.
#'
#' @noRd
read_alevin_metadata <- function(quant_dir, tech_version){
  cmd_info_path <- file.path(quant_dir, "cmd_info.json")
  permit_json_path <- file.path(quant_dir, "generate_permit_list.json")
  # Unused file, but leaving for future reference
  # collate_json_path <- file.path(quant_dir, "collate.json")
  quant_json_path <- file.path(quant_dir, "quant.json")
  aux_meta_path <- file.path(quant_dir, "aux_info", "meta_info.json")

  if(!file.exists(quant_json_path)){
    # file for alevin-fry < 0.4.1
    quant_json_path <- file.path(quant_dir, "meta_info.json")
  }

  # get cmd_info and aux_info/meta_info.json, which should always be present
  if (file.exists(cmd_info_path)){
    cmd_info <- jsonlite::read_json(cmd_info_path)
  } else {
    stop("cmd_info.json is missing")
  }
  if (file.exists(aux_meta_path)){
    aux_meta <- jsonlite::read_json(aux_meta_path)
  } else {
    stop("meta_info.json in aux_info folder is missing")
  }

  # Read other info files if they exist. Otherwise, create dummy values
  if (file.exists(permit_json_path)){
    permit_info <- jsonlite::read_json(permit_json_path)
  } else {
    permit_info <- list()
  }
  if (file.exists(quant_json_path)){
    quant_info <- jsonlite::read_json(quant_json_path)
  } else {
    quant_info <- list()
  }

  # Create a metadata list
  meta <- list(salmon_version = cmd_info$salmon_version,
               reference_index = cmd_info[['index']],
               total_reads = aux_meta[['num_processed']],
               mapped_reads = aux_meta[['num_mapped']])
  # using $ notation  for `salmon_version` to get partial matching due to salmon 1.5.2 bug
  # see https://github.com/COMBINE-lab/salmon/issues/691

  # if we have permit_info data, we used alevin-fry, otherwise alevin
  if (length(permit_info) == 0){
    meta$mapping_tool <- "alevin"
  } else {
    meta$mapping_tool <- "alevin-fry"
  }

  # Add other metadata
  # assume all alevin-fry tool versions are the same
  meta$alevinfry_version <- permit_info[['version_str']]
  meta$af_permit_type <- permit_info[['permit-list-type']]
  meta$af_resolution <- quant_info[['resolution_strategy']]
  meta$af_tx2gene <- cmd_info[['tgMap']]
  meta$usa_mode <- quant_info[['usa_mode']]
  meta$af_num_cells <- quant_info[['num_quantified_cells']]
  meta$tech_version <- tech_version


  return(meta)
}

