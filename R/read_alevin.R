#' Read in counts data processed with Alevin or Alevin-fry
#'
#' @param quant_dir Path to directory where output files are located.
#' @param intron_mode Logical indicating if the files included alignment to intronic regions.
#'   Default is FALSE.
#' @param usa_mode Logical indicating if Alevin-fry was used, if USA mode was invoked.
#'   Default is FALSE.
#' @param which_counts Which type of counts should be included,
#'   only counts aligned to spliced cDNA ("spliced") or all spliced and unspliced cDNA ("unspliced").
#'   Applies if `intron_mode` or `usa_mode` is TRUE.
#'   Default is "spliced".
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
                        intron_mode = FALSE,
                        usa_mode = FALSE,
                        which_counts = c("spliced", "unspliced")){

  which_counts <- match.arg(which_counts)

  # checks for intron_mode and usa_mode
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
  meta <- read_alevin_metadata(quant_dir, which_counts)

  # Read the count data
  if(usa_mode) {
    if(meta$usa_mode != TRUE){
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
#' @param quant_dir Path to directory where output files are located.
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
  counts <- Matrix::readMM(file = file.path(quant_dir, "alevin", "quants_mat.mtx"))%>%
    Matrix::t() %>%
    as("dgCMatrix")
  cols <- readLines(file.path(quant_dir, "alevin", "quants_mat_cols.txt"))
  rows <- readLines(file.path(quant_dir, "alevin", "quants_mat_rows.txt"))
  dimnames(counts) <- list(cols, rows)
  return(counts)
}

#' Read in counts data processed with Alevin or alevin-fry in tximport-compatible formats
#'
#' @param quant_dir Path to directory where output files are located.
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
#' @param quant_dir Path alevin output directory
#' @param transcript_type A string indicating output types
#'
#' @return A list containing alevin run metadata,
#'   with NULL values for missing elements.
#'
#' @noRd
read_alevin_metadata <- function(quant_dir){
  cmd_info_path <- file.path(quant_dir, "cmd_info.json")
  permit_json_path <- file.path(quant_dir, "generate_permit_list.json")
  collate_json_path <- file.path(quant_dir, "collate.json")
  quant_json_path <- file.path(quant_dir, "quant.json")
  if(!file.exists(quant_json_path)){
    # file for alevin-fry < 0.4.1
    quant_json_path <- file.path(quant_dir, "meta_info.json")
  }

  # get cmd_info, which should always be present
  if (file.exists(cmd_info_path)){
    cmd_info <- jsonlite::read_json(cmd_info_path)
  } else {
    stop("cmd_info.json is missing")
  }

  # Read other info files if they exist. Otherwise, create dummy values
  if (file.exists(permit_json_path)){
    permit_info <- jsonlite::read_json(permit_json_path)
  } else {
    permit_info <- list()
  }
  if (file.exists(quant_json_path)){
    collate_info <- jsonlite::read_json(collate_json_path)
  } else {
    collate_info <- list()
  }
  if (file.exists(quant_json_path)){
    quant_info <- jsonlite::read_json(quant_json_path)
  } else {
    quant_info <- list()
  }

  # Create a metadata list
  meta <- list(salmon_version = cmd_info[['salmon_version']],
               reference_index = cmd_info[['index']])

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
  meta$usa_mode <- quant_info[['usa_mode']]

  return(meta)
}

