#' Read in counts data processed with Alevin or Alevin-fry
#'
#' @param quant_dir Path to directory where output files are located.
#' @param intron_mode Logical indicating if the files included alignment to intronic regions.
#'   Default is FALSE.
#' @param usa_mode Logical indicating if Alevin-fry was used, if the USA mode was invoked.
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

  # check for cmd_info
  cmd_info_file <- file.path(quant_dir, "cmd_info.json")
  if(!file.exists(cmd_info_file)){
    stop("Missing cmd_info.json in Alevin output directory")
  } else {
    cmd_info <- jsonlite::read_json(cmd_info_file)
  }
  salmon_version <- cmd_info$salmon_version

  if(usa_mode) {
    # read in counts using read_usa mode
    counts <- read_usa_mode(quant_dir)
  } else {
    counts <- read_tximport(quant_dir)
  }

  if (intron_mode | usa_mode) {
    counts <- collapse_intron_counts(counts, which_counts)
  }
  sce <- SingleCellExperiment(list(counts = counts))

  # set metadata fields
  metadata(sce)$pipeline <- "alevin"
  metadata(sce)$salmon_version <- cmd_info$salmon_version
  return(sce)
}

#' Read in counts data processed with Alevin-fry in USA mode
#'
#' @param quant_dir Path to directory where output files are located.
#'
#' @return unfiltered and uncollapsed gene x cell counts matrix
#'
read_usa_mode <- function(quant_dir){

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

  quant_json_path <- file.path(quant_dir, "quant.json")
  if(!file.exists(quant_json_path)){
    # file for alevin-fry < 0.4.1
    quant_json_path <- file.path(quant_dir, "meta_info.json")
    if(!file.exists(quant_json_path)){
      stop("Missing quant.json (or meta_info.json) in Alevin output directory")
    }
  }

  # check that USA mode is true in JSON file
  quant_json <- jsonlite::fromJSON(quant_json_path)
  if(quant_json$usa_mode != "TRUE"){
    stop("Output files not in USA mode")
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
