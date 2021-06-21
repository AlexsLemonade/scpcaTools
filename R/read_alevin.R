#' Read in counts data processed with Alevin or Alevin-fry
#'
#' @param quant_dir Full path to directory where output files are located.
#' @param intron_mode Boolean indicating if the files included alignment to intronic regions.
#'        Default is FALSE.
#' @param usa_mode Boolean indicating if Alevin-fry was used, if the USA mode was invoked.
#'        Default is FALSE.
#' @param which_counts If intron_mode is TRUE, which type of counts should be included,
#'        only counts aligned to spliced cDNA ("spliced") or all spliced and unspliced cDNA ("unspliced").
#'        Default is "spliced".
#'
#' @return SingleCellExperiment of unfiltered gene x cell counts matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' read_alevin(quant_dir,
#' intron_mode = TRUE,
#' usa_mode = TRUE,
#' which_counts = "unspliced")
#' }
read_alevin <- function(quant_dir, intron_mode = FALSE, usa_mode = FALSE,
                        which_counts = c("spliced", "unspliced")){

  which_counts <- match.arg(which_counts)

  # checks for intron_mode and usa_mode
  if(!is.boolean(intron_mode)){
    stop("intron_mode must be set as TRUE or FALSE")
  }
  if(!is.boolean(usa_mode)){
    stop("usa_mode must be set as TRUE or FALSE")
  }

  ## check that intron_metadata is provided
  # if(intron_mode == TRUE & usa_mode == FALSE){
  #   if(!intron_metadata){
  #     stop("Missing intron metadata table")
  #   } else {
  #     if(colnames(intron_metadata != c("spliced", "intron"))) {
  #       stop("Incorrect column names for intron metadata")
  #     }
  #   }
  # }

  if(usa_mode == TRUE) {
    # read in counts using read_usa mode
    counts <- read_usa_mode(quant_dir, which_counts)
    sce <- SingleCellExperiment(list(counts = counts))

  } else {
    # use tximport for all non-usa mode
    alevin_files <- list("quants_mat_cols.txt", "quants_mat_rows.txt", "quants_mat.gz", "alevin.log")

    # check that all files exist in quant directory
    if(!dir.exists(file.path(quant_dir, "alevin"))){
      stop("Missing alevin directory with output files")
    }
    for (file in alevin_files){
      if(!file.exists(file.path(quant_dir, "alevin", file))){
        error_message <- paste("Missing alevin output file", file, sep = " ")
        stop(error_message)
      }
    }

    if(!file.exists(file.path(quant_dir, "cmd_info.json"))){
      stop("Missing cmd_info.json in Alevin output directory")
    }

    txi <- tximport::tximport(file.path(quant_dir, "alevin", "quants_mat.gz"), type = "alevin")
    sce <- SingleCellExperiment(list(counts = txi$counts))

    # collapse intron counts for intron_mode = TRUE
    if (intron_mode == TRUE) {
      counts < - collapse_intron_counts(counts(sce), which_counts)
      sce <- SingleCellExperiment(list(counts = counts))
    }
  }
  return(sce)
}

#' Read in counts data processed with Alevin-fry in USA mode
#'
#' @param quant_dir Full path to directory where output files are located.
#' @param which_counts If intron_mode is TRUE, which type of counts should be included,
#'        only counts aligned to spliced cDNA ("spliced") or all spliced and unspliced cDNA ("unspliced").
#'        Default is "spliced".
#'
#' @return unfiltered gene x cell counts matrix
#'
#' @examples
read_usa_mode <- function(quant_dir, which_counts = c("spliced", "unspliced")){

  which_counts <- match.arg(which_counts)

  # check that all files exist in quant_dir
  alevin_files <- list("quants_mat_cols.txt", "quants_mat_rows.txt", "quants_mat.mtx", "alevin.log")

  # check that all files exist in quant directory
  if(!dir.exists(file.path(quant_dir, "alevin"))){
    stop("Missing alevin directory with output files")
  }
  for (file in alevin_files){
    if(!file.exists(file.path(quant_dir, "alevin", file))){
      error_message <- paste("Missing alevin output file", file, sep = " ")
      stop(error_message)
    }
  }
  meta_json_path <- file.path(quant_dir, "meta_info.json")
  if(!file.exists(meta_json_path)){
    stop("Missing meta_info.json in Alevin output directory")
  }

  # check that USA mode is true in JSON file
  meta_json <- rjson::fromJSON(file = meta_json_path)
  if(meta_json$usa_mode != "TRUE"){
    stop("Output files not in USA mode")
  }

  # read in .mtx files
  mtx <- Matrix::readMM(file = file.path(quant_dir, "alevin", "quants_mat.mtx"))%>%
    Matrix::t() %>%
    as("dgCMatrix")
  cols <- readr::read_csv(file = file.path(quant_dir, "alevin", "quants_mat_cols.txt"),
                          col_names = "gene_id")
  rows <- readr::read_csv(file = file.path(quant_dir, "alevin", "quants_mat_rows.txt"),
                          col_names = "barcodes")
  dimnames(mtx) <- list(cols$gene_id, rows$barcodes)

  if (which_counts == "spliced") {
    # only combine counts from S + A
    intron_genes <- str_subset(rownames(mtx), "-U$")
    splice_mtx <- mtx[!(rownames(mtx) %in% intron_genes),]
    # remove A from rowname
    rownames(splice_mtx) <- str_remove(rownames(splice_mtx), "-A$")
    # combine counts based on gene name
    counts <- Matrix.utils::aggregate.Matrix(splice_mtx, rownames(splice_mtx))

  } else if (which_counts == "unspliced") {
    # combine counts from U, S, and A
    # remove A from rowname
    rownames(mtx) <- str_remove(rownames(mtx), "-A$")
    rownames(mtx) <- str_remove(rownames(mtx), "-U$")
    # combine counts based on gene name
    counts <- Matrix.utils::aggregate.Matrix(mtx, rownames(mtx))

  }
  return(counts)
}

