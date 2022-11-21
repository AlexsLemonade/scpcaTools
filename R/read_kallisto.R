#' Read in counts data processed with Kallisto
#'
#' @param quant_dir Path to directory where output files are located.
#' @param include_unspliced Whether or not to include the unspliced reads in the counts matrix.
#'   If TRUE, the main "counts" assay will contain unspliced reads and spliced reads and an additional "spliced"
#'   assay will contain spliced reads only. If TRUE, requires that data has been aligned to a reference contianing
#'   spliced and unspliced reads.
#'   Default is TRUE.
#' @param round_counts Logical indicating in the count matrix should be rounded to integers on import.
#'   Default is TRUE.
#'
#' @return SingleCellExperiment of unfiltered gene x cell counts matrix
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # import output files processed with kallisto with alignment to cDNA + introns,
#' # including all unspliced cDNA counts in final counts matrix
#' read_kallisto(quant_dir,
#'               include_unspliced = TRUE,
#'               which_counts = "unspliced")
#' }
read_kallisto <- function(quant_dir,
                          include_unspliced = TRUE,
                          round_counts = TRUE) {

  if(!is.logical(include_unspliced)){
    stop("include_unspliced must be set as TRUE or FALSE")
  }

  if(!is.logical(round_counts)){
    stop("round_counts must be set as TRUE or FALSE")
  }

  kallisto_files <- c("gene_count.mtx", "gene_count.genes.txt", "gene_count.barcodes.txt")
  kallisto_dir <- file.path(quant_dir, "counts")

  if(!dir.exists(file.path(kallisto_dir))){
    stop("Missing kallisto directory with output files")
  }

  missing <- !file.exists(file.path(kallisto_dir, kallisto_files))
  if(any(missing)) {
    missing_files <- paste(kallisto_files[missing], collapse = ", ")
    stop(paste0("Missing Kallisto output file(s): ", missing_files))
  }

  counts <- Matrix::readMM(file.path(kallisto_dir, "gene_count.mtx"))|>
    BiocGenerics::t() |> # transpose to gene x cell orientation
    as("CsparseMatrix") # compress sparse matrix
  dimnames(counts) <- list(readLines(file.path(kallisto_dir, "gene_count.genes.txt")),
                           readLines(file.path(kallisto_dir, "gene_count.barcodes.txt")))

  # get kallisto version
  run_info_path <- file.path(quant_dir, "run_info.json")
  if(file.exists(run_info_path)){
    run_info <- jsonlite::read_json(run_info_path)
  } else {
    stop("run_info.json is missing")
  }

  # initiate metadata
  meta <- list(
    mapping_tool = "kallisto",
    kallisto_version = run_info[["kallisto_version"]],
    include_unspliced = include_unspliced
  )

  # generate the SCE object containing either counts and spliced assays or just counts assay
  sce <- build_sce(counts,
                   include_unspliced,
                   round_counts)

  # set transcript type based on including unspliced or not
  if(include_unspliced){
    meta$transcript_type <- c("total", "spliced")

  } else {
    meta$transcript_type <- "spliced"
  }

  # add metadata to sce
  metadata(sce) <- meta

  return(sce)
}
