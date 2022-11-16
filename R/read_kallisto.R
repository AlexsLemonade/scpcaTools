#' Read in counts data processed with Kallisto
#'
#' @param quant_dir Path to directory where output files are located.
#' @param intron_mode Logical indicating if the files included alignment to intronic regions.
#'   Default is FALSE.
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
#'               intron_mode = TRUE,
#'               which_counts = "unspliced")
#' }
read_kallisto <- function(quant_dir,
                          intron_mode = FALSE,
                          round_counts = TRUE) {

  if(!is.logical(intron_mode)){
    stop("intron_mode must be set as TRUE or FALSE")
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

  # generate the SCE object containing either counts and spliced assays or just counts assay
  sce <- build_sce(counts,
                   intron_mode,
                   round_counts)

  return(sce)
}
