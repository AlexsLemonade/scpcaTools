#' Read in counts data processed with Kallisto
#'
#' @param quant_dir Full path to directory where output files are located.
#' @param intron_mode Boolean indicating if the files included alignment to intronic regions. Default is FALSE.
#' @param which_counts If intron_mode is TRUE, which type of counts should be included,
#'        only counts aligned to spliced cDNA ("spliced") or all spliced and unspliced cDNA ("unspliced").
#'        Default is "spliced".
#'
#' @return SingleCellExperiment of unfiltered gene x cell counts matrix
#' @export
#'
#' @examples
#' \dontrun{
#' read_kallisto(quant_dir,
#' intron_mode = TRUE,
#' which_counts = "intron")
#' }
read_kallisto <- function(quant_dir, intron_mode = FALSE, which_counts = c("spliced", "unspliced")) {

  which_counts <- match.arg(which_counts)

  if(!is.boolean(intron_mode)){
    stop("intron_mode must be set as TRUE or FALSE")
  }

  kallisto_files <- list("gene_count.mtx", "gene_count.genes.txt", "gene_count.barcodes.txt")
  kallisto_dir <- file.path(quant_dir, "counts")

  if(!dir.exists(file.path(kallisto_dir))){
    stop("Missing kallisto directory with output files")
  }
  for (file in kallisto_files){
    if(!file.exists(file.path(kallisto_dir, file))){
      error_message <- paste("Missing kallisto output file", file, sep = " ")
      stop(error_message)
    }
  }

  base_file <- file.path(kallisto_dir, "gene_count")
  counts <- Matrix::readMM(paste0(base_file,".mtx"))%>%
    t() %>% # transpose to gene x cell orientation
    as("dgCMatrix") # compress sparse matrix
  dimnames(counts) <- list(readLines(paste0(base,".genes.txt")),
                           readLines(paste0(base,".barcodes.txt")))

  if(intron_mode == TRUE) {
    collapsed_counts <- collapse_intron_counts(counts, which_counts)
    sce <- SingleCellExperiment(list(counts = collapsed_counts))
  } else {
    sce <- SingleCellExperiment(list(counts = counts))
  }

  return(sce)
}
