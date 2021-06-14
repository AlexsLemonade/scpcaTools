#' Title
#'
#' @param quant_dir
#' @param intron_mode
#' @param which_counts
#' @param intron_metadata
#'
#' @return
#' @export
#'
#' @examples
read_kallisto <- function(quant_dir, intron_mode = TRUE, which_counts = c("cDNA", "intron"), intron_metadata) {

  base_file <- file.path(quant_dir, "counts", "gene_count")
  counts <- Matrix::readMM(paste0(base_file,".mtx"))%>%
    t() %>% # transpose to gene x cell orientation
    as("dgCMatrix") # compress sparse matrix
  dimnames(counts) <- list(readLines(paste0(base,".", "genes.txt")),
                           readLines(paste0(base,".barcodes.txt")))

  if(intron_mode == TRUE) {
    sce <- collapse_intron_counts(counts, which_counts, intron_metadata)
  } else {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  }

  return(sce)
}
