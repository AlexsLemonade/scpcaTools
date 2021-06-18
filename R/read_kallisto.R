#' Read in counts data processed with Kallisto
#'
#' @param quant_dir Full path to directory where output files are located.
#' @param intron_mode Boolean indicating if the files included alignment to intronic regions. Default is FALSE.
#' @param which_counts The type of counts (cDNA or intron) to include if alignment to intronic regions is TRUE. Default is FALSE.
#' @param intron_metadata_path Full path to a two column tsv file containing gene names for both spliced and intronic regions.
#'           Only required if intron_mode = TRUE and usa_mode = FALSE.
#'
#' @return SingleCellExperiment of unfiltered gene x cell counts matrix
#' @export
#'
#' @examples
#' \dontrun{
#' read_kallisto(quant_dir,
#' intron_mode = TRUE,
#' which_counts = "intron",
#' intron_metadata_path = "<path to intron metadata.tsv>")
#' }
read_kallisto <- function(quant_dir, intron_mode = FALSE, which_counts = c("cDNA", "intron"), intron_metadata_path) {

  base_file <- file.path(quant_dir, "counts", "gene_count")
  counts <- Matrix::readMM(paste0(base_file,".mtx"))%>%
    t() %>% # transpose to gene x cell orientation
    as("dgCMatrix") # compress sparse matrix
  dimnames(counts) <- list(readLines(paste0(base,".", "genes.txt")),
                           readLines(paste0(base,".barcodes.txt")))

  if(intron_mode == TRUE) {
    sce <- collapse_intron_counts(counts, which_counts, intron_metadata_path)
  } else {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  }

  return(sce)
}
