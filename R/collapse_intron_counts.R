#' Merge counts from intron reads with corresponding cDNA reads
#'
#' @param counts Counts matrix with rownames corresponding to gene names and colnames corresponding to cell barcodes.
#'   Can be a counts matrix with intron counts specified by -I or an alevin-fry "USA" matrix,
#'   with intron counts marked by "-U" and ambiguous counts "-A".
#' @param which_counts Which type of counts should be included:
#'   Only counts aligned to spliced cDNA ("spliced") or all spliced and unspliced cDNA ("unspliced" or "total").
#'   Ambiguous counts in USA mode are always included.
#'   Default is "spliced".
#'
#' @return unfiltered gene x cell counts matrix
#'
#' @examples
#' \dontrun{
#'
#' # only keep counts from spliced cDNA in final counts matrix
#' collapse_intron_counts(counts,
#'                        which_counts = "spliced")
#'
#' # include unspliced cDNA in final counts matrix
#' collapse_intron_counts(counts,
#'                        which_counts = "total")
#' }
#'
#'
#'
collapse_intron_counts <- function(counts,
                                   which_counts = c("spliced", "unspliced", "total")){

  which_counts <- match.arg(which_counts)

  introns <- str_detect(rownames(counts), "-I$")
  usa_introns <- str_detect(rownames(counts), "-U$")
  usa_ambiguous <- str_detect(rownames(counts), "-A$")

  if(any(introns) & any(usa_introns)){
    stop("Gene ids with both -U and -I suffixes, can't determine read matrix type.")
  }
  if(any(introns)){
    spliced_genes <- !introns
  } else if(any(usa_introns)){
    spliced_genes <- !usa_introns & !usa_ambiguous
  } else {
    stop('No counts corresponding to intronic reads detected,
         must have "-I" or "-U" at the end of gene name to signify intronic reads.')
  }
  if(all(!spliced_genes)){
    stop("Missing spliced genes in counts matrix.")
  }


  if(which_counts == "spliced") {
    counts <- counts[spliced_genes | usa_ambiguous,]
  }
  # remove -I -U -A at the end of the gene names
  genes <- str_remove(rownames(counts), "-[IUA]$")

  # aggregate Matrix counts by gene name using matrix math
  # get unique gene names, encode genes as factors
  ugenes <- unique(genes)
  fgenes <- factor(genes, levels = ugenes)
  # create a transformer matrix, with 1s at positions that correspond to old/new names, zeros elsewhere
  transformer <- Matrix::sparseMatrix(as.integer(fgenes), 1:length(fgenes), x =1)
  rownames(transformer) <- ugenes
  # matrix multiplication to collapse
  counts <- transformer %*% counts
  # drop extra slots silently
  counts <- suppressWarnings(BiocGenerics::updateObject(counts))
  return(counts)

}
