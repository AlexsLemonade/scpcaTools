#' Merge counts from intron reads with corresponding cDNA reads
#'
#' @param counts Counts matrix with rownames corresponding to gene names and colnames corresponding to cell barcodes.
#' @param which_counts If intron_mode is TRUE, which type of counts should be included:
#'   only counts aligned to spliced cDNA ("spliced") or all spliced and unspliced cDNA ("unspliced").
#'   Default is "spliced".
#'
#' @return unfiltered gene x cell counts matrix
#'
#' @examples
#' \dontrun{
#' collapse_intron_counts(counts, which_counts = "spliced")
#' }
#'
#' @noRd
collapse_intron_counts <- function(counts, which_counts = c("spliced", "unspliced")){

  which_counts <- match.arg(which_counts)

  intron_genes <- str_subset(rownames(counts), "-I$")
  if(is.na(intron_genes)){
    stop('No counts corresponding to intronic reads detected,
         must have tag "-I" at the end of gene name to signify intronic read.')
  }

  if(which_counts == "spliced") {
    spliced_genes <- str_subset(rownames(counts), "-I$", negate = TRUE)
    counts <- counts[spliced_genes,]

  } else if (which_counts == "unspliced") {
    # remove -I at the end of the gene name
    rownames(counts) <- str_remove(rownames(counts), "-I$")
    # aggregate Matrix counts by gene name
    counts <- Matrix.utils::aggregate.Matrix(counts, rownames(counts))
  }
  return(counts)

}
