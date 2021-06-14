#' Title
#'
#' @param counts
#' @param which_counts
#' @param intron_metadata
#'
#' @return
#' @export
#'
#' @examples
collapse_intron_counts <- function(counts, which_counts = c("splice", "intron"), intron_metadata){

  which_counts <- match.arg(which_counts)

  if(which_counts == "intron") {
    shared_genes <- intersect(row.names(counts), rownames(intron_metadata))
     # replace row names with -I appended with corresponding spliced gene
    row.names(counts)[which(row.names(counts) %in% shared_genes)] <- intron_metadata[shared_genes, "spliced"]
    # aggregate Matrix counts by gene name
    counts <- Matrix.utils::aggregate.Matrix(counts, row.names(counts))

  } else if (which_counts == "splice") {
    intron_rows <- counts[grep("-I", row.names(counts)),]
    counts <- counts[-intron_rows,]
  }
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  return(sce)

}
