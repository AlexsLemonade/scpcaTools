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
collapse_intron_counts <- function(counts, which_counts = c("cDNA", "intron"), intron_metadata_path){

  which_counts <- match.arg(which_counts)
  intron_metadata <- readr::read_tsv(intron_metadata_path)

  if(which_counts == "intron") {
    shared_genes <- intersect(row.names(counts), rownames(intron_metadata))
     # replace row names with -I appended with corresponding spliced gene
    row.names(counts)[which(row.names(counts) %in% shared_genes)] <- intron_metadata[shared_genes, "spliced"]
    # aggregate Matrix counts by gene name
    counts <- Matrix.utils::aggregate.Matrix(counts, row.names(counts))

  } else if (which_counts == "cDNA") {
    rows_to_keep <- intersect(row.names(counts), intron_metadata$spliced)
    counts <- counts[rows_to_keep,]
  }
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  return(sce)

}
