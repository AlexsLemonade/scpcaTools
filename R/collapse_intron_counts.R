#' Merge counts from intron reads with corresponding cDNA reads
#'
#' @param counts Counts matrix in sparse matrix format.
#' @param which_counts The type of counts (cDNA or intron) to include if alignment to intronic regions is TRUE. Default is FALSE.
#' @param intron_metadata Full path to a two column tsv file containing gene names for both spliced and intronic regions.
#'
#' @return SingleCellExperiment of unfiltered gene x cell counts matrix
#' @export
#'
#' @examples
#' \dontrun{
#' collapse_intron_counts(counts, which_counts = "cDNA",
#' intron_metadata_path = "<path to intron metadata.tsv>")
#' }
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
