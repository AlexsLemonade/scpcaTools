#' Title
#'
#' @param quant_dir
#' @param tool
#' @param intron_mode
#' @param usa_mode
#'
#' @return
#' @export
#'
#' @examples
read_alevin <- function(quant_dir, intron_mode = FALSE, usa_mode = FALSE, which_counts, intron_metadata){

  if(usa_mode == TRUE) {
    # if usa mode is true then read counts using read_usa
    sce <- read_usa_mode(quant_dir, which_counts)

  } else {
    sce <- tximport::tximport(file.path(quant_dir, "alevin", "quants_mat.gz")) %>%
      SingleCellExperiment::SingleCellExperiment(list(counts = counts))

    if (intron_mode == TRUE) {
      sce < - collapse_intron_counts(counts(sce), which_counts, intron_metadata)
    }
  }
  return(sce)
}

read_usa_mode <- function(quant_dir, which_counts){

  # read in .mtx files
  mtx <- Matrix::readMM(file = file.path(quant_dir, "alevin", "quants_mat.mtx"))
  cols <- readr::read_csv(file = file.path(quant_dir, "alevin", "quants_mat_cols.txt"),
                          col_names = "gene_id")
  rows <- readr::read_csv(file = file.path(quant_dir, "alevin", "quants_mat_rows.txt"),
                          col_names = "barcodes")

  if (which_counts = c("S")) {
    # only combine counts from S + A

  } else if ("U" %in% which_counts) {
    # combine counts from U, S, and A

  }
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  return(sce)
}

collapse_intron_counts <- function(counts, which_counts, intron_metadata){

  if(which_counts === c("U", "S")) {
    shared_genes <- intersect(row.names(counts), rownames(intron_metadata))
    # replace row names with -I appended with corresponding spliced gene
    row.names(counts)[which(row.names(counts) %in% shared_genes)] <- intron_metadata[shared_genes, "spliced"]
    # aggregate Matrix counts by gene name
    counts <- Matrix.utils::aggregate.Matrix(counts, row.names(counts))

  } else if (which_counts == c("S")) {
    intron_rows <- grep("-I", row.names(counts))
    counts <- counts[-intron_rows,]
  }
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  return(sce)

}
