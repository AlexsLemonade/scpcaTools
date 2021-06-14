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
read_alevin <- function(quant_dir, intron_mode = FALSE, usa_mode = FALSE,
                        which_counts = c("cDNA", "intron"), intron_metadata_path){

  which_counts <- match.arg(which_counts)

  if(usa_mode == TRUE) {
    # read in counts using read_usa mode
    sce <- read_usa_mode(quant_dir, which_counts)

  } else {
    # use tximport for all non-usa mode
    txi <- tximport::tximport(file.path(quant_dir, "alevin", "quants_mat.gz"), type = "alevin")
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = txi$counts))

    # collapse intron counts for intron_mode = TRUE
    if (intron_mode == TRUE) {
      sce < - collapse_intron_counts(counts(sce), which_counts, intron_metadata_path)
    }
  }
  return(sce)
}

read_usa_mode <- function(quant_dir, which_counts = c("cDNA", "intron")){

  which_counts <- match.arg(which_counts)

  # read in .mtx files
  mtx <- Matrix::readMM(file = file.path(quant_dir, "alevin", "quants_mat.mtx"))%>%
    Matrix::t() %>%
    as("dgCMatrix")
  cols <- readr::read_csv(file = file.path(quant_dir, "alevin", "quants_mat_cols.txt"),
                          col_names = "gene_id")
  rows <- readr::read_csv(file = file.path(quant_dir, "alevin", "quants_mat_rows.txt"),
                          col_names = "barcodes")
  dimnames(mtx) <- list(cols$gene_id, rows$barcodes)

  if (which_counts == "cDNA") {
    # only combine counts from S + A
    intron_genes <- rownames(mtx)[grep("-U", rownames(mtx))]
    splice_mtx <- mtx[!(rownames(mtx) %in% intron_genes),]
    # remove A from rowname
    rownames(splice_mtx) <- gsub("-A", "", rownames(splice_mtx))
    # combine counts based on gene name
    counts <- Matrix.utils::aggregate.Matrix(splice_mtx, row.names(splice_mtx))

  } else if (which_counts == "intron") {
    # combine counts from U, S, and A
    # remove A from rowname
    rownames(mtx) <- gsub("-A", "", rownames(mtx))
    rownames(mtx) <- gsub("-U", "", rownames(mtx))
    # combine counts based on gene name
    counts <- Matrix.utils::aggregate.Matrix(mtx, row.names(mtx))

  }
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  return(sce)
}

