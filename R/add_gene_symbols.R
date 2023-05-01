#' Add Gene symbols to a SingleCellExperiment from a GTF file
#'
#' @param sce SingleCellExperiment object
#' @param gene_info  A GRanges object or data frame with annotation information corresponding to the SCE.
#'   Most importantly should contain a `gene_id` column that matches the row names of the SCE
#'   and a `gene_name` column with the gene symbols should be included.
#'
#' @return SingleCellExperiment with gene_symbols added (or replaced) to the rowData table.
#'
#' @import SummarizedExperiment
#' @importFrom rlang .data
#'
#' @export
#'

add_gene_symbols <- function(sce, gene_info) {
  if (!is(sce, "SingleCellExperiment")) {
    stop("sce must be a SingleCellExperiment object.")
  }
  if (is.null(gene_info$gene_id)) {
    stop("`gene_info` must contain a `gene_id` column.")
  }
  if (is.null(gene_info$gene_name)) {
    stop("`gene_info` must contain a `gene_name` column.")
  }

  gene_symbols <- gene_info |>
    as.data.frame() |>
    dplyr::select(
      "gene_id" = "gene_id",
      "gene_symbol" = "gene_name"
    ) |>
    dplyr::filter(.data$gene_symbol != "NA") |>
    tidyr::drop_na("gene_symbol") |>
    dplyr::distinct() |>
    ## in case there are any duplicate gene_ids (there shouldn't be!)
    dplyr::group_by(.data$gene_id) |>
    dplyr::summarise(gene_symbol = paste(.data$gene_symbol, collapse = ";")) |>
    dplyr::pull("gene_symbol", name = .data$gene_id)

  rowData(sce)$gene_symbol <- unname(gene_symbols[rownames(sce)])
  return(sce)
}
