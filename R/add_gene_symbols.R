#' Add Gene symbols to a SingleCellExperiment from a GTF file
#'
#' @param sce SingleCellExperiment object
#' @param gtf  A GRanges object or data frame with annotation information corresponding to the SCE.
#'   Most importantly, the `gene_id` column should match the row names of the SCE
#'   and a `gene_name` column with the gene symbols should be included.
#'
#' @return SingleCellExperiment with gene_symbols added (or replaced) to the rowData table.
#'
#' @import SummarizedExperiment
#' @importFrom rlang .data

add_gene_symbols <- function(sce, gtf){
  gene_symbols <- gtf |>
    as.data.frame() |>
    dplyr::select("gene_id" = "gene_id",
                  "gene_symbol" = "gene_name") |>
    dplyr::filter(.data$gene_symbol != "NA") |>
    tidyr::drop_na(.data$gene_symbol) |>
    dplyr::distinct() |>
    ## in case there are any duplicate gene_ids (there shouldn't be!)
    dplyr::group_by(.data$gene_id) |>
    dplyr::summarise(gene_symbol = paste(.data$gene_symbol, collapse = ";")) |>
    dplyr::pull(.data$gene_symbol, name = .data$gene_id)

  rowData(sce)$gene_symbol = unname(gene_symbols[rownames(sce)])
  return(sce)
}
