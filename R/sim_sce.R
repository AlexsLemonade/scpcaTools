#' Create a random SingleCellExperiment object
#'
#' @param n_genes The number of genes to simulate
#' @param n_cells The number of cells to create.
#' @param n_empty How many "empty" droplets to create.
#'   This will generally be much larger than the number of cells.
#' @param n_groups How many expression groups (cell clusters, roughly) to create.
#'
#' @return A SingleCellExperiment object.
#'
#' @importFrom stats rexp rpois
#' @importFrom utils head
#'
#' @examples
#' \dontrun{
#' sim_sce(n_genes = 100, n_cells = 100, n_empty = 2000, n_groups = 4)
#' }
#'
sim_sce <- function(n_genes = 200, n_cells = 100, n_empty = 1000, n_groups = 3) {
  if (n_genes < 1) {
    stop("n_genes must be a positive number.")
  }
  if (n_cells < 1) {
    stop("n_cells must be a postive number.")
  }
  # enforce some minimums
  n_empty <- max(n_empty, 0)
  n_groups <- max(n_groups, 1)

  if (n_cells + n_empty > 1000000) {
    warning("Did you really want to simulate more than a million droplets?")
  }

  # generate a random distribution of gene expression for each group
  gene_expr <- rexp(n_genes * n_groups, 0.2)
  # expression for each gene x cell
  cell_expr <- matrix(rpois(n_genes * n_cells, lambda = gene_expr), nrow = n_genes)

  # empty drops will have expression calculated from the overall mean
  empty_expr <- rowMeans(cell_expr) * 0.01
  empty_expr <- matrix(rpois(n_genes * n_empty, lambda = empty_expr), nrow = n_genes)

  # generate gene names
  genes <- sprintf("GENE%04d", seq_len(n_genes))

  # random cell barcodes
  cell_barcodes <- replicate(
    (n_cells + n_empty) * 1.1, # account for (rare, but possible) duplicates
    paste0(sample(c("A", "T", "G", "C"), 12, replace = TRUE), collapse = "")
  ) |>
    unique() |>
    head(n_cells + n_empty)

  # combine matrices and build SCE
  all_expr <- as(cbind(cell_expr, empty_expr), "CsparseMatrix")
  dimnames(all_expr) <- list(genes, cell_barcodes)
  SingleCellExperiment(list(counts = all_expr))
}
