test_that("adding gene_symbols to sce works", {
  sce <- sim_sce(n_cells = 5, n_genes = 5, n_empty = 0)
  gene_info <- data.frame(gene_id = rownames(sce)[3:5],
                          gene_name = c("a", "b", "c"))
  sce <- add_gene_symbols(sce, gene_info)
  expect_equal(rowData(sce)$gene_symbol, c(NA, NA, "a", "b", "c"))
})
