test_that("simulated SCE objects are created", {
  sce <- sim_sce(n_genes = 20, n_cells = 10, n_empty = 100)
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(dim(sce), c(20, 110))
  sce <- sim_sce(n_genes = 20, n_cells = 10, n_empty = 0)
  expect_equal(dim(sce), c(20, 10))
})
