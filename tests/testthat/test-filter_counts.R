set.seed(1665)
sce <- sim_sce(n_cells = 1000, n_genes = 200, n_empty = 100)

test_that("Cell filtering with emptyDropsCellRanger() is as expected", {
  filtered_sce <- filter_counts(sce)
  expect_lt(ncol(filtered_sce), ncol(sce))
  expect_true(all(colnames(filtered_sce) %in% colnames(sce)))
  expect_equal(nrow(filtered_sce), nrow(sce))
})

test_that("Cell filtering with emptyDrops() is as expected", {
  filtered_sce <- filter_counts(sce, cr_like = FALSE)
  expect_lt(ncol(filtered_sce), ncol(sce))
  expect_true(all(colnames(filtered_sce) %in% colnames(sce)))
  expect_equal(nrow(filtered_sce), nrow(sce))
})

test_that("Cell filtering removes row stats", {
  sce <- scuttle::addPerFeatureQCMetrics(sce)
  filtered_sce <- filter_counts(sce)
  expect_null(rowData(filtered_sce)$mean)
  expect_null(rowData(filtered_sce)$detected)
})
