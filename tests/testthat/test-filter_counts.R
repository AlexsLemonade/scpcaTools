set.seed(1665)
sce <- sim_sce(n_cells = 100, n_genes = 200, n_empty = 1000)

test_that("Cell filtering is as expected", {
  filtered_sce <- filter_counts(sce)
  expect_equal(ncol(filtered_sce), 100)
  expect_equal(colnames(filtered_sce), colnames(sce)[1:100])
  expect_equal(nrow(filtered_sce), 200)
})

test_that("Cell filtering removes row stats", {
  sce <- scuttle::addPerFeatureQCMetrics(sce)
  filtered_sce <- filter_counts(sce)
  expect_null(rowData(filtered_sce)$mean)
  expect_null(rowData(filtered_sce)$detected)
})
