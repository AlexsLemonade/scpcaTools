set.seed(1665)
sce <- sim_sce(n_cells = 200, n_genes = 200, n_empty = 100000)

test_that("Cell filtering with emptyDropsCellRanger() is as expected", {
  filtered_sce <- filter_counts(sce)
  expect_lt(ncol(filtered_sce), ncol(sce))
  expect_true(all(colnames(filtered_sce) %in% colnames(sce)))
  expect_equal(nrow(filtered_sce), nrow(sce))
  expect_true(metadata(filtered_sce)$filtering_method == "emptyDropsCellRanger")
})

test_that("Cell filtering with emptyDrops() is as expected", {
  filtered_sce <- filter_counts(sce, cr_like = FALSE)
  expect_lt(ncol(filtered_sce), ncol(sce))
  expect_true(all(colnames(filtered_sce) %in% colnames(sce)))
  expect_equal(nrow(filtered_sce), nrow(sce))
  expect_true(metadata(filtered_sce)$filtering_method == "emptyDrops")
})

test_that("Cell filtering with UMI cutoff is as expected", {
  low_cells_sce <- sim_sce(n_cells = 100, n_genes = 200, n_empty = 10000)
  low_cells_filtered_sce <- filter_counts(low_cells_sce)
  expect_lt(ncol(low_cells_filtered_sce), ncol(low_cells_sce))
  expect_true(all(colnames(low_cells_filtered_sce) %in% colnames(low_cells_sce)))
  expect_equal(nrow(low_cells_filtered_sce), nrow(low_cells_sce))
  expect_true(metadata(low_cells_filtered_sce)$filtering_method == "UMI cutoff")
})

test_that("Cell filtering removes row stats", {
  sce <- scuttle::addPerFeatureQCMetrics(sce)
  filtered_sce <- filter_counts(sce)
  expect_null(rowData(filtered_sce)$mean)
  expect_null(rowData(filtered_sce)$detected)
})
