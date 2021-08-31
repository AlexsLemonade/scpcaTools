test_that("merging alternative experiments works as expected", {
  # create two SCEs with random data
  set.seed(1665)
  sce <- sim_sce(n_cells = 100, n_genes = 200, n_empty = 0)
  alt_sce <- sim_sce(n_cells = 40, n_genes = 10, n_empty = 0)
  # rename the alt_exp so it shares half of its cells with the base
  colnames(alt_sce) <- c(colnames(alt_sce)[1:20], colnames(sce)[1:20])
  # merge alt into main
  alt_name = "alt"
  merged_sce <- merge_altexp(sce, alt_sce, alt_name)
  # test general expectations
  expect_true(alt_name %in% altExpNames(merged_sce))
  expect_equal(colnames(sce), colnames(altExp(merged_sce, alt_name)))

  # the cells not in the original alt_exp should all be zero counts
  missing_cells <- which(!colnames(sce) %in% colnames(alt_sce))
  missing_counts <- counts(altExp(merged_sce, alt_name)[, missing_cells])
  expect_true(all(missing_counts == 0))
})
