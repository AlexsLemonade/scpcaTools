# create a SCE with random data
set.seed(1665)
sce <- sim_sce(n_cells = 100, n_genes = 200, n_empty = 0)

test_that("Check QC addition", {
  mito <- c(rownames(sce)[1:10], "ZZZZZZ") # include a gene that is not in the sce
  sce <- add_cell_mito_qc(sce, mito = mito, miQC = TRUE)
  expected_cols <- c("sum", "detected", "total",
                     "subsets_mito_sum","subsets_mito_detected", "subsets_mito_percent",
                     "prob_compromised")
  # check column names
  expect_true(all(expected_cols %in% names(colData(sce))))
  # make sure we don't get all zeros, which would indicate match failure
  expect_false(is.null(metadata(sce)$miQC_model))
  expect_gt(mean(sce$subsets_mito_percent), 0)

})
