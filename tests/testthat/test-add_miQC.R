# create a SCE with random data
set.seed(1665)
sce <- sim_sce(n_cells = 100, n_genes = 200, n_empty = 0)

test_that("miQC addition is correct", {
  mito <- rownames(sce)[1:10]
  sce <- scuttle::addPerCellQCMetrics(sce, subsets = list(mito = mito))
  sce <- add_miQC(sce)
  expect_true("prob_compromised" %in% colnames(colData(sce)))
  expect_true("miQC_model" %in% names(metadata(sce)))
})

test_that("miQC failures are handled", {
  sce$subsets_mito_percent <- 0
  sce <- add_miQC(sce)
  expect_true(all(is.na(sce$prob_compromised)))
  expect_null(metadata(sce)$miQC_model)
})
