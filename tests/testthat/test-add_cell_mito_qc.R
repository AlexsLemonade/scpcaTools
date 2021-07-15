# create a SCE with random data
set.seed(1665)
ncells <- 100
# generate barcodes
cell_barcodes <- replicate(
  ncells + 10, # account for (rare) potential duplicates
  paste0(sample(c("A","T","G","C"), 10, replace = TRUE), collapse = '')
) %>%
  unique() %>%
  head(ncells)
genes <- sprintf("GENE%04d", 1:200)

counts <- matrix(rpois(ncells * length(genes), 2),
                 ncol = ncells,
                 dimnames = list(genes, cell_barcodes))

sce <- SingleCellExperiment(list(counts = counts))

test_that("Check QC addition", {
  mito <- c(genes[1:10], "GENE9999") # include a gene that is not in the sce
  sce <- add_cell_mito_qc(sce, mito = mito)
  expected_cols <- c("sum", "detected", "total",
                     "mito_sum","mito_detected", "mito_percent")
  expect_true(all(expected_cols %in% names(colData(sce))))
  expect_gt(mean(sce$mito_percent), 0)

})
