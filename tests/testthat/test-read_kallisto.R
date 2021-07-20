dir <- system.file("extdata", package="scpcaData")
kallisto_dir <- file.path(dir, "Breast_Cancer_3p_LT/kallisto_txome")

# expected matrix size
sce_size <- c(60666, 30339)

test_that("Check that kallisto import works", {
  sce <- read_kallisto(kallisto_dir)
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(dim(sce), sce_size)
  # check that column names are barcodes
  col_barcode <- str_detect(colnames(sce), "^[ACGT]+$")
  expect_true(all(col_barcode))
})
