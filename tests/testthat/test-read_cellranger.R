dir <- system.file("extdata", package="scpcaData")
cellranger_dir <- file.path(dir, "Breast_Cancer_3p_LT/cellranger_cdna")

# expected matrix size
sce_size <- c(19970, 687)

test_that("Check that cellranger import works", {
  sce <- read_cellranger(cellranger_dir)
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(dim(sce), sce_size)
  # check that column names are barcodes
  col_barcode <- stringr::str_detect(colnames(sce), "^[ACGT]+$")
  expect_true(all(col_barcode))
})
