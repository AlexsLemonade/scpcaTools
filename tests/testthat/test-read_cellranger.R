dir <- system.file("extdata", package = "scpcaData")
cellranger_dir <- file.path(dir, "Breast_Cancer_3p_LT", "cellranger_cdna", "outs")
h5_file <- file.path(cellranger_dir, "filtered_feature_bc_matrix.h5")

# expected matrix size
sce_size <- c(19970, 687)

test_that("Check that cellranger import works", {
  sce <- read_cellranger(h5_file)
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(dim(sce), sce_size)
  # check that column names are barcodes
  col_barcode <- str_detect(colnames(sce), "^[ACGT]+$")
  expect_true(all(col_barcode))
})
