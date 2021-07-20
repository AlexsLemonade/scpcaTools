dir <- system.file("extdata", package="scpcaData")
usa_dir <- file.path(dir, "Breast_Cancer_3p_LT/alevinfry_usa_filtered")
intron_dir <- file.path(dir, "Breast_Cancer_3p_LT/alevinfry_intron_filtered")

# expected matrix size
sce_size <- c(60321, 614)

test_that("Check that reading USA mode works", {
  sce <- read_alevin(usa_dir,
                     usa_mode = TRUE,
                    which_counts = "spliced")
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(dim(sce), sce_size)
  # check that column names are barcodes
  col_barcode <- str_detect(colnames(sce), "^[ACGT]+$")
  expect_true(all(col_barcode))
  unmerged_genes <- str_subset(rownames(sce), "-[IUA]$")
  expect_length(unmerged_genes, 0)
})


test_that("Check that reading intron mode works", {
  sce <- read_alevin(intron_dir,
                     intron_mode = TRUE,
                     which_counts = "unspliced")
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(dim(sce), sce_size)
  # check that column names are barcodes
  col_barcode <- stringr::str_detect(colnames(sce), "^[ACGT]+$")
  expect_true(all(col_barcode))
  unmerged_genes <- str_subset(rownames(sce), "-[IUA]$")
  expect_length(unmerged_genes, 0)

})
