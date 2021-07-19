dir <- system.file("extdata", package="scpcaData")
alevin_dir <- file.path(dir, "Breast_Cancer_3p_LT/alevin_txome")
usa_dir <- file.path(dir, "Breast_Cancer_3p_LT/alevinfry_usa_filtered")
intron_dir <- file.path(dir, "Breast_Cancer_3p_LT/alevinfry_intron_filtered")

# expected alevin matrix size
sce_alevin_size <- c(60275, 10378)
# expected alevin-fry matrix size
sce_af_size <- c(60321, 614)

test_that("reading salmon alevin data works", {
  sce <- read_alevin(alevin_dir)
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(dim(sce), sce_alevin_size)
  # check that column names are barcodes
  col_barcode <- str_detect(colnames(sce), "^[ACGT]+$")
  expect_true(all(col_barcode))
})

test_that("reading alevin-fry USA mode works", {
  sce <- read_alevin(usa_dir,
                     usa_mode = TRUE,
                     which_counts = "spliced")
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(dim(sce), sce_af_size)
  # check that column names are barcodes
  col_barcode <- str_detect(colnames(sce), "^[ACGT]+$")
  expect_true(all(col_barcode))
  # no remaining unspliced
  unmerged_genes <- str_subset(rownames(sce), "-[IUA]$")
  expect_length(unmerged_genes, 0)
})


test_that("reading alevin-fry intron mode works", {
  sce <- read_alevin(intron_dir,
                     intron_mode = TRUE,
                     which_counts = "unspliced")
  expect_s4_class(sce, "SingleCellExperiment")
  expect_equal(dim(sce), sce_af_size)
  # check that column names are barcodes
  col_barcode <- stringr::str_detect(colnames(sce), "^[ACGT]+$")
  expect_true(all(col_barcode))
  # no remaining unspliced
  unmerged_genes <- str_subset(rownames(sce), "-[IUA]$")
  expect_length(unmerged_genes, 0)

})
