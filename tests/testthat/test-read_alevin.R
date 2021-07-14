dir <- system.file("extdata", package="scpcaData")
usa_dir <- file.path(dir, "Breast_Cancer_3p_LT/alevinfry_usa_filtered")
intron_dir <- file.path(dir, "Breast_Cancer_3p_LT/alevinfry_intron_filtered")

# expected matrix size
sce_size <- c(60321, 614)

test_that("Check that reading USA mode works", {
  spliced <- read_alevin(usa_dir,
                         usa_mode = TRUE,
                         which_counts = "spliced")
  expect_s4_class(spliced, "SingleCellExperiment")
  expect_equal(dim(spliced), sce_size)
})


test_that("Check that reading intron mode works", {
  spliced <- read_alevin(intron_dir,
                         intron_mode = TRUE,
                         which_counts = "unspliced")
  expect_s4_class(spliced, "SingleCellExperiment")
  expect_equal(dim(spliced), sce_size)
})
