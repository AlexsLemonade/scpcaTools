dir <- system.file("extdata", package="scpcaData")
usa_dir <- file.path(dir, "alevin_usa_Breast_Cancer_3p_LT")
sce_size <- c(60321, 9319)


test_that("Read USA mode works", {
  spliced <- read_alevin(usa_dir,
                         usa_mode = TRUE,
                         which_counts = "spliced")
  expect_s4_class(spliced, "SingleCellExperiment")
  expect_equal(dim(spliced), sce_size)
})
