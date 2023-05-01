set.seed(1665)
sce <- sim_sce(n_cells = 100, n_genes = 200, n_empty = 0)
colData(sce) <- DataFrame("test_column" = sample(0:10, 100, rep= TRUE))
rowData(sce) <- DataFrame("test_row" = sample(0:10, 200, rep = TRUE))

test_that("Conversion of SCE to AnnData works as expected", {
  anndata_file <- "test_anndata.h5"
  # quiet messages that will come with zellkonverter
  suppressPackageStartupMessages(library(zellkonverter))

  # test that the H5 file is created
  expect_snapshot_file({
    sce_to_anndata(sce, anndata_file)
    anndata_file
  })

  # some tests that the converted object contains col/rowData found in original SCE
  converted_sce <- zellkonverter::readH5AD(anndata_file)
  expect_equal(dim(sce), dim(converted_sce))
  expect_equal(colnames(colData(sce)), colnames(colData(converted_sce)))
  expect_equal(colnames(rowData(sce)), colnames(rowData(converted_sce)))

  # check that inputting an improper file name causes a failure
  anndata_bad_file <- "test_anndata"
  expect_error(sce_to_anndata(sce, anndata_bad_file))
})
