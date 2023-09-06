# quiet messages that will come with zellkonverter
suppressPackageStartupMessages(library(zellkonverter))

set.seed(1665)
sce <- sim_sce(n_cells = 100, n_genes = 200, n_empty = 0)
colData(sce) <- DataFrame("test_column" = sample(0:10, 100, rep = TRUE))
rowData(sce) <- DataFrame("test_row" = sample(0:10, 200, rep = TRUE))

# define anndata output
anndata_file <- tempfile(fileext = ".h5")

test_that("Conversion of SCE to AnnData works as expected", {

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

  # test that H5 file is created with new assay name
  # add logcounts
  logcounts(sce) <- counts(sce)
  new_anndata_file <- tempfile(fileext = ".h5")
  expect_snapshot_file({
    sce_to_anndata(sce, new_anndata_file, x_assay_name = "logcounts")
    new_anndata_file
  })

})

test_that("Conversion of SCE to AnnData fails as expected", {

  # check that inputting an improper file name causes a failure
  anndata_bad_file <- tempfile(fileext = ".rds")
  expect_error(sce_to_anndata(sce, anndata_bad_file))

  # improper assay name causes a failure
  expect_error(sce_to_anndata(sce, anndata_file, x_assay_name = "not an assay"))
})
