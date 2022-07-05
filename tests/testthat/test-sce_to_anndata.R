set.seed(1665)
sce <- sim_sce(n_cells = 100, n_genes = 200, n_empty = 0)

test_that("Conversion of SCE to AnnData works as expected", {
  anndata_file <- "test_anndata.h5"

  # test that the H5 file is created
  expect_snapshot(sce_to_anndata(sce, anndata_file))

  # check that inputting an improper file name causes a failure
  anndata_bad_file <- "test_anndata"
  expect_error(sce_to_anndata(sce, anndata_bad_file))
})
