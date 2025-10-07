# quiet messages that will come with zellkonverter
suppressPackageStartupMessages(library(zellkonverter))

set.seed(1665)
sce <- sim_sce(n_cells = 100, n_genes = 200, n_empty = 0)
# need to pull out barcodes to add to colData, otherwise colnames get replaced with NULL
barcodes <- colnames(sce)
colData(sce) <- DataFrame(
  "test_column" = sample(0:10, 100, rep = TRUE),
  "na_column" = NA,
  "na_char_column" = NA_character_,
  "some_na" = c("a", NA),
  row.names = barcodes
)
rowData(sce) <- DataFrame("test_row" = sample(0:10, 200, rep = TRUE))

# add some metadata
metadata(sce) <- list(
  library_id = "library1",
  sample_metadata = data.frame(
    library_id = "library1",
    sample_id = "sample1"
  ),
  metadata_dataframe = tibble::tibble(
    "char_col" = c("name1", NA),
    "num_col" = 1:2,
    "na_col" = c(NA, NA),
    "list_col" = list(list(a = 1, b = "foo"), list(a = 2, b = "bar"))
  ),
  metadata_S4_DataFrame = DataFrame(head(iris)),
  # this should get removed:
  metadata_list = list(value = "value")
)

# define anndata output
tempdir <- tempdir()
anndata_file <- file.path(tempdir, "anndata.h5")

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


  # expect that sample metadata has been removed from converted SCE
  expect_setequal(
    c("library_id", "sample_metadata", "metadata_dataframe", "metadata_S4_DataFrame"),
    names(metadata(converted_sce))
  )

  # check that list columns aren't present
  expect_setequal(
    names(metadata(converted_sce)$metadata_dataframe),
    c("char_col", "num_col", "na_col")
  )

  # check that the NA-only column remains logical
  expect_type(
    metadata(converted_sce)$metadata_dataframe$na_col,
    "logical"
  )


  # test that H5 file is created with new assay name
  # add logcounts
  logcounts(sce) <- counts(sce)
  new_anndata_file <- file.path(tempdir, "new_anndata.h5")
  expect_snapshot_file({
    sce_to_anndata(sce, new_anndata_file, x_assay_name = "logcounts")
    new_anndata_file
  })

  new_sce <- zellkonverter::readH5AD(new_anndata_file, X_name = "X")
  # check that counts is in assay names and no logcounts
  expect_setequal(
    assayNames(new_sce),
    c("X", "counts")
  )

  # check that X is equal to logcounts
  expect_equal(
    logcounts(sce),
    assay(new_sce, "X")
  )
})

test_that("Conversion of SCE to AnnData works with additional arguments", {
  # test that the H5 file is created with additional options
  expect_snapshot_file({
    sce_to_anndata(sce, anndata_file, verbose = FALSE)
    anndata_file
  })
})

test_that("Conversion of SCE to AnnData fails as expected", {
  # check that inputting an improper file name causes a failure
  anndata_bad_file <- tempfile(fileext = ".rds")
  expect_error(sce_to_anndata(sce, anndata_bad_file))

  # improper assay name causes a failure
  expect_error(sce_to_anndata(sce, anndata_file, x_assay_name = "not an assay"))

  # conversion fails if < 2 cells
  small_sce <- sce[, 1]
  expect_error(sce_to_anndata(small_sce, anndata_file))

  # check that conversion fails with wrong compression type
  expect_error(sce_to_anndata(sce, anndata_file, compression = "not a compression"))
})
