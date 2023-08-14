# establish test sce
sce <- sim_sce(n_genes = 20, n_cells = 100, n_empty = 0)

# create sample data frame
sample_metadata_df <- data.frame(
  sample_id = "sample_id",
  library_id = "library_id"
)

test_that("`add_sample_metadata` works as expected", {

  updated_sce <- add_sample_metadata(sce,
                                     metadata_df = sample_metadata_df)

  expect_equal(
    metadata(updated_sce)$sample_metadata,
    sample_metadata_df
  )

})

test_that("`add_sample_metadata` fails as exepected", {

  expect_error(add_sample_metadata(sce = "not an sce",
                                   metadata_df = sample_metadata_df))

  expect_error(add_sample_metadata(sce,
                                   metadata_df = "not a data frame"))

})
