# establish single test sce
sce <- sim_sce(n_genes = 20, n_cells = 100, n_empty = 0)
metadata(sce)$sample_id <- "sample_id"
metadata(sce)$library_id <- "library_id"

# create sample data frame
sample_metadata_df <- data.frame(
  sample_id = "sample_id",
  library_id = "library_id"
)

# add sample metadata to object
sce <- add_sample_metadata(sce,
                           metadata_df = sample_metadata_df)

# add batch column to sce object
sce$library_id <- "library_id"

# create a merged object for testing the function still works with a merged object
set.seed(1665)
total_cells <- 24
total_genes <- 12
# create three individual objects
sce1 <- sim_sce(n_cells = total_cells / 3, n_genes = total_genes, n_empty = 0)
sce2 <- sim_sce(n_cells = total_cells / 3, n_genes = total_genes, n_empty = 0)
sce3 <- sim_sce(n_cells = total_cells / 3, n_genes = total_genes, n_empty = 0)

# combine into a list
sce_list <- purrr::imap(
  list(
    "sce1" = sce1,
    "sce2" = sce2,
    "sce3" = sce3
  ),
  # add some colData and metadata to each SCE object
  \(sce, batch){

    # add some colData so merge_sce_list doesn't give a warning
    colData(sce)[["sum"]] <- runif(ncol(sce))

    # add some metadata
    library_id <- batch
    sample_id <- batch

    metadata(sce)$library_id <- library_id
    metadata(sce)$sample_id <- sample_id
    metadata(sce)$sample_metadata <- data.frame(
      sample_id = sample_id,
      library_id = library_id
    )

    return(sce)
  }
)

# create merged sce
merged_sce <- merge_sce_list(sce_list = sce_list,
                             retain_coldata_cols = "sum")


test_that("`metadata_to_coldata` works as expected with a single object", {

  sce_with_metadata <- metadata_to_coldata(sce,
                                           batch_column = "library_id")

  # check that colData now contains sample id column
  expect_contains(
    colnames(colData(sce_with_metadata)),
    "sample_id"
  )

  # check that sample id column contains expected value
  expect_equal(
    unique(colData(sce_with_metadata)$sample_id),
    "sample_id"
  )

})

test_that("`metadata_to_coldata` works as expected with a merged object", {


  merged_sce_with_metadata <- metadata_to_coldata(merged_sce,
                                                  batch_column = "library_id")

  # check that colData now contains sample id column
  expect_contains(
    colnames(colData(merged_sce_with_metadata)),
    "sample_id"
  )

  # check that sample id column contains expected value
  expect_setequal(
    unique(colData(merged_sce_with_metadata)$sample_id),
    c("sce1", "sce2", "sce3")
  )

})

test_that("`metadata_to_coldata` fails as expected", {

  # sce is required
  expect_error(metadata_to_coldata(sce = "not an sce"))

  # column name is missing from colData
  expect_error(metadata_to_coldata(sce,
                                   batch_column = "not a column"))

  # column is missing from sample metadata
  # add a new column that is only in the colData and not in the sample metadata
  sce$batch_column <- "library_id"
  expect_error(metadata_to_coldata(sce,
                                   batch_column = "batch_column"))

  # no sample metadata present
  metadata(sce) <- metadata(sce)[!names(metadata(sce)) %in% "sample_metadata"]
  expect_error(metadata_to_coldata(sce,
                                   batch_column = "library_id"))

})

test_that("`metadata_to_coldata` fails with mis-matching library ids", {

  # replace existing sample metadata with one that has a library id
  # missing from the sce$library_id column
  sample_metadata_df <- data.frame(
    sample_id = "sample_id",
    library_id = "not a matching library"
  )
  metadata(sce)$sample_metadata <- sample_metadata_df

  expect_error(metadata_to_coldata(sce,
                                   batch_column = "library_id"))

})
