# generate testing data
set.seed(1665)
merged_sce <- sim_sce(n_cells = 303, n_genes = 100, n_empty = 0)
batches <- rep(c("a", "b", "c"), each = 101)
barcodes <- rownames(colData(merged_sce))

# set up colData rownames as batch-barcode
new_rownames <- glue::glue("{batches}-{barcodes}")
rownames(colData(merged_sce)) <- new_rownames

# Add "sample" to colData
colData(merged_sce)$sample <- batches

# add PCs for testing (exact numbers don't matter)
# make a 300x100 matrix
reducedDim(merged_sce, "PCA") <- matrix(runif(303 * 100, min = 0, max = 100), nrow = 303)

# grab pcs and batches
pcs <- reducedDim(merged_sce, "PCA")
batches <- merged_sce$sample

test_that("`filter_pcs` works as expected", {

  # expect rownames of returned pcs to be the same as the batch column
  named_pcs <- filter_pcs(pcs, batches)
  expect_true(all.equal(rownames(named_pcs), batches))

})

test_that("`filter_pcs` removes NAs as expected", {

  batches[1:5] <- NA # set a subset of batch labels to NAs

  # expect rownames to be equal to only the batches without NA
  named_pcs <- filter_pcs(pcs, batches)
  expect_true(all.equal(rownames(named_pcs), batches[!is.na(batches)]))
  expect_false(length(rownames(named_pcs)) == length(batches))
})

test_that("`filter_pcs` fails without batch labels", {

  all_na_batches <- NA
  expect_error(filter_pcs(pcs, all_na_batches))

})

test_that("`downsample_pcs` works as expected", {

  test_frac = 0.8
  downsampled <- downsample_pcs(pcs, frac_cells = test_frac)
  downsampled_pcs <-downsampled
  downsampled_batch_labels <- rownames(downsampled)

  # check that the correct number of cells were downsampled
  expect_true(length(rownames(downsampled_pcs)) == round(test_frac*length(rownames(pcs))))
  expect_true(length(downsampled_batch_labels) == round(test_frac*length(batches)))

  # check that the columns in downsampled were in the original pcs
  expect_true(length(colnames(downsampled_pcs)) == length(colnames(pcs)))

})

test_that("`downsample_pcs` fails as expected", {

  # check that downsampling fails when providing a fraction outside of 0-1 range
  expect_error(downsample_pcs(pcs, frac_cells = 2))
  expect_error(downsample_pcs(pcs, frac_cells = 0))

})

test_that("`calculate_silhouette_width` works as expected", {

  nreps <- 5
  asw <- calculate_silhouette_width(integrated_sce = merged_sce,
                                    pc_name = "PCA",
                                    nreps = 5,
                                    batch_column = "sample")

  expected_cols <- c(
    "rep", "silhouette_width", "silhouette_cluster",
    "other_cluster", "pc_name"
  )
  # check column names
  expect_true(all(expected_cols %in% colnames(asw)))

  # check that pc name is PCA
  expect_true(all(asw$pc_name == "PCA"))

  # check that number of reps is correct
  expect_true(all(sort(unique(asw$rep)) == 1:nreps))

  # check that width is in expected range
  expect_true(all(asw$silhouette_width >= -1 |
                    asw$silhouette_width <= 1))

})

test_that("`calculate_silhouette_width`fails as expected", {

  # missing pc name
  expect_error(calculate_silhouette_width(merged_sce,
                                          pc_name = "test_PC",
                                          batch_column = "sample"))

  # error in nreps
  expect_error(calculate_silhouette_width(merged_sce,
                                          pc_name = "PCA",
                                          nreps = "twenty",
                                          batch_column = "sample"))

  # incorrect batch labels
  expect_error(calculate_silhouette_width(merged_sce,
                                          pc_name = "PCA",
                                          batch_column = "batch"))

})
