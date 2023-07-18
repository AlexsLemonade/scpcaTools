# generate testing data
set.seed(1665)
merged_sce <- sim_sce(n_cells = 300, n_genes = 100, n_empty = 0)
batches <- rep(c("a", "b", "c"), each = 100)
barcodes <- rownames(colData(merged_sce))

# set up colData rownames as batch-barcode
new_rownames <- glue::glue("{batches}-{barcodes}")
rownames(colData(merged_sce)) <- new_rownames

# Add "sample" to colData
colData(merged_sce)$sample <- batches

# add PCs for testing (exact numbers don't matter)
# make a 300x100 matrix
reducedDim(merged_sce, "PCA") <- matrix(runif(300 * 100, min = 0, max = 100), nrow = 300)

# grab pcs and batches
pcs <- reducedDim(merged_sce, "PCA")
batches <- merged_sce$sample

test_that("`set_pc_rownames` works as expected", {

  # expect rownames of returned pcs to be the same as the batch column
  named_pcs <- set_pc_rownames(pcs, batches)
  expect_true(all.equal(rownames(named_pcs), batches))

})

test_that("`set_pc_rownames` removes NAs as expected", {

  batches[1:5] <- NA # set a subset of batch labels to NAs

  # expect rownames to be equal to only the batches without NA
  named_pcs <- set_pc_rownames(pcs, batches)
  expect_true(all.equal(rownames(named_pcs), batches[!is.na(batches)]))
  expect_false(length(rownames(named_pcs)) == length(batches))
})

test_that("`set_pc_rownames` fails without batch labels", {

  short_batches <- c("a", "b", "c")
  expect_error(scpcaTools:::set_pc_rownames(pcs, short_batches))

})

test_that("`downsample_pcs` works as expected", {

  downsampled <- downsample_pcs(pcs, frac_cells = 0.8)
  downsampled_pcs <- downsampled$pcs
  downsampled_batch_labels <- downsampled$batch_labels

  # check that the correct number of cells were downsampled
  expect_true(length(rownames(downsampled_pcs)) == 0.8*length(rownames(pcs)))
  expect_true(length(downsampled_batch_labels) == 0.8*length(batches))

  # check that the columns in downsampled were in the original pcs
  expect_true(length(colnames(downsampled_pcs)) == length(colnames(pcs)))

})

test_that("`downsample_pcs` fails as expected", {

  # check that downsampling fails when providing a fraction outside of 0-1 range
  expect_error(downsample_pcs(pcs, frac_cells = 2))
  expect_error(downsample_pcs(pcs, frac_cells = 0))

})

test_that("`calculate_silhouette_width` works as expected", {

  asw <- calculate_silhouette_width(merged_sce,
                                    pc_name = "PCA",
                                    batch_column = "samples")

  expected_cols <- c(
    "rep", "silhouette_width", "silhouette_cluster",
    "other_cluster", "pc_name"
  )
  # check column names
  expect_true(all(expected_cols %in% colnames(asw)))

  # check that pc name is PCA
  expect_true(asw$pc_name == "PCA")

  # check that nreps equal 1-20
  expect_true(asw$nreps %in% 1:20)

})

test_that("`calculate_silhouette_width`fails as expected", {

  # missing pc name
  expect_error(calculate_silhouette_width(merged_sce,
                                          pc_name = "test_PC",
                                          batch_column = "samples"))

  # error in nreps
  expect_error(calculate_silhouette_width(merged_sce,
                                          pc_name = "PCA",
                                          nreps = "twenty",
                                          batch_column = "samples"))

  # incorrect batch labels
  expect_error(calculate_silhouette_width(merged_sce,
                                          pc_name = "PCA",
                                          batch_column = "batch"))

})
