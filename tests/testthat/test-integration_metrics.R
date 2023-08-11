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
# add barcode to colData
colData(merged_sce)$cell_id <- barcodes

# add PCs for testing (exact numbers don't matter)
# make a 300x100 matrix
reducedDim(merged_sce, "PCA") <- matrix(runif(303 * 100, min = 0, max = 100), nrow = 303)
rownames(reducedDim(merged_sce, "PCA")) <- colnames(merged_sce)

# grab pcs and batches
pcs <- reducedDim(merged_sce, "PCA")
batches <- merged_sce$sample

# create list of sce objects
batch_ids <- unique(batches)
sce_list <- purrr::map(batch_ids,
                       \(batch){
                         individual_sce <- merged_sce[, which(colData(merged_sce)$sample == batch)]
                         # individual sce object columns should be named without batches
                         colnames(individual_sce) <- stringr::word(colnames(individual_sce), -1, sep = "-")
                         return(individual_sce)
                       }) |>
  purrr::set_names(batch_ids)

test_that("`filter_pcs` works as expected", {
  # expect rownames of returned pcs to be the same as the batch column
  named_pcs <- filter_pcs(pcs, batches)
  expect_true(all.equal(rownames(named_pcs), batches))
})

test_that("`filter_pcs` works as expected when rename_pcs is FALSE", {
  # set input PCA names to be something else
  named_pcs <- filter_pcs(pcs, batches, rename_pcs = FALSE)
  expect_true(all.equal(rownames(named_pcs), rownames(pcs)))
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
  test_frac <- 0.8
  downsampled <- downsample_pcs(pcs, frac_cells = test_frac)
  downsampled_pcs <- downsampled
  downsampled_batch_labels <- rownames(downsampled)

  # check that the correct number of cells were downsampled
  expect_true(length(rownames(downsampled_pcs)) == round(test_frac * length(rownames(pcs))))
  expect_true(length(downsampled_batch_labels) == round(test_frac * length(batches)))

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
  asw <- calculate_silhouette_width(
    merged_sce = merged_sce,
    pc_names = "PCA",
    nreps = 5,
    batch_column = "sample"
  )

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
    pc_names = "test_PC",
    batch_column = "sample"
  ))

  # error in nreps
  expect_error(calculate_silhouette_width(merged_sce,
    pc_names = "PCA",
    nreps = "twenty",
    batch_column = "sample"
  ))

  # incorrect batch labels
  expect_error(calculate_silhouette_width(merged_sce,
    pc_names = "PCA",
    batch_column = "batch"
  ))
})

test_that("`within_batch_ari_from_pcs` works as expected", {

  ari_from_pcs <- within_batch_ari_from_pcs(individual_sce_list = sce_list,
                                            merged_sce = merged_sce,
                                            pc_name = "PCA",
                                            batch_column = "sample",
                                            cell_id_column = "cell_id")

  expected_cols <- c(
    "ari", "batch_id", "pc_name"
  )
  # check column names
  expect_true(all(expected_cols %in% colnames(ari_from_pcs)))

  # check that pc name is PCA
  expect_true(all(ari_from_pcs$pc_name == "PCA"))

  # check that ari is in expected range
  expect_true(all(ari_from_pcs$ari >= 0 |
                    ari_from_pcs$ari <= 1))
})

test_that("`within_batch_ari_from_pcs`fails as expected", {

  # missing pc name
  expect_error(within_batch_ari_from_pcs(individual_sce_list = sce_list,
                                         merged_sce = merged_sce,
                                         pc_name = "test_PC",
                                         batch_column = "sample",
                                         cell_id_column = "cell_id"))


  # incorrect batch labels
  expect_error(within_batch_ari_from_pcs(individual_sce_list = sce_list,
                                         merged_sce = merged_sce,
                                         pc_name = "PCA",
                                         batch_column = "batch",
                                         cell_id_column = "cell_id"))
  # incorrect barcode label
  expect_error(within_batch_ari_from_pcs(individual_sce_list = sce_list,
                                         merged_sce = merged_sce,
                                         pc_name = "PCA",
                                         batch_column = "sample",
                                         cell_id_column = "not a barcode"))

  ## missing names for sce list
  # save to a new variable so we can use the sce list again later
  unnamed_sce_list <- sce_list
  names(unnamed_sce_list) <- NULL

  expect_error(within_batch_ari_from_pcs(individual_sce_list = unnamed_sce_list,
                                         merged_sce = merged_sce,
                                         pc_name = "PCA",
                                         batch_column = "sample",
                                         cell_id_column = "cell_id"))

})

test_that("`calculate_within_batch_ari` works as expected", {

  # add fastmnn pca to test multiple pcs
  reducedDim(merged_sce, "fastMNN_PCA") <- matrix(runif(303 * 100, min = 0, max = 100), nrow = 303)

  ari <- calculate_within_batch_ari(individual_sce_list = sce_list,
                                    merged_sce = merged_sce,
                                    pc_names = c("PCA", "fastMNN_PCA"),
                                    batch_column = "sample",
                                    cell_id_column = "cell_id")

  expected_cols <- c(
    "ari", "batch_id", "pc_name"
  )
  # check column names
  expect_true(all(expected_cols %in% colnames(ari)))

  # check that pc name is PCA
  expect_true(all(ari$pc_name %in% c("PCA", "fastMNN_PCA")))

  # check that ari is in expected range
  expect_true(all(ari$ari >= 0 |
                    ari$ari <= 1))
})

test_that("`calculate_within_batch_ari`fails as expected", {

  # missing pc name
  expect_error(calculate_within_batch_ari(individual_sce_list = sce_list,
                                          merged_sce = merged_sce,
                                          pc_names = "test_PC",
                                          batch_column = "sample",
                                          cell_id_column = "cell_id"))


  # incorrect batch labels
  expect_error(within_batch_ari_from_pcs(individual_sce_list = sce_list,
                                         merged_sce = merged_sce,
                                         pc_names = c("PCA", "fastMNN_PCA"),
                                         batch_column = "batch",
                                         cell_id_column = "cell_id"))

  # incorrect barcode labels
  expect_error(within_batch_ari_from_pcs(individual_sce_list = sce_list,
                                         merged_sce = merged_sce,
                                         pc_names = c("PCA", "fastMNN_PCA"),
                                         batch_column = "sample",
                                         cell_id_column = "not a barcode"))

})

test_that("`calculate_ilisi` works as expected", {

  ilisi <- calculate_ilisi(merged_sce = merged_sce,
                           pc_name = "PCA",
                           batch_column = "sample")

  expected_cols <- c(
    "ilisi_score", "cell_barcode", "batch_id", "ilisi_score_norm"
  )
  # check column names
  expect_true(all(expected_cols %in% colnames(ilisi)))

  # check that normalized ilisi scores are in in expected range
  expect_true(
    all(
      dplyr::between(ilisi$ilisi_score_norm, 0, 1)
    )
  )
})

test_that("`calculate_ilisi` fails as expected", {

  # fails with missing PC name
  expect_error(calculate_ilisi(merged_sce = merged_sce,
                               pc_name = "test_PC",
                               batch_column = "sample"))

  # fails with incorrect batch column
  expect_error(calculate_ilisi(merged_sce = merged_sce,
                               pc_name = "PCA",
                               batch_column = "library_id"))

})
